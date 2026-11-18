"""Bempp on the unit cube, using the three mirror planes of the body, as the
MFS does (cube_mfs.py). The answer to "what would a boundary element code
gain from the symmetry?", measured instead of estimated.

The mesh (graded towards the edges, like bem_cube_graded.py) is built to be
invariant under the eight reflections. A reflection acts on the RWG
coefficients as a signed permutation T_R. A current of parity class c lies in
the span of the symmetrised functions B_c e_k = sum_R sgn_c(R) T_R e_k, one
for each orbit of degrees of freedom, and the Galerkin equations need to be
tested only with one function per orbit. So only the rows of the EFIE matrix
that belong to one octant are assembled (an eighth of the work), folded with
B_c, and each class is solved by a small LU factorisation. The result is the
Galerkin solution of the full mesh, to rounding.

    python3 bem_cube_sym.py 8,12,16 2 [reg sing] [--wave=60,25,35] [--check]

--check also solves the full system and prints the difference.
"""
import itertools
import sys
import time

import bem_env  # noqa: F401  (memory for the OpenCL driver; must come before bempp)
import bempp_cl.api as bempp
import numpy as np
from numba import njit, prange
import scipy.linalg as sla
import scipy.sparse as sp

import bem_reference as br
from cube_mfs import ROT, diagonal_orbits

K = 2 * np.pi
H = 0.5
D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])
REFL = np.array(list(itertools.product((1.0, -1.0), repeat=3)))
CLASSES = [tuple(float(v) for v in c) for c in itertools.product((1, -1), repeat=3)]


def sign(cls, R):
    return float(np.prod([s for s, r in zip(cls, R) if r < 0]))


def sym_grid(n, g):
    """Cube mesh with n intervals per edge, graded towards the edges with
    exponent g, invariant under the eight reflections. Elements are ordered
    by reflection, so element e of the octant (+,+,+) has the image
    e + r * T under reflection number r (T elements per octant)."""
    m = n // 2
    s = 0.5 * (1 - (1 - np.arange(m + 1) / m) ** g)
    pts, tris = [], []
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        base = len(pts)
        for i in range(m + 1):
            for j in range(m + 1):
                x = np.zeros(3)
                x[a], x[b], x[c] = H, s[i], s[j]
                pts.append(x)
        idx = lambda i, j: base + i * (m + 1) + j
        for i in range(m):
            for j in range(m):
                tris.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)))
                tris.append((idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)))
    pts, tris = np.array(pts), np.array(tris)
    key, verts, elems = {}, [], []
    for R in REFL:
        loc = []
        for x in pts * R:
            kk = tuple(np.round(x, 12) + 0.0)
            if kk not in key:
                key[kk] = len(verts)
                verts.append(x)
            loc.append(key[kk])
        loc = np.array(loc)
        t = loc[tris]
        elems.append(t if np.prod(R) > 0 else t[:, ::-1])
    verts = np.array(verts)
    vmap = np.array([[key[tuple(np.round(x * R, 12) + 0.0)] for x in verts] for R in REFL])
    return bempp.Grid(verts.T, np.concatenate(elems).T.astype(np.uint32)), vmap, len(tris)


def reflections(space, vmap, T):
    """Signed permutations of the RWG degrees of freedom under the eight
    reflections: perm[r, j], sgn[r, j] with T_R e_j = sgn * e_perm."""
    grid = space.grid
    ne, N = grid.number_of_elements, space.global_dof_count
    el, vt = grid.elements.T, grid.vertices.T
    cen = vt[el].mean(axis=1)
    opp = np.zeros((ne, 3), np.int64)
    sig = np.zeros((ne, 3))
    lc = np.array([[1.0 / 3], [1.0 / 3]])
    for e in range(ne):                                    # which vertex each local function flows away from, and its sign
        w = space.evaluate(e, lc)[:, :, 0].T               # (3 functions, 3 components)
        dvec = cen[e][None, :] - vt[el[e]]                 # (3 vertices, 3)
        cs = (w @ dvec.T) / (np.linalg.norm(w, axis=1)[:, None] * np.linalg.norm(dvec, axis=1)[None, :])
        v = np.argmax(np.abs(cs), axis=1)
        opp[e], sig[e] = el[e][v], np.sign(cs[np.arange(3), v])
    l2g = space.local2global
    first = np.full(N, -1)
    for e in range(ne - 1, -1, -1):
        first[l2g[e]] = e
    rsign = REFL
    ridx = {tuple(r): i for i, r in enumerate(rsign)}
    perm = np.zeros((8, N), np.int64)
    sgn = np.zeros((8, N))
    j = np.arange(N)
    e = first
    l = np.argmax(l2g[e] == j[:, None], axis=1)
    r_of = e // T                                          # reflection number of the octant the element lies in
    for q in range(8):
        r2 = np.array([ridx[tuple(rsign[r] * rsign[q])] for r in range(8)])[r_of]
        e2 = e + (r2 - r_of) * T
        p2 = vmap[q][opp[e, l]]
        l2 = np.argmax(opp[e2] == p2[:, None], axis=1)
        assert np.all(opp[e2, l2] == p2)
        perm[q], sgn[q] = l2g[e2, l2], sig[e, l] * sig[e2, l2]
    return perm, sgn


@bempp.callable(complex=True, parameterized=True)
def class_trace(x, n, domain_index, result, parameters):
    """(part of p exp(i k d.x) with parity (sx, sy, sz)) x n; parameters = (k, d, p, s)."""
    e0 = e1 = e2 = 0j
    for a in (1.0, -1.0):
        for b in (1.0, -1.0):
            for c in (1.0, -1.0):
                sg = 1.0 + 0j
                if a < 0:
                    sg = sg * parameters[7]
                if b < 0:
                    sg = sg * parameters[8]
                if c < 0:
                    sg = sg * parameters[9]
                ph = sg * np.exp(1j * parameters[0] * (a * x[0] * parameters[1] + b * x[1] * parameters[2] + c * x[2] * parameters[3])) / 8
                e0 += a * parameters[4] * ph
                e1 += b * parameters[5] * ph
                e2 += c * parameters[6] * ph
    result[0] = e1 * n[2] - e2 * n[1]
    result[1] = e2 * n[0] - e0 * n[2]
    result[2] = e0 * n[1] - e1 * n[0]


@njit(parallel=True, cache=True)
def fold(A, rows, cols, coef):
    """out[i, k] = sum over the eight images r of coef[r, k] * A[rows[i], cols[r, k]]."""
    out = np.zeros((len(rows), cols.shape[1]), np.complex128)
    for i in prange(len(rows)):
        a = A[rows[i]]
        for r in range(cols.shape[0]):
            for k in range(cols.shape[1]):
                out[i, k] += coef[r, k] * a[cols[r, k]]
    return out


def classes_of(d, p):
    x = np.random.default_rng(0).uniform(-0.5, 0.5, (40, 3))
    out = []
    for c in CLASSES:
        v = sum(sign(c, R) * (p * R)[None, :] * np.exp(1j * K * ((x * R) @ d))[:, None] for R in REFL)
        if np.abs(v).max() > 1e-10:
            out.append(c)
    return out


def solve(n, g, points, d=D, p=P, rotate=True):
    """Field at points and the times (set-up, assembly, fold and LU, evaluation)."""
    t0 = time.time()
    grid, vmap, T = sym_grid(n, g)
    div = bempp.function_space(grid, "RWG", 0)
    N = div.global_dof_count
    perm, sgn = reflections(div, vmap, T)
    rep = np.nonzero(perm.min(axis=0) == np.arange(N))[0]            # one degree of freedom per orbit
    l2g = div.local2global
    in_rep = np.zeros(N, bool)
    in_rep[rep] = True
    support = np.nonzero(in_rep[l2g].any(axis=1))[0].astype(np.uint32)
    test = bempp.function_space(grid, "SNC", 0, support_elements=support, include_boundary_dofs=False)
    full_of = np.full(test.global_dof_count, -1)                    # test function -> degree of freedom of the full space
    tl = test.local2global[support]
    ok = test.local_multipliers[support] != 0                      # local functions that are not degrees of freedom have multiplier 0
    full_of[tl[ok]] = l2g[support][ok]
    rows = np.nonzero(in_rep[full_of])[0]
    order = np.argsort(full_of[rows])
    rows = rows[order]
    assert np.array_equal(full_of[rows], rep)
    t1 = time.time()
    efie = bempp.operators.boundary.maxwell.electric_field(div, div, test, K)
    A = np.asarray(efie.weak_form().to_dense())                     # (test functions, N)
    t2 = time.time()
    # The rotation about the body diagonal permutes the mirrors (and maps the mesh to itself), so the part of
    # parity c turned j times is the solution of class c for the wave turned back j times, turned forward:
    # at most four of the eight classes are folded and factorised (cube_mfs.diagonal_orbits; rotate=False: all eight).
    orbit = diagonal_orbits(classes_of(d, p)) if rotate else {c: [0] for c in classes_of(d, p)}
    J = np.zeros((3, N), complex)
    sizes = []
    for c, turns in orbit.items():
        sg = np.array([sign(c, R) for R in REFL])
        B = sp.coo_matrix(((sg[:, None] * sgn[:, rep]).ravel(), (perm[:, rep].ravel(), np.tile(np.arange(len(rep)), 8))),
                          shape=(N, len(rep))).tocsc()
        keep = np.nonzero(np.abs(B).sum(axis=0).A1 > 1e-12)[0]      # orbits on a mirror plane drop out of some classes
        B = B[:, keep]
        cols = rep[keep]
        Ac = fold(A, rows[keep], perm[:, cols], sg[:, None] * sgn[:, cols])   # A[rows[keep]] @ B
        lu = sla.lu_factor(Ac, overwrite_a=True, check_finite=False)
        for j in turns:
            Q = np.linalg.matrix_power(ROT, j)
            par = np.concatenate(([K], Q.T @ d, Q.T @ p, c))
            b = bempp.GridFunction(div, fun=class_trace, dual_space=test, function_parameters=par).projections(test)[rows][keep]
            J[j] += B @ sla.lu_solve(lu, b, check_finite=False)
            sizes.append(len(keep))
    t3 = time.time()
    E = 0
    for j in range(3):
        if np.any(J[j]):
            Q = np.linalg.matrix_power(ROT, j)
            lam = bempp.GridFunction(div, coefficients=J[j])
            E = E - (bempp.operators.potential.maxwell.electric_field(div, np.ascontiguousarray((points @ Q).T), K) * lam).T @ Q.T
    t4 = time.time()
    return E, N, sizes, (t1 - t0, t2 - t1, t3 - t2, t4 - t3), grid


if __name__ == "__main__":
    wave = [a[7:] for a in sys.argv if a.startswith("--wave=")]
    args = [a for a in sys.argv if not a.startswith("--")]
    ns = [int(a) for a in args[1].split(",")]
    g = float(args[2])
    reg = int(args[3]) if len(args) > 3 else 4
    sing = int(args[4]) if len(args) > 4 else 6
    if wave:
        from cube_mfs import wave as wave_dp
        d, p = wave_dp(*(float(v) for v in wave[0].split(",")))
    else:
        d, p = D, P
    bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
    bempp.GLOBAL_PARAMETERS.quadrature.singular = sing
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    ref = np.load("cube_oblique.npz")["Eref"] if wave else np.load("fig_cost_cube.npz")["ref_E"]
    solve(4, 1.0, pts, d, p)                                        # warm-up: compilation
    for n in ns:
        E, N, sizes, (ts, ta, tl, te), grid = solve(n, g, pts, d, p)
        err = np.linalg.norm(E - ref) / np.linalg.norm(ref)
        msg = (f"cube, symmetric BEM n={n} g={g:g} quad={reg}/{sing} dofs={N} class sizes={sizes} setup={ts:.1f}s assembly={ta:.1f}s "
               f"fold+lu={tl:.1f}s eval={te:.1f}s total={ts + ta + tl + te:.1f}s diff from the MFS reference: {err:.2e}")
        if "--check" in sys.argv:
            t0 = time.time()
            lam, div, _ = br.efie_solve(grid, K, d, p)
            Ef = -(bempp.operators.potential.maxwell.electric_field(div, pts.T, K) * lam).T
            msg += f"  | full solve {time.time() - t0:.1f}s, difference of the fields {np.linalg.norm(E - Ef) / np.linalg.norm(Ef):.1e}"
        print(msg, flush=True)
        tag = f"n{n}_g{g:g}" + ("" if (reg, sing) == (4, 6) else f"_q{reg}_{sing}") + ("_obl" if wave else "")
        np.savez(f"bem_cube_sym_{tag}.npz", pts=pts, E=E, k=K, n=n, g=g, ndof=N, d=d, p=p, t_setup=ts, t_assembly=ta,
                 t_lu=tl, t_eval=te, t_solve=ts + ta + tl + te)
