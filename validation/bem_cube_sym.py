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

import bem_mirror as bm
import bem_reference as br
from cube_mfs import ROT, diagonal_orbits

K = 2 * np.pi
H = 0.5
D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])
REFL = bm.group()


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
    return bm.mirror_grid(np.array(pts), np.array(tris), REFL)


def solve(n, g, points, d=D, p=P, rotate=True):
    """Field at points and the times (set-up, assembly, fold and LU, evaluation)."""
    t0 = time.time()
    sym = bm.Mirror(*sym_grid(n, g), REFL, K)
    t2 = time.time()
    # The rotation about the body diagonal permutes the mirrors (and maps the mesh to itself), so the part of
    # parity c turned j times is the solution of class c for the wave turned back j times, turned forward:
    # at most four of the eight classes are folded and factorised (cube_mfs.diagonal_orbits; rotate=False: all eight).
    cl = bm.classes_of(d, p, K, REFL)
    orbit = diagonal_orbits(cl) if rotate else {c: [0] for c in cl}
    J, sizes = {}, []
    for j in range(3):
        only = [c for c in orbit if j in orbit[c]]
        if only:
            Q = np.linalg.matrix_power(ROT, j)
            J[j], sz = sym.solve(Q.T @ d, Q.T @ p, only)
            sizes += sz
    t3 = time.time()
    E = 0
    for j in J:
        Q = np.linalg.matrix_power(ROT, j)
        E = E + sym.field(J[j], points @ Q) @ Q.T
    t4 = time.time()
    return E, sym.N, sizes, (t2 - t0 - sym.t_assembly, sym.t_assembly, t3 - t2, t4 - t3), sym.div.grid


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
