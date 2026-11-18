"""Unit cube, any incident plane wave: MFS with edge sources and singular
tails, solved by the eight parity classes of the three mirror planes.

The cube and the source layout are symmetric under x -> -x, y -> -y and
z -> -z. Any incident field splits into eight parts of definite parity
(sx, sy, sz), E(R x) = s R E(x) for each mirror R; each part is scattered into
a field of the same parity and is represented by the sources of one octant:
the three edges that meet the vertex (+, +, +) with half of their modes (the
mirror along an edge maps the edge to itself and selects the parity of the
density), one vertex, and an eighth of the base points. Collocation points
are needed in one octant only. The class matrices are signed sums of the same
eight reflected blocks, so the line integrals are shared between the classes
(line_tree.py writes them straight into the class matrices) and each class
has its own small QR factorisation, which then serves every incident wave.

A wave travelling along a coordinate axis and polarised along another one
excites two classes; a general wave excites all eight.

    python3 cube_mfs.py J [theta_deg phi_deg psi_deg]

runs the standard wave of the other scripts (along -z, E along y) or the wave
arriving from (theta, phi) with polarisation angle psi, prints the time split,
the boundary residual and the optical-theorem mismatch, and compares with the
stored reference fields. For a general wave six of the eight classes are
images of two others under the rotation about a body diagonal, so four
matrices suffice (diagonal_orbits; test_cube_rotation.py).
"""
import itertools
import sys
import time

import numpy as np

import line_tree as lt
import pec_cube_lightning3d as m
from lsq import factor

COL = dict(col_p=3, col_q=0.5, beta=0.5)
REFL = [np.array(r, float) for r in itertools.product((1, -1), repeat=3)]
CLASSES = [tuple(float(v) for v in c) for c in itertools.product((1, -1), repeat=3)]


def balanced(J):
    """All resolution parameters grow with the depth J: modes per line
    (n0 + 2 log(1/d)), base grid, and the collocation that samples it. The
    tail exponents are the solver defaults (derived, DERIVATIONS.md B); the
    vertex tail (one exponent, DERIVATIONS.md D) is on."""
    n0 = max(4, 2 * J - 4)
    n_base, hmax = 8, 0.1          # one rule for all levels (NOTES.md item 68). Until 17 Nov 2026: (4, 0.25) for J <= 5 and
    #                                (6, 0.15) above, which gave the same errors 20 to 80% faster but failed from J = 12 on.
    return dict(q=0.5, gamma=2.0, kinds="em", n0=n0, n_base=n_base, vtail_n=9), hmax


def wave(theta, phi, psi):
    """Propagation direction d and polarisation p of a plane wave arriving from
    the direction (theta, phi) (degrees); psi rotates p from theta_hat towards phi_hat."""
    t, f, s = np.deg2rad([theta, phi, psi])
    rhat = np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
    that = np.array([np.cos(t) * np.cos(f), np.cos(t) * np.sin(f), -np.sin(t)])
    fhat = np.array([-np.sin(f), np.cos(f), 0.0])
    return -rhat, np.cos(s) * that + np.sin(s) * fhat


def sign(cls, R):
    return float(np.prod([s for s, r in zip(cls, R) if r < 0]))


def class_part(x, d, p, cls):
    """Part of the plane wave p exp(i k d.x) with parity cls, at the points x."""
    return sum(sign(cls, R) * (p * R)[None, :] * np.exp(1j * m.K * ((x * R) @ d))[:, None] for R in REFL) / 8


def classes_of(d, p):
    """The parity classes that the wave (d, p) excites."""
    x = np.random.default_rng(0).uniform(-0.5, 0.5, (40, 3))
    return [c for c in CLASSES if np.abs(class_part(x, d, p, c)).max() > 1e-12]


ROT = np.eye(3)[[1, 2, 0]]                             # a third of a turn about the body diagonal: (ROT v)_i = v_{i+1}


def diagonal_orbits(classes):
    """The cube is also symmetric under the rotation ROT about a body
    diagonal, which permutes the three mirrors. If E solves the problem of
    the wave (d, p), then x -> Q E(Q^T x) solves that of (Q d, Q p), and if E
    has the parities c, the turned field has the parities c[i+1]. So the
    part of parity c' = c turned j times is obtained from the matrix of
    class c with the wave turned back j times, and of the eight classes
    only four need a matrix: (+,+,+), (-,-,-) and one of each triple.
    Returns {class to assemble: [j, ...]}, with j = 0 for the class itself."""
    orbit, todo = {}, list(classes)
    while todo:
        c = todo.pop(0)
        orbit[c] = [0]
        for j in (1, 2):
            cj = tuple(c[(i + j) % 3] for i in range(3))
            if cj in todo:
                todo.remove(cj)
                orbit[c].append(j)
    return orbit


class Cube:
    """Assembles and factorises the given parity classes once; solve() then
    takes any plane wave whose classes are among them."""

    def __init__(self, J, kw=None, hmax=None, tail=9, classes=CLASSES, col=COL, lsq="auto"):
        if kw is None:
            kw, hmax = balanced(J)
        kw = dict(kw)
        if not tail:
            kw.pop("vtail_n", None)                    # "no tail" means neither the edge nor the vertex tail
        full = m.Sources(J, tail_n=tail, **kw)         # the full layout; only its bookkeeping is used
        octant = lambda y: y[np.all(y > 0, axis=1)] if len(y) else y
        self.pe, self.pm = octant(full.points), octant(full.mpoints)
        assert np.abs(full.points).min() > 1e-9, "a base point lies on a mirror plane"
        self.lines = [ln for ln in full.lines if ln.yb > 0 and ln.yc > 0
                      and not (isinstance(ln, m.VertexTail) and ln.pts[0, 2] < 0)]
        self.classes = [tuple(float(v) for v in c) for c in classes]
        self.depths = full.depths
        # column bookkeeping: unknowns are grouped in triples ("slots")
        self.slot0 = np.cumsum([0, len(self.pe), len(self.pm)] + [ln.n for ln in self.lines])
        self.keep = {c: self._keep(c) for c in self.classes}
        self.unknowns = {c: int(k.sum()) for c, k in self.keep.items()}
        # collocation in the octant, rows weighted by sqrt(area) * (distance to the nearest edge)^beta
        d_min = self.depths[-1] if J >= 0 else 0.2
        lv = max(2, int(np.ceil(np.log(d_min / (2 * m.H)) / np.log(col["col_q"]))) + 1)
        x, w, T, _, _ = m.cube_cells(lv, col["col_p"], col["col_q"], hmax)
        sel = np.all(x > -1e-12, axis=1)
        self.x, self.T = x[sel], T[sel]
        self.rw = np.sqrt(w[sel]) * m.edge_dist(self.x) ** col["beta"]
        self.rows = 2 * len(self.x)
        t0 = time.time()
        A = self.matrix(self.x, self.T * self.rw[:, None, None], self.classes)
        self.t_assemble = time.time() - t0
        t0 = time.time()
        self.fac = {}
        for i, c in enumerate(self.classes):
            self.fac[c] = factor(np.ascontiguousarray(A[i][:, :self.unknowns[c]]), lsq)
            if self.fac[c].dR.min() < 1e-12 * self.fac[c].dR.max():
                print(f"warning: class {c} is numerically rank deficient "
                      f"({self.fac[c].dR.min() / self.fac[c].dR.max():.1e})")
        self.t_factor = time.time() - t0

    def _keep(self, cls):
        """Which unknowns (slot, component) a class keeps."""
        keep = [np.ones((len(self.pe) + len(self.pm), 3), bool)]
        for ln in self.lines:
            if isinstance(ln, m.VertexTail):
                keep.append(np.ones((ln.n, 3), bool))
                continue
            mag = ln.kernel is m.magnetic_pairs
            par = (-1.0) ** getattr(ln, "deg", np.arange(ln.n))[:, None]
            sc = np.where(np.arange(3) == getattr(ln, "axis_index", ln.ax), -1.0, 1.0)[None, :] * (-1.0 if mag else 1.0)
            keep.append((cls[ln.ax] * par * sc > 0) & getattr(ln, "colmask", True))
        return np.concatenate(keep).reshape(-1)

    def matrix(self, x, T, classes):
        """Matrices of the given classes, (classes, rows, columns): rows are
        (target, component of T), columns the unknowns of the class (padded
        to the largest class). T is (M, nc, 3), or None for the three
        Cartesian components. The sum over the eight mirror images is taken
        inside the compiled kernels."""
        ncls = len(classes)
        nmax = max(self.unknowns[c] for c in classes)
        cmap = np.full((ncls, 3 * self.slot0[-1]), -1, np.int64)
        for i, c in enumerate(classes):
            cmap[i, self.keep[c]] = np.arange(self.unknowns[c])
        Tk = np.eye(3)[None] if T is None else np.ascontiguousarray(T)
        nc = Tk.shape[1]
        A = np.zeros((ncls, len(x) * nc, nmax), complex)
        sgn = np.array([[sign(c, R) for R in REFL] for c in classes])
        refl = np.array(REFL)
        x = np.ascontiguousarray(x)
        if len(self.pe) + len(self.pm):                # point dipoles: electric, then magnetic (all components kept)
            lt.points_add(x, Tk, refl, self.pe, self.pm if len(self.pm) else np.zeros((0, 3)), m.K, A,
                          np.ascontiguousarray(cmap[:, 0]), sgn)
        i = 0
        while i < len(self.lines):
            ln = self.lines[i]
            lo = 3 * self.slot0[2 + i]
            if isinstance(ln, m.VertexTail):           # a few tied point dipoles; all components kept
                for r, R in enumerate(REFL):
                    B = ln.block(x * R, np.broadcast_to(Tk * R, (len(x), nc, 3))).reshape(len(x) * nc, -1)
                    for j in range(ncls):
                        A[j][:, cmap[j, lo]:cmap[j, lo] + 3 * ln.n] += sgn[j, r] * B
                i += 1
                continue
            mate = self.lines[i + 1] if i + 1 < len(self.lines) else None
            pair = mate is not None and m.same_place(ln, mate)
            tail = isinstance(ln, m.TailLine)
            b = ln.base if tail else ln
            mag = b.kernel is m.magnetic_pairs
            cm_first = np.ascontiguousarray(cmap[:, lo:lo + 3 * ln.n])
            cm_mate = np.ascontiguousarray(cmap[:, lo + 3 * ln.n:lo + 3 * (ln.n + mate.n)]) if pair else None
            lt.tree_add(x, Tk, refl, b.ax, b.b, b.c, np.array(ln.pos) if tail else (b.yb, b.yc), b.lam, b.d, b.n, m.K, m.H,
                        A, None if mag else cm_first, cm_first if mag else cm_mate, sgn,
                        rot=ln.rot if tail else None,
                        we=ln.w if tail and not mag else None,
                        wm=(mate.w if pair else ln.w) if tail and (pair or mag) else None)
            i += 2 if pair else 1
        return A

    def solve(self, d, p, only=None):
        """Coefficients per class for the wave (d, p), and the relative
        boundary residual. `only`: solve just these classes (the others must
        then be handled by another Cube)."""
        coef, r2, b2, c2 = {}, 0.0, 0.0, 0.0
        for c in classes_of(d, p):
            if only is not None and c not in only:
                continue
            if c not in self.fac:
                raise ValueError(f"the wave excites class {c}, which was not assembled")
            b = -np.einsum("mab,mb->ma", self.T, class_part(self.x, d, p, c)).reshape(-1) * np.repeat(self.rw, 2)
            coef[c], res = self.fac[c].solve(b, residual=True)
            r2 += float(res) ** 2
            b2 += np.linalg.norm(b) ** 2
            c2 += np.linalg.norm(coef[c] * self.fac[c].scale) ** 2
        self.coef_norm, self.rhs_norm2 = float(np.sqrt(c2)), float(b2)
        return coef, np.sqrt(r2 / b2)

    def far_pattern(self, coef, d):
        """Far-field pattern in the direction d and on the quadrature sphere
        (for the cross-sections; patterns of class groups add)."""
        F = self.far(coef, np.vstack((d[None, :], _SPHERE_DIRS)))
        return F[0], F[1:]

    def field(self, coef, pts):
        cl = list(coef)
        A = self.matrix(pts, None, cl)
        return sum(A[i][:, :self.unknowns[c]] @ coef[c] for i, c in enumerate(cl)).reshape(-1, 3)

    def far(self, coef, dirs, rfar=1.0e6):
        return self.field(coef, rfar * dirs) * rfar * np.exp(-1j * m.K * rfar)

    def cross_sections(self, coef, d, p, nt=40, nphi=80):
        """Extinction (optical theorem) and scattering cross-sections."""
        g, w = np.polynomial.legendre.leggauss(nt)
        ph = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
        ct, phh = np.meshgrid(g, ph, indexing="ij")
        st = np.sqrt(1 - ct ** 2)
        dirs = np.column_stack(((st * np.cos(phh)).ravel(), (st * np.sin(phh)).ravel(), ct.ravel()))
        F = self.far(coef, np.vstack((dirs, d[None, :])))
        ext = 4 * np.pi / m.K * float(np.imag(np.conj(p) @ F[-1]))
        return ext, float(np.sum(np.repeat(w, nphi) * 2 * np.pi / nphi * np.sum(np.abs(F[:-1]) ** 2, axis=1)))


def dense(J):
    """The balanced family with two source depths per octave instead of one
    (ratio 2^(-1/2), 2 J + 1 depths down to the same deepest depth
    0.2 * 2^(-J)) and four collocation nodes per panel. The family with one
    depth per octave stalls near 1e-4 (NOTES.md, item 48); this one does not.
    Returns the keywords for Cube."""
    kw, hmax = balanced(J)
    kw.update(q=2 ** -0.5)
    return dict(kw=kw, hmax=hmax, col=dict(COL, col_p=4))


def tapered(J, n=None, col_p=4):
    """The balanced family with tapered depths, as in lightning solvers:
    n depths d_j = 0.2 exp(-sigma (sqrt(n) - sqrt(j))), j = n..1, dense near
    0.2 and sparse near the deepest depth 0.2 * 2^(-J) (sigma is fixed by
    the two ends). Default n = J + 1, the number of depths of the family
    with one depth per octave. Returns the keywords for Cube."""
    n = J + 1 if n is None else n
    sigma = J * np.log(2.0) / (np.sqrt(n) - 1.0)
    kw, hmax = balanced(J)
    kw.update(depths=0.2 * np.exp(-sigma * (np.sqrt(n) - np.sqrt(np.arange(n, 0, -1)))))
    return dict(kw=kw, hmax=hmax, col=dict(COL, col_p=col_p))


def run(J, d=m.D, p=m.P, pts=None, family="balanced", max_gb=6.0, **kw):
    """One wave: field at pts, total time, and a dictionary of diagnostics.
    family: "balanced" (one depth per octave) or "dense" (two). If the
    matrices of all classes together would exceed max_gb, the classes are
    assembled and solved in groups (the kernel work is then repeated)."""
    t0 = time.time()
    if family == "dense":
        kw = dict(dense(J), **kw)
        J = 2 * J
    elif family == "tapered":
        kw = dict(tapered(J, kw.pop("n", None), kw.pop("col_p", 4)), **kw)
        J = len(kw["kw"]["depths"]) - 1
    orbit = diagonal_orbits(classes_of(d, p))          # {class to assemble: the rotations that give the other classes}
    todo = list(orbit)
    E, info, r2, b2, c2 = None, None, 0.0, 0.0, 0.0
    ext = sca_F = 0.0
    t = 0.0
    dirs = np.vstack((d[None, :], _SPHERE_DIRS))
    while todo:
        cube = Cube(J, classes=todo[:1], **kw) if len(todo) == 1 else None
        if cube is None:
            n1 = _size(J, todo[0], kw)
            g = max(1, min(len(todo), int(max_gb * 2 ** 30 / (16.0 * n1[0] * n1[1]))))
            cube = Cube(J, classes=todo[:g], **kw)
        group, todo = cube.classes, todo[len(cube.classes):]
        sols = []
        for j in range(3):                             # the wave turned back by j thirds of a turn about the body diagonal
            only = [c for c in group if j in orbit[c]]
            if only:
                Q = np.linalg.matrix_power(ROT, j)
                coef, res = cube.solve(Q.T @ d, Q.T @ p, only=only)
                Eg = cube.field(coef, pts @ Q) @ Q.T if pts is not None else None
                E = Eg if E is None else E + Eg
                sols.append((Q, coef, res, cube.rhs_norm2, cube.coef_norm))
        t, t0 = t + time.time() - t0, None             # the cross-sections below are not part of the timing
        for Q, coef, res, b, cn in sols:
            r2, b2, c2 = r2 + res ** 2 * b, b2 + b, c2 + cn ** 2
            F = cube.far(coef, dirs @ Q) @ Q.T
            ext, sca_F = ext + F[0], sca_F + F[1:]
        part = dict(unknowns=sum(cube.unknowns.values()), per_class=list(cube.unknowns.values()), rows=cube.rows,
                    t_assemble=cube.t_assemble, t_factor=cube.t_factor)
        if info is None:
            info = part
        else:
            for key in ("unknowns", "t_assemble", "t_factor"):
                info[key] += part[key]
            info["per_class"] += part["per_class"]
        t0 = time.time()
    info.update(residual=float(np.sqrt(r2 / b2)), coef_norm=float(np.sqrt(c2)))
    info["ext"] = 4 * np.pi / m.K * float(np.imag(np.conj(p) @ ext))
    info["sca"] = float(np.sum(_SPHERE_W * np.sum(np.abs(sca_F) ** 2, axis=1)))
    info["optical"] = abs(info["ext"] - info["sca"]) / info["sca"]
    return E, t, info


def _size(J, cls, kw):
    """(rows, unknowns) of one class, without assembling anything."""
    kw = dict(kw)
    col = kw.pop("col", COL)
    src_kw, hmax, tail = kw.get("kw"), kw.get("hmax"), kw.get("tail", 9)
    if src_kw is None:
        src_kw, hmax = balanced(J)
    src_kw = dict(src_kw)
    if not tail:
        src_kw.pop("vtail_n", None)
    full = m.Sources(J, tail_n=tail, **src_kw)
    d_min = full.depths[-1] if J >= 0 else 0.2
    lv = max(2, int(np.ceil(np.log(d_min / (2 * m.H)) / np.log(col["col_q"]))) + 1)
    s, _ = m.graded_nodes_capped(lv, col["col_p"], col["col_q"], hmax)
    rows = 2 * 3 * (len(s) // 2) ** 2
    return rows, full.n_unknowns // 8 + 1               # an octant of the points, a quarter of the lines, half their modes


def _sphere_rule(nt=40, nphi=80):
    g, w = np.polynomial.legendre.leggauss(nt)
    ph = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    ct, phh = np.meshgrid(g, ph, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    dirs = np.column_stack(((st * np.cos(phh)).ravel(), (st * np.sin(phh)).ravel(), ct.ravel()))
    return dirs, np.repeat(w, nphi) * 2 * np.pi / nphi


_SPHERE_DIRS, _SPHERE_W = _sphere_rule()


if __name__ == "__main__":
    Js = [int(v) for v in sys.argv[1].split(",")]
    standard = len(sys.argv) < 5
    d, p = (m.D, m.P) if standard else wave(*(float(v) for v in sys.argv[2:5]))
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    run(2, d, p, pts)                                  # warm-up (compilation, threads)
    for J in Js:
        E, t, info = run(J, d, p, pts)
        print(f"cube J={J} d={np.round(d, 4)} p={np.round(p, 4)}: classes {len(info['per_class'])}, unknowns {info['per_class']}, "
              f"rows {info['rows']}, assembly {info['t_assemble']:.1f}s, factorisations {info['t_factor']:.1f}s, total {t:.1f}s")
        print(f"  boundary residual {info['residual']:.2e}, sigma_ext {info['ext']:.6f}, sigma_sca {info['sca']:.6f}, "
              f"optical mismatch {info['optical']:.1e}, |c| {info['coef_norm']:.1e}", flush=True)
        if standard:
            ref = np.load("fig_cost_cube.npz")
            for Jr in (J, 9):
                if f"tail_E{Jr}" in ref:
                    print(f"  difference from the stored field of the two-mirror solver, J = {Jr}: "
                          f"{np.linalg.norm(E - ref[f'tail_E{Jr}']) / np.linalg.norm(E):.2e}")
        else:
            import glob
            np.savez(f"cube_general_t_J{J}_obl.npz", E=E, pts=pts, d=d, p=p, ext=info["ext"], sca=info["sca"], t=t)
            for f in sorted(glob.glob("bem_cube_graded_*_obl.npz")):
                b = np.load(f)
                if np.allclose(b["d"], d) and np.allclose(b["p"], p):
                    print(f"  difference from {f} ({int(b['ndof'])} dofs): {np.linalg.norm(E - b['E']) / np.linalg.norm(E):.2e}; "
                          f"Bempp sigma_ext {float(b['ext']):.6f}")
