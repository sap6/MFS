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
stored reference fields. The older solver for the standard wave with two
mirrors (pec_cube_lightning3d.SymSources) gives the same field to rounding.
"""
import itertools
import sys
import time

import numpy as np

import line_tree as lt
import pec_cube_lightning3d as m
from lsq import QR

COL = dict(col_p=3, col_q=0.5, beta=0.5)
REFL = [np.array(r, float) for r in itertools.product((1, -1), repeat=3)]
CLASSES = [tuple(float(v) for v in c) for c in itertools.product((1, -1), repeat=3)]


def balanced(J):
    """All resolution parameters grow with the depth J: modes per line
    (n0 + 2 log(1/d)), base grid, and the collocation that samples it. The
    tail exponents are the solver defaults (derived, DERIVATIONS.md B); the
    vertex tail (one exponent, DERIVATIONS.md D) is on."""
    n0 = max(4, 2 * J - 4)
    n_base, hmax = (4, 0.25) if J <= 5 else (6, 0.15) if J == 6 else (8, 0.1)
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


class Cube:
    """Assembles and factorises the given parity classes once; solve() then
    takes any plane wave whose classes are among them."""

    def __init__(self, J, kw=None, hmax=None, tail=9, classes=CLASSES, col=COL):
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
        A, off = self.matrix(self.x, self.T * self.rw[:, None, None], self.classes)
        self.t_assemble = time.time() - t0
        t0 = time.time()
        self.fac, self.scale = {}, {}
        for i, c in enumerate(self.classes):
            Ac = np.asfortranarray(A[:, off[i]:off[i + 1]])
            self.scale[c] = np.linalg.norm(Ac, axis=0)
            Ac /= self.scale[c]
            self.fac[c] = QR(Ac)
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
        """Rows: (target, component of T); columns: the unknowns of the given
        classes, one class after the other. T is (M, nc, 3), or None for the
        three Cartesian components. Returns the matrix and the column offsets."""
        ncls = len(classes)
        off = np.cumsum([0] + [self.unknowns[c] for c in classes])
        cmap = np.full((ncls, 3 * self.slot0[-1]), -1, np.int64)
        for i, c in enumerate(classes):
            cmap[i, self.keep[c]] = off[i] + np.arange(self.unknowns[c])
        nc = 3 if T is None else T.shape[1]
        A = np.zeros((len(x) * nc, off[-1]), complex)
        npt = len(self.pe) + len(self.pm)
        for R in REFL:
            sgn = np.array([sign(c, R) for c in classes])
            xr = x * R
            Tr = np.diag(R)[None] if T is None else T * R
            Tfull = np.broadcast_to(Tr, (len(x), nc, 3))
            if npt:                                    # point dipoles: electric, then magnetic
                G = [m.dyadic(xr, self.pe, m.K)]
                if len(self.pm):
                    G.append(m.magnetic_pairs(xr[:, None, :] - self.pm[None, :, :]))
                B = np.einsum("mab,mnbc->manc", Tfull, np.concatenate(G, axis=1)).reshape(len(x) * nc, -1)
                for i in range(ncls):
                    A[:, cmap[i, :3 * npt]] += sgn[i] * B
            i = 0
            while i < len(self.lines):
                ln = self.lines[i]
                lo = 3 * self.slot0[2 + i]
                if isinstance(ln, m.VertexTail):
                    B = ln.block(xr, Tfull).reshape(len(x) * nc, -1)
                    for j in range(ncls):
                        A[:, cmap[j, lo:lo + 3 * ln.n]] += sgn[j] * B
                    i += 1
                    continue
                mate = self.lines[i + 1] if i + 1 < len(self.lines) else None
                pair = mate is not None and m.same_place(ln, mate)
                tail = isinstance(ln, m.TailLine)
                b = ln.base if tail else ln
                mag = b.kernel is m.magnetic_pairs
                cm_first = np.ascontiguousarray(cmap[:, lo:lo + 3 * ln.n])
                cm_mate = np.ascontiguousarray(cmap[:, lo + 3 * ln.n:lo + 3 * (ln.n + mate.n)]) if pair else None
                lt.tree_add(xr, Tr, b.ax, b.b, b.c, np.array(ln.pos) if tail else (b.yb, b.yc), b.lam, b.d, b.n, m.K, m.H,
                            A, None if mag else cm_first, cm_first if mag else cm_mate, sgn,
                            rot=ln.rot if tail else None,
                            we=ln.w if tail and not mag else None,
                            wm=(mate.w if pair else ln.w) if tail and (pair or mag) else None)
                i += 2 if pair else 1
        return A, off

    def solve(self, d, p):
        """Coefficients per class for the wave (d, p), and the relative boundary residual."""
        coef, r2, b2 = {}, 0.0, 0.0
        for c in classes_of(d, p):
            if c not in self.fac:
                raise ValueError(f"the wave excites class {c}, which was not assembled")
            b = -np.einsum("mab,mb->ma", self.T, class_part(self.x, d, p, c)).reshape(-1) * np.repeat(self.rw, 2)
            y, _, res = self.fac[c].solve(b, residual=True)
            r2 += float(res) ** 2
            b2 += np.linalg.norm(b) ** 2
            coef[c] = y / self.scale[c]
        self.coef_norm = float(np.sqrt(sum(np.linalg.norm(v * self.scale[c]) ** 2 for c, v in coef.items())))
        return coef, np.sqrt(r2 / b2)

    def field(self, coef, pts):
        cl = list(coef)
        A, off = self.matrix(pts, None, cl)
        return (A @ np.concatenate([coef[c] for c in cl])).reshape(-1, 3)

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


def run(J, d=m.D, p=m.P, pts=None, **kw):
    """One wave: field at pts, total time, and a dictionary of diagnostics."""
    t0 = time.time()
    cube = Cube(J, classes=classes_of(d, p), **kw)
    coef, res = cube.solve(d, p)
    E = cube.field(coef, pts) if pts is not None else None
    t = time.time() - t0
    ext, sca = cube.cross_sections(coef, d, p)
    return E, t, dict(unknowns=sum(cube.unknowns.values()), per_class=list(cube.unknowns.values()), rows=cube.rows,
                      t_assemble=cube.t_assemble, t_factor=cube.t_factor, residual=res, ext=ext, sca=sca,
                      optical=abs(ext - sca) / sca, coef_norm=cube.coef_norm)


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
