"""Unit cube, any incident plane wave: MFS with edge sources and singular tail,
by parity classes.

The source layout has the mirrors x -> -x and y -> -y. Any incident field
splits into four parts of definite parity (sx, sy) under them, each part is
scattered into a field of the same parity, and each parity class is
represented by a quarter of the sources (pec_cube_lightning3d.SymSources).
The four class matrices are sums, with different signs, of the same four
reflected blocks, so the assembly is shared; each class then needs its own QR
factorisation, which serves every incident wave.

    python3 cube_general.py J [theta_deg phi_deg psi_deg]

runs one wave arriving from the direction (theta, phi) with polarisation
angle psi (0: theta-polarised), prints the optical-theorem mismatch, and
compares with bem_cube_graded_*_obl.npz if present. Without angles it
reproduces the standard case of the other scripts (wave along -z, E along y),
for which only the class (+1, -1) is excited.
"""
import sys
import time

import numpy as np
from lsq import QR

import pec_cube_lightning3d as m
from fig_cost_cube import COL, balanced

CLASSES = ((1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0))
REFL = (np.array([1.0, 1, 1]), np.array([-1.0, 1, 1]), np.array([1.0, -1, 1]), np.array([-1.0, -1, 1]))


def wave(theta, phi, psi):
    """Propagation direction d and polarisation p of a plane wave arriving from
    the direction (theta, phi) (degrees); psi rotates p from theta_hat towards phi_hat."""
    t, f, s = np.deg2rad([theta, phi, psi])
    rhat = np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
    that = np.array([np.cos(t) * np.cos(f), np.cos(t) * np.sin(f), -np.sin(t)])
    fhat = np.array([-np.sin(f), np.cos(f), 0.0])
    return -rhat, np.cos(s) * that + np.sin(s) * fhat


def plane(x, d, p):
    return p[None, :] * np.exp(1j * m.K * (x @ d))[:, None]


def class_part(x, d, p, cls):
    """Part of the plane wave with parity cls = (sx, sy), at the points x."""
    sig = (1.0, cls[0], cls[1], cls[0] * cls[1])
    return sum(s * plane(x * R, d, p) * R[None, :] for R, s in zip(REFL, sig)) / 4


class CubeSolver:
    """Factorises the four parity classes once; solve() then takes any plane waves."""

    def __init__(self, J, tail=9, classes=CLASSES, batch=1200):
        kw, hmax = balanced(J)
        self.src = {c: m.SymSources(J, tail_n=tail, parity=c, **kw) for c in classes}
        any_src = next(iter(self.src.values()))
        lv = max(2, int(np.ceil(np.log(any_src.depths[-1] / (2 * m.H)) / np.log(COL["col_q"]))) + 1)
        x, w, T, _, _ = m.cube_cells(lv, COL["col_p"], COL["col_q"], hmax)
        sel = m.quarter(x)
        self.x, self.T = x[sel], T[sel]
        self.rw = np.repeat(np.sqrt(w[sel]) * m.edge_dist(self.x) ** COL["beta"], 2)
        t0 = time.time()
        A = {c: np.empty((2 * len(self.x), s.n_unknowns), complex) for c, s in self.src.items()}
        n_keep, any_src.n_unknowns = any_src.n_unknowns, any_src._n_all
        try:                                               # the four reflected blocks, shared by all classes
            for i in range(0, len(self.x), batch):
                sl = slice(i, i + batch)
                blocks = [m.Sources.matrix(any_src, self.x[sl] * R, self.T[sl] * R, batch) for R in REFL]
                for c, s in self.src.items():
                    sig = (1.0, c[0], c[1], c[0] * c[1])
                    A[c][2 * i:2 * i + 2 * len(self.x[sl])] = sum(g * B[:, s.keep] for g, B in zip(sig, blocks))
        finally:
            any_src.n_unknowns = n_keep
        self.t_assemble = time.time() - t0
        t0 = time.time()
        self.QR, self.scale = {}, {}
        for c in self.src:
            A[c] *= self.rw[:, None]
            self.scale[c] = np.linalg.norm(A[c], axis=0)
            A[c] /= self.scale[c]
            self.QR[c] = QR(A[c])
            dR = self.QR[c].dR
            if dR.min() < 1e-12 * dR.max():
                print(f"warning: class {c} is numerically rank deficient ({dR.min() / dR.max():.1e})")
        self.t_factor = time.time() - t0
        self.unknowns = {c: s.n_unknowns for c, s in self.src.items()}
        self.rows = 2 * len(self.x)

    def solve(self, d, p):
        """Coefficients per class for the wave (d, p), and the relative boundary residual."""
        coef, r2, b2 = {}, 0.0, 0.0
        for c, fac in self.QR.items():
            b = -np.einsum("mab,mb->ma", self.T, class_part(self.x, d, p, c)).reshape(-1) * self.rw
            if np.linalg.norm(b) < 1e-13 * np.linalg.norm(self.rw):
                coef[c] = np.zeros(fac.n, complex)
                continue
            y, _, res = fac.solve(b, residual=True)
            r2 += res[()] ** 2
            b2 += np.linalg.norm(b) ** 2
            coef[c] = y / self.scale[c]
        return coef, np.sqrt(r2 / b2)

    def field(self, coef, pts):
        return sum((s.matrix(pts) @ coef[c]).reshape(-1, 3) for c, s in self.src.items() if np.any(coef[c]))

    def far(self, coef, dirs, rfar=1.0e6):
        return self.field(coef, rfar * dirs) * rfar * np.exp(-1j * m.K * rfar)

    def cross_sections(self, coef, d, p, nt=40, nphi=80):
        from bem_optical import sphere_rule
        dirs, w = sphere_rule(nt, nphi)
        F = self.far(coef, np.vstack((dirs, d[None, :])))
        ext = 4 * np.pi / m.K * float(np.imag(np.conj(p) @ F[-1]))
        return ext, float(np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1)))


if __name__ == "__main__":
    J = int(sys.argv[1])
    standard = len(sys.argv) < 5
    d, p = (m.D, m.P) if standard else wave(*(float(v) for v in sys.argv[2:5]))
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    t0 = time.time()
    cs = CubeSolver(J)
    coef, res = cs.solve(d, p)
    E = cs.field(coef, pts)
    t = time.time() - t0
    ext, sca = cs.cross_sections(coef, d, p)
    print(f"cube J={J} d={np.round(d, 4)} p={np.round(p, 4)}: unknowns per class {list(cs.unknowns.values())}, rows {cs.rows}, "
          f"assembly {cs.t_assemble:.0f}s, four factorisations {cs.t_factor:.0f}s, total {t:.0f}s")
    print(f"  boundary residual {res:.2e}, sigma_ext {ext:.6f}, sigma_sca {sca:.6f}, optical mismatch {abs(ext - sca) / sca:.1e}")
    np.savez(f"cube_general_t_J{J}_{'std' if standard else 'obl'}.npz", E=E, pts=pts, d=d, p=p, ext=ext, sca=sca, t=t)
    if standard:
        ref = np.load("fig_cost_cube.npz")
        for Jr in (J, 9):
            if f"tail_E{Jr}" in ref:
                print(f"  difference from the symmetric solver (fig_cost_cube.npz), J = {Jr}: "
                      f"{np.linalg.norm(E - ref[f'tail_E{Jr}']) / np.linalg.norm(E):.2e}")
    else:
        import glob
        for f in sorted(glob.glob("bem_cube_graded_*_obl.npz")):
            b = np.load(f)
            if np.allclose(b["d"], d) and np.allclose(b["p"], p):
                print(f"  difference from {f} ({int(b['ndof'])} dofs): {np.linalg.norm(E - b['E']) / np.linalg.norm(E):.2e}; "
                      f"Bempp sigma_ext {float(b['ext']):.6f}")
