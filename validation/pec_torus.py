"""MFS for exterior scattering of a plane wave by a perfectly conducting torus.

Surface: (rho - R)^2 + z^2 = a^2 (axis of revolution = z). Sources: electric
dipoles (3 unknowns each) on rings; the rings sit on a closed curve in the
meridian half-plane (default: a circle of radius a_s about rho = R_s, i.e. an
inner torus). Boundary condition: tangential E_sc = -tangential E_inc at
collocation points on rings of the surface, least squares, rows weighted by
sqrt(area).

Both point sets are tensor grids (ring index, angle about the axis), so the
matrix is block-circulant and the least-squares problem splits exactly into one
small problem per Fourier mode about the axis (`solve`). `solve_dense` builds
the full matrix and is only there to check that the two agree. The incident
wave may come from any direction: the splitting uses the symmetry of the body,
not of the wave.

Error measures (no closed-form solution exists for the torus):
  bc    tangential boundary residual at independent test points (relative L2)
  opt   optical theorem mismatch |sigma_ext - sigma_sca| / sigma_sca
  known a field with known exterior solution (dipoles inside the tube)
and comparison with Bempp (`bem_torus.py`).

Usage: python3 pec_torus.py            # self-checks, a few seconds
"""
import time

import numpy as np
from scipy.linalg import lstsq

RFAR = 1.0e6


def wave(theta_deg=0.0, phi_deg=0.0, pol_deg=90.0):
    """Plane wave arriving *from* the direction (theta, phi). pol is the angle
    of the polarisation from the theta-unit vector towards the phi-unit vector.
    Default: travelling along -z, polarised along y (the sphere and cube case)."""
    th, ph, al = np.radians([theta_deg, phi_deg, pol_deg])
    src = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    e_th = np.array([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)])
    e_ph = np.array([-np.sin(ph), np.cos(ph), 0.0])
    return -src, np.cos(al) * e_th + np.sin(al) * e_ph


def plane_wave(x, k, d, p):
    return p[None, :] * np.exp(1j * k * (x @ d))[:, None]


def dyadic(x, y, k):
    """(I + grad grad / k^2) e^{ikr}/r, shape (len(x), len(y), 3, 3)."""
    dv = x[:, None] - y[None]
    r = np.linalg.norm(dv, axis=2)
    kr = k * r
    g = np.exp(1j * kr) / r
    fI = (1 + 1j / kr - 1 / kr ** 2) * g
    fR = (1 + 3j / kr - 3 / kr ** 2) * g
    rh = dv / r[..., None]
    return fI[..., None, None] * np.eye(3) - fR[..., None, None] * rh[..., :, None] * rh[..., None, :]


# ---------------------------------------------------------------- geometry

def circle(n, centre, radius, offset=0.5):
    """n points on a circle in the meridian half-plane: (rho, z, angle)."""
    t = 2 * np.pi * (np.arange(n) + offset) / n
    return centre + radius * np.cos(t), radius * np.sin(t), t


def ring_points(rho, z, nphi, shift=0.0):
    """Rings of nphi points about the z axis. Returns (len(rho), nphi, 3) and the angles."""
    ph = 2 * np.pi * (np.arange(nphi) + shift) / nphi
    x = np.stack((rho[:, None] * np.cos(ph), rho[:, None] * np.sin(ph), z[:, None] + 0 * ph), axis=2)
    return x, ph


def cyl_frame(ph):
    """Unit vectors e_rho, e_phi, e_z at angles ph: (len(ph), 3, 3), rows are the vectors."""
    c, s, o = np.cos(ph), np.sin(ph), 0 * ph
    return np.stack((np.stack((c, s, o), 1), np.stack((-s, c, o), 1), np.stack((o, o, o + 1), 1)), 1)


class Torus:
    def __init__(self, R=1.0, a=0.5):
        self.R, self.a = R, a
        self.c = np.sqrt(R * R - a * a)       # focal ring of the toroidal coordinates

    def surface(self, ntheta, nphi, shift=0.0, offset=0.5):
        """Collocation rings: points (ntheta, nphi, 3), tangents (ntheta, nphi, 2, 3),
        area weight per ring."""
        rho, z, t = circle(ntheta, self.R, self.a, offset)
        x, ph = ring_points(rho, z, nphi, shift)
        F = cyl_frame(ph)                                      # (nphi, 3, 3)
        t_pol = -np.sin(t)[:, None, None] * F[None, :, 0] + np.cos(t)[:, None, None] * F[None, :, 2]
        t_tor = np.broadcast_to(F[None, :, 1], t_pol.shape)
        w = self.a * rho * (2 * np.pi / ntheta) * (2 * np.pi / nphi)
        return x, np.stack((t_pol, t_tor), axis=2), w

    def test_points(self, n=1500, seed=3):
        """Independent points on the surface (area-uniform), with tangents."""
        rng = np.random.default_rng(seed)
        t = np.empty(0)
        while len(t) < n:                                       # rejection: density ~ rho
            c = rng.uniform(0, 2 * np.pi, 2 * n)
            keep = rng.uniform(0, self.R + self.a, 2 * n) < self.R + self.a * np.cos(c)
            t = np.concatenate((t, c[keep]))
        t = t[:n]
        ph = rng.uniform(0, 2 * np.pi, n)
        rho = self.R + self.a * np.cos(t)
        x = np.column_stack((rho * np.cos(ph), rho * np.sin(ph), self.a * np.sin(t)))
        F = cyl_frame(ph)
        t_pol = -np.sin(t)[:, None] * F[:, 0] + np.cos(t)[:, None] * F[:, 2]
        return x, np.stack((t_pol, F[:, 1]), axis=1)


# ---------------------------------------------------------------- solver

class Solution:
    """Dipole positions y (N, 3) and moments alpha (N, 3), Cartesian."""

    def __init__(self, y, alpha, k, info):
        self.y, self.alpha, self.k, self.info = y, alpha, k, info

    def field(self, pts, batch=200):
        out = np.empty((len(pts), 3), complex)
        for i in range(0, len(pts), batch):
            out[i:i + batch] = np.einsum("pnab,nb->pa", dyadic(pts[i:i + batch], self.y, self.k), self.alpha)
        return out

    def far(self, dirs):
        """F with E_sc ~ F exp(ikr)/r in the directions dirs (unit vectors)."""
        ph = np.exp(-1j * self.k * (dirs @ self.y.T))            # (D, N)
        s = ph @ self.alpha
        return s - dirs * np.sum(dirs * s, axis=1, keepdims=True)


def sphere_rule(nt=48, nphi=96):
    g, w = np.polynomial.legendre.leggauss(nt)
    ph = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    ct, phh = np.meshgrid(g, ph, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    dirs = np.column_stack(((st * np.cos(phh)).ravel(), (st * np.sin(phh)).ravel(), ct.ravel()))
    return dirs, np.repeat(w, nphi) * 2 * np.pi / nphi


def cross_sections(sol, d, p):
    dirs, w = sphere_rule()
    F = sol.far(np.vstack((dirs, d[None])))
    ext = 4 * np.pi / sol.k * np.imag(np.conj(p) @ F[-1])
    sca = np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1))
    return float(ext), float(sca)


def solve(torus, k, rhs_fun, src_rho, src_z, nphi, ntheta_col=None, over_phi=1, rcond=None):
    """Least squares, one small problem per Fourier mode about the axis.

    src_rho, src_z: the source rings (meridian coordinates). nphi: dipoles per
    ring. Collocation: ntheta_col rings (default 2 x number of source rings) of
    over_phi * nphi points. rhs_fun(x) -> prescribed E (minus the incident field).
    """
    t0 = time.time()
    ns = len(src_rho)
    nc = ntheta_col or 2 * ns
    mphi = over_phi * nphi
    x, T, w = torus.surface(nc, mphi)
    sw = np.sqrt(w)
    # one source per ring (angle 0, staggered rings do not matter for accuracy here)
    y0 = np.column_stack((src_rho, 0 * src_rho, src_z))
    G = dyadic(x.reshape(-1, 3), y0, k).reshape(nc, mphi, ns, 3, 3)
    # rows: tangential components; columns: dipole along e_rho, e_phi, e_z at angle 0 = x, y, z
    B = np.einsum("ipca,ipjab->ipcjb", T, G) * sw[:, None, None, None, None]
    Bh = np.fft.fft(B, axis=1)                                   # (nc, mphi, 2, ns, 3) over the shift
    b = np.einsum("ipca,ipa->ipc", T, rhs_fun(x.reshape(-1, 3)).reshape(nc, mphi, 3)) * sw[:, None, None]
    bh = np.fft.fft(b, axis=1)                                   # (nc, mphi, 2)
    ah = np.empty((nphi, ns, 3), complex)
    smin, smax = np.inf, 0.0
    for m in range(nphi):
        mu = m + nphi * np.arange(over_phi)
        A = Bh[:, mu].transpose(1, 0, 2, 3, 4).reshape(over_phi * nc * 2, ns * 3)
        r = bh[:, mu].transpose(1, 0, 2).reshape(-1)
        sol, _, _, sv = lstsq(A, r, cond=rcond)
        ah[m] = sol.reshape(ns, 3)
        smin, smax = min(smin, sv[-1]), max(smax, sv[0])
    # back to dipoles: alpha_q = (1/nphi) sum_m ah_m exp(2 pi i m q / nphi), cylindrical components
    acyl = np.fft.ifft(ah, axis=0)                               # (nphi, ns, 3)
    yq, ph = ring_points(src_rho, src_z, nphi)                   # (ns, nphi, 3)
    F = cyl_frame(ph)                                            # (nphi, 3, 3)
    alpha = np.einsum("qjc,qca->jqa", acyl, F)                   # Cartesian, (ns, nphi, 3)
    info = {"unknowns": 3 * ns * nphi, "rows": 2 * nc * mphi, "t": time.time() - t0,
            "alpha_norm": float(np.linalg.norm(alpha)), "cond": float(smax / smin)}
    sol = Solution(yq.reshape(-1, 3), alpha.reshape(-1, 3), k, info)
    sol.modes = ah                                               # (mode about the axis, ring, [rho, phi, z])
    sol.rhs_modes = np.linalg.norm(bh[:, :nphi], axis=(0, 2))    # size of the data in each mode
    return sol


def solve_dense(torus, k, rhs_fun, src_rho, src_z, nphi, ntheta_col=None, over_phi=1):
    """The same least-squares problem with the full matrix (check only)."""
    ns = len(src_rho)
    nc = ntheta_col or 2 * ns
    x, T, w = torus.surface(nc, over_phi * nphi)
    sw = np.repeat(np.sqrt(w), over_phi * nphi)[:, None]
    x, T = x.reshape(-1, 3), T.reshape(-1, 2, 3)
    y = ring_points(src_rho, src_z, nphi)[0].reshape(-1, 3)
    A = np.einsum("mca,mnab->mcnb", T, dyadic(x, y, k)) * sw[:, :, None, None]
    b = np.einsum("mca,ma->mc", T, rhs_fun(x)) * sw
    alpha = lstsq(A.reshape(2 * len(x), 3 * len(y)), b.reshape(-1))[0].reshape(-1, 3)
    return Solution(y, alpha, k, {"unknowns": 3 * len(y)})


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def residual(torus, sol, rhs_fun, n=1500):
    xt, Tt = torus.test_points(n)
    got = np.einsum("mca,ma->mc", Tt, sol.field(xt))
    want = np.einsum("mca,ma->mc", Tt, rhs_fun(xt))
    return rel(got, want)


def scatter(torus, k, d, p, ns, nphi, a_s=None, R_s=None, **kw):
    """Plane-wave scattering with sources on the inner torus (R_s, a_s)."""
    a_s = 0.5 * torus.a if a_s is None else a_s
    R_s = torus.R if R_s is None else R_s
    rho, z, _ = circle(ns, R_s, a_s)
    rhs = lambda x: -plane_wave(x, k, d, p)
    sol = solve(torus, k, rhs, rho, z, nphi, **kw)
    sol.info["bc"] = residual(torus, sol, rhs)
    ext, sca = cross_sections(sol, d, p)
    sol.info.update(sigma_ext=ext, sigma_sca=sca, opt=abs(ext - sca) / sca)
    return sol


def dipole_field(x, k, y, mom):
    return np.einsum("pnab,nb->pa", dyadic(x, y, k), mom)


def main():
    k = 2 * np.pi
    tor = Torus(1.0, 0.5)
    d, p = wave()
    # 1. the mode-by-mode solver against the full matrix
    rho, z, _ = circle(8, 1.0, 0.25)
    rhs = lambda x: -plane_wave(x, k, d, p)
    s1 = solve(tor, k, rhs, rho, z, 24)
    s2 = solve_dense(tor, k, rhs, rho, z, 24)
    pts = 2.5 * sphere_rule(8, 16)[0]
    print(f"mode-by-mode against full matrix (576 unknowns): field differs by {rel(s1.field(pts), s2.field(pts)):.1e}")
    # 2. known exterior solution: dipoles inside the tube, off the centre circle
    yk = np.array([[1.08, 0.0, 0.05], [-0.3, 0.95, -0.06], [0.0, -0.93, 0.02]])
    mk = np.array([[1.0, 0.5j, -0.3], [0.2, 1.0, 1.0j], [-1.0j, 0.3, 0.7]])
    known = lambda x: dipole_field(x, k, yk, mk)
    for ns, nphi in ((12, 40), (20, 64), (28, 96), (36, 128)):
        rho, z, _ = circle(ns, 1.0, 0.25)
        s = solve(tor, k, known, rho, z, nphi)
        print(f"known solution: unknowns={s.info['unknowns']:6d} bc={residual(tor, s, known):.1e} "
              f"field error at |x|=2.5: {rel(s.field(pts), known(pts)):.1e}  t={s.info['t']:.2f}s")
    # 3. plane wave
    for ns, nphi in ((12, 40), (20, 64), (28, 96), (36, 128)):
        s = scatter(tor, k, d, p, ns, nphi, a_s=0.25)
        i = s.info
        print(f"plane wave: unknowns={i['unknowns']:6d} bc={i['bc']:.1e} optical={i['opt']:.1e} "
              f"sigma_ext={i['sigma_ext']:.10f} |alpha|={i['alpha_norm']:.1e} t={i['t']:.2f}s")


if __name__ == "__main__":
    main()
