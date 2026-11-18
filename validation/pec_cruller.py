"""MFS for a perfectly conducting "cruller": a torus whose tube radius varies,

    x(theta, phi) = ((1 + H cos theta) cos phi, (1 + H cos theta) sin phi, H sin theta),
    H = 0.5 + 0.1 cos(5 phi + 3 theta),

a standard test surface of high-order integral equation codes (Barnett,
Greengard and Hagstrom, arXiv:1904.00076). It has no mirror plane and is not
a body of revolution. Sources: electric dipoles (3 unknowns each) on the
surface with H replaced by eta * H, on a tensor grid in (theta, phi).
Collocation on a finer tensor grid of the surface, two tangential components,
rows weighted by sqrt(area), dense least squares. No symmetry is used.

    python3 pec_cruller.py                 # known solution and plane wave, a few sizes
    python3 pec_cruller.py eta 0.3,0.5,0.7 # the source surface
"""
import sys
import time

import numpy as np

from lsq import factor
from pec_torus import Solution, cross_sections, dipole_field, plane_wave, rel, sphere_rule, wave

K = 2 * np.pi


def surface(th, ph, eta=1.0, amp=0.1):
    """Points, unit normals and area elements at the parameters (th, ph)."""
    u = 5 * ph + 3 * th
    H, Ht, Hp = eta * (0.5 + amp * np.cos(u)), -3 * eta * amp * np.sin(u), -5 * eta * amp * np.sin(u)
    ct, st, cp, sp = np.cos(th), np.sin(th), np.cos(ph), np.sin(ph)
    rho = 1 + H * ct
    x = np.stack((rho * cp, rho * sp, H * st), axis=-1)
    rt, rp = Ht * ct - H * st, Hp * ct                           # d rho / d theta, d rho / d phi
    xt = np.stack((rt * cp, rt * sp, Ht * st + H * ct), axis=-1)
    xp = np.stack((rp * cp - rho * sp, rp * sp + rho * cp, Hp * st), axis=-1)
    n = np.cross(xp, xt)                                         # outward
    J = np.linalg.norm(n, axis=-1)
    return x, n / J[..., None], J, xt


def grid(nt, nphi, eta=1.0, offset=0.5):
    th, ph = np.meshgrid(2 * np.pi * (np.arange(nt) + offset) / nt, 2 * np.pi * (np.arange(nphi) + offset) / nphi, indexing="ij")
    x, n, J, xt = surface(th.ravel(), ph.ravel(), eta)
    t1 = xt / np.linalg.norm(xt, axis=1)[:, None]
    T = np.stack((t1, np.cross(n, t1)), axis=1)                  # two orthonormal tangents
    return x, T, J * (2 * np.pi / nt) * (2 * np.pi / nphi)


def test_points(n=1500, seed=3):
    rng = np.random.default_rng(seed)
    th, ph = rng.uniform(0, 2 * np.pi, (2, n))
    x, nrm, J, xt = surface(th, ph)
    t1 = xt / np.linalg.norm(xt, axis=1)[:, None]
    return x, np.stack((t1, np.cross(nrm, t1)), axis=1), np.sqrt(J)


def tangential_matrix(x, T, y, k, batch=400):
    """Rows: the two tangential components of E at x; columns: the three components of each dipole at y."""
    A = np.empty((len(x), 2, len(y), 3), complex)
    for i in range(0, len(x), batch):
        dv = x[i:i + batch, None] - y[None]
        r = np.linalg.norm(dv, axis=2)
        kr = k * r
        g = np.exp(1j * kr) / r
        fI, fR = (1 + 1j / kr - 1 / kr ** 2) * g, (1 + 3j / kr - 3 / kr ** 2) * g
        Tr = np.einsum("mca,mna->mcn", T[i:i + batch], dv / r[..., None])
        A[i:i + batch] = fI[:, None, :, None] * T[i:i + batch, :, None, :] - (fR[:, None, :] * Tr)[..., None] * (dv / r[..., None])[:, None]
    return A.reshape(2 * len(x), 3 * len(y))


def solve(k, rhs_fun, nt, nphi, eta=0.5, over=2):
    """nt x nphi dipoles on the surface shrunk by eta; over^2 collocation points per dipole."""
    t0 = time.time()
    y = grid(nt, nphi, eta)[0]
    x, T, w = grid(over * nt, over * nphi)
    sw = np.repeat(np.sqrt(w), 2)
    A = tangential_matrix(x, T, y, k) * sw[:, None]
    b = np.einsum("mca,ma->mc", T, rhs_fun(x)).reshape(-1) * sw
    t1 = time.time()
    f = factor(A, "qr")
    alpha, res = f.solve(b, residual=True)
    info = dict(unknowns=3 * len(y), rows=len(b), t=time.time() - t0, t_assemble=t1 - t0, fit=float(res / np.linalg.norm(b)),
                alpha_norm=float(np.linalg.norm(alpha)))
    return Solution(y, alpha.reshape(-1, 3), k, info)


def residual(sol, rhs_fun, n=1500):
    xt, Tt, sw = test_points(n)
    got = np.einsum("mca,ma->mc", Tt, sol.field(xt)) * sw[:, None]
    want = np.einsum("mca,ma->mc", Tt, rhs_fun(xt)) * sw[:, None]
    return rel(got, want)


def scatter(k, d, p, nt, nphi, eta=0.5, over=2):
    rhs = lambda x: -plane_wave(x, k, d, p)
    sol = solve(k, rhs, nt, nphi, eta, over)
    sol.info["bc"] = residual(sol, rhs)
    ext, sca = cross_sections(sol, d, p)
    sol.info.update(sigma_ext=ext, sigma_sca=sca, opt=abs(ext - sca) / sca)
    return sol


SIZES = ((8, 24), (12, 36), (16, 48), (20, 60), (24, 72))

if __name__ == "__main__":
    d, p = wave(60.0, 25.0, 35.0)
    pts = 2.5 * sphere_rule(10, 20)[0]
    if len(sys.argv) > 1 and sys.argv[1] == "eta":
        for eta in (float(v) for v in sys.argv[2].split(",")):
            for nt, nphi in SIZES[:int(sys.argv[3]) if len(sys.argv) > 3 else 3]:
                i = scatter(K, d, p, nt, nphi, eta).info
                print(f"eta={eta:.2f} unknowns={i['unknowns']:6d} bc={i['bc']:.1e} optical={i['opt']:.1e} |alpha|={i['alpha_norm']:.1e} t={i['t']:.1f}s", flush=True)
        sys.exit()
    yk = np.array([[1.08, 0.0, 0.05], [-0.3, 0.95, -0.06], [0.0, -0.93, 0.02]])
    mk = np.array([[1.0, 0.5j, -0.3], [0.2, 1.0, 1.0j], [-1.0j, 0.3, 0.7]])
    known = lambda x: dipole_field(x, K, yk, mk)
    for nt, nphi in SIZES[:4]:
        s = solve(K, known, nt, nphi)
        print(f"known solution: unknowns={s.info['unknowns']:6d} bc={residual(s, known):.1e} field error at |x|=2.5: {rel(s.field(pts), known(pts)):.1e} t={s.info['t']:.1f}s", flush=True)
    for nt, nphi in SIZES[:4]:
        i = scatter(K, d, p, nt, nphi).info
        print(f"plane wave: unknowns={i['unknowns']:6d} bc={i['bc']:.1e} optical={i['opt']:.1e} sigma_ext={i['sigma_ext']:.9f} |alpha|={i['alpha_norm']:.1e} t={i['t']:.1f}s", flush=True)
