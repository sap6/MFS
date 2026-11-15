"""2D prototype for the almond: plane-wave scattering by a 3:1 ellipse.

Compares two source layouts:
  confocal  monopoles on a confocal ellipse inside the scatterer
  skeleton  monopoles and dipoles on the interfocal segment (the "skeleton"),
            clustered towards the two foci, optionally extending past them
Errors: relative L2 boundary residual at independent test points.
"""
import sys
import numpy as np
from scipy.linalg import lstsq
from scipy.special import hankel1

A, B = 1.0, 1.0 / 3.0
C = np.sqrt(A * A - B * B)


def boundary(n):
    th = (np.arange(n) + 0.5) / n * 2 * np.pi
    x = np.column_stack((A * np.cos(th), B * np.sin(th)))
    nr = np.column_stack((B * np.cos(th), A * np.sin(th)))
    jac = np.linalg.norm(nr, axis=1)
    return x, nr / jac[:, None], jac * 2 * np.pi / n


def fields(x, y, k, nrm, dip):
    dv = x[:, None, :] - y[None, :, :]
    r = np.linalg.norm(dv, axis=2)
    h0, h1 = hankel1(0, k * r), hankel1(1, k * r)
    cols = []
    if nrm is None:
        cols.append(h0)
        if dip:
            cols += [h1 * dv[..., c] / r for c in range(2)]
    else:
        rn = np.einsum("mnc,mc->mn", dv, nrm) / r
        cols.append(-k * h1 * rn)
        if dip:
            dh1 = h0 - h1 / (k * r)
            for c in range(2):
                e = dv[..., c] / r
                cols.append(k * dh1 * rn * e + h1 * (nrm[:, None, c] - rn * e) / r)
    return np.concatenate(cols, axis=1)


def solve(pol, k, y, dip, n_col, th_inc=np.deg2rad(30)):
    d = np.array([np.cos(th_inc), np.sin(th_inc)])
    def rhs(x, n):
        u = np.exp(1j * k * (x @ d))
        return -u if pol == "TM" else -1j * k * (n @ d) * u
    x, n, w = boundary(n_col)
    Aop = fields(x, y, k, None if pol == "TM" else n, dip) * np.sqrt(w)[:, None]
    b = rhs(x, n) * np.sqrt(w)
    sc = np.linalg.norm(Aop, axis=0)
    c, *_ = lstsq(Aop / sc, b, lapack_driver="gelsy", cond=1e-14)
    xt, nt, wt = boundary(3 * n_col + 7)
    res = fields(xt, y, k, None if pol == "TM" else nt, dip) @ (c / sc) - rhs(xt, nt)
    return np.sqrt(np.sum(wt * np.abs(res) ** 2) / np.sum(wt * np.abs(rhs(xt, nt)) ** 2)), Aop.shape[1]


def confocal(n, tau):
    th = (np.arange(n) + 0.5) / n * 2 * np.pi
    return np.column_stack((np.sqrt(C * C + (tau * B) ** 2) * np.cos(th), tau * B * np.sin(th)))


def skeleton(n_mid, J, q, beyond=0.0):
    """Points on the x-axis: Chebyshev-like in the middle, then clustered
    geometrically towards each focus from inside; optionally a few beyond."""
    gap = A - C
    pts = list(C * np.cos((np.arange(n_mid) + 0.5) / n_mid * np.pi) * (1 - 2 * gap))
    for j in range(J + 1):
        dlt = 2 * gap * q ** j
        pts += [C - dlt, -(C - dlt)]
    for f in np.atleast_1d(beyond):
        if f > 0:
            pts += [C + f * gap, -(C + f * gap)]
    pts = np.unique(np.round(pts, 14))
    return np.column_stack((pts, np.zeros(len(pts))))


if __name__ == "__main__":
    for pol in ("TM", "TE"):
        for k in (1.2, 10.0):
            print(f"--- {pol}, k a = {k}")
            for n, tau in ((40, 0.5), (80, 0.5), (160, 0.5), (80, 0.2), (160, 0.2), (320, 0.2)):
                e, unk = solve(pol, k, confocal(n, tau), False, 4 * n)
                print(f"confocal tau={tau} unk={unk:4d}: {e:.1e}")
            for n_mid, J, q in ((10, 4, 0.5), (20, 8, 0.5), (30, 12, 0.5), (40, 20, 0.5), (40, 20, 0.7)):
                for dip in (False, True):
                    e, unk = solve(pol, k, skeleton(n_mid, J, q), dip, 600)
                    print(f"skeleton n_mid={n_mid} J={J} q={q} dip={dip} unk={unk:4d}: {e:.1e}")
