"""Compiled assembly of the field of a line dipole distribution.

For every target the line integral is evaluated with a quadrature adapted to
that target (panels in the stretched variable u growing dyadically away from
the nearest point of the line), and accumulated against Legendre modes of the
dipole density. Equivalent to Line._block in pec_cube_lightning3d.py but
without temporaries and with one panel count per target.
"""
import numpy as np
from numba import njit, prange


@njit(parallel=True, fastmath=True, cache=True)
def line_block(x, T, use_T, ax, b, c, yb, yc, te, lam, d, n, g, gw, cap, magnetic, k, H, out):
    M = x.shape[0]
    ng = g.shape[0]
    n_cap = int(np.ceil(2.0 / cap))
    nc = 2 if use_T else 3
    for m in prange(M):
        xa, xb, xc = x[m, ax], x[m, b], x[m, c]
        t0 = min(max(xa, -te), te)
        db, dc = xb - yb, xc - yc
        r = np.sqrt(db * db + dc * dc + (xa - t0) ** 2)
        u0 = np.log((H + t0) / (H - t0)) / lam
        w0 = min(r * (2 * H / (H * H - t0 * t0)) / lam, cap)
        n_dy = max(1, int(np.ceil(np.log2(cap / w0))) + 1)
        n_pan = n_dy + n_cap
        P = np.empty(n)
        G = np.zeros((3, 3), np.complex128)
        GT = np.zeros((3, 3), np.complex128)
        dv = np.empty(3)
        acc = np.zeros((nc, n, 3), np.complex128)
        for side in range(2):
            sgn = -1.0 if side == 0 else 1.0
            prev = 0.0
            for ip in range(n_pan):
                if ip < n_dy:
                    cur = min(w0 * 2.0 ** ip, cap)
                else:
                    cur = cap * (2 + ip - n_dy)
                cur = max(cur, prev)
                a0 = min(max(u0 + sgn * prev, -1.0), 1.0)
                b0 = min(max(u0 + sgn * cur, -1.0), 1.0)
                prev = cur
                half = 0.5 * abs(b0 - a0)
                if half < 1e-300:
                    continue
                mid = 0.5 * (a0 + b0)
                for iq in range(ng):
                    u = mid + half * g[iq]
                    ch = np.cosh(0.5 * u * lam)
                    t = H * np.tanh(0.5 * u * lam)
                    wt = half * gw[iq] * H * 0.5 * lam / (ch * ch) / d
                    dv[ax] = xa - t
                    dv[b] = db
                    dv[c] = dc
                    rr = np.sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2])
                    kr = k * rr
                    ph = (np.cos(kr) + 1j * np.sin(kr)) / rr
                    if magnetic:
                        f = (1j - 1.0 / kr) * ph
                        e0, e1, e2 = dv[0] / rr, dv[1] / rr, dv[2] / rr
                        G[0, 0] = 0.0; G[0, 1] = -f * e2; G[0, 2] = f * e1
                        G[1, 0] = f * e2; G[1, 1] = 0.0; G[1, 2] = -f * e0
                        G[2, 0] = -f * e1; G[2, 1] = f * e0; G[2, 2] = 0.0
                    else:
                        fI = (1.0 + 1j / kr - 1.0 / (kr * kr)) * ph
                        fR = (1.0 + 3j / kr - 3.0 / (kr * kr)) * ph
                        for i in range(3):
                            for j in range(3):
                                G[i, j] = -fR * dv[i] * dv[j] / (rr * rr)
                            G[i, i] += fI
                    if use_T:
                        for a in range(2):
                            for j in range(3):
                                GT[a, j] = T[m, a, 0] * G[0, j] + T[m, a, 1] * G[1, j] + T[m, a, 2] * G[2, j]
                    else:
                        for a in range(3):
                            for j in range(3):
                                GT[a, j] = G[a, j]
                    P[0] = 1.0
                    if n > 1:
                        P[1] = u
                    for kk in range(2, n):
                        P[kk] = ((2 * kk - 1) * u * P[kk - 1] - (kk - 1) * P[kk - 2]) / kk
                    for a in range(nc):
                        for kk in range(n):
                            pw = wt * P[kk]
                            for j in range(3):
                                acc[a, kk, j] += pw * GT[a, j]
        for a in range(nc):
            for kk in range(n):
                for j in range(3):
                    out[m, a, kk, j] = acc[a, kk, j]
