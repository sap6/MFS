"""2D prototype of the edge problem: plane-wave scattering by a unit square.

TM: u = E_z, Dirichlet (u_sc = -u_inc), sources = monopoles H0(k r).
TE: u = H_z, Neumann (d_n u_sc = -d_n u_inc); this is the analogue of fitting
    the tangential electric field, which is singular like rho^(-1/3).

Sources: a base ring of points on an inner square plus points clustered
towards each corner along one or more rays. Used to find a source layout and
row weighting that converge, before porting to the 3D cube.
"""
import itertools
import sys

import numpy as np
from scipy.linalg import lstsq
from scipy.special import hankel1

K = 2 * np.pi
H = 0.5
TH = np.deg2rad(200.0)                      # incidence direction
DIR = np.array([np.cos(TH), np.sin(TH)])


def boundary(levels, p, q=0.5, hmax=0.1):
    dist = np.concatenate(([0.0], H * q ** np.arange(levels, -1, -1)))
    brk = [0.0]
    for a, b in zip(dist[:-1], dist[1:]):
        n = int(np.ceil((b - a) / hmax))
        brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
    brk = np.array(brk)
    full = np.concatenate((-H + brk, (H - brk[::-1])[1:]))
    g, w = np.polynomial.legendre.leggauss(p)
    s = np.concatenate([0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(full[:-1], full[1:])])
    ws = np.concatenate([0.5 * (b - a) * w for a, b in zip(full[:-1], full[1:])])
    pts, nrm = [], []
    for ax, sg in itertools.product((0, 1), (-1, 1)):
        x = np.zeros((len(s), 2))
        x[:, ax], x[:, 1 - ax] = sg * H, s
        n = np.zeros((len(s), 2))
        n[:, ax] = sg
        pts.append(x)
        nrm.append(n)
    x = np.concatenate(pts)
    return x, np.tile(ws, 4), np.concatenate(nrm), H - np.abs(x).min(axis=1) * 0 - np.sort(np.abs(x), axis=1)[:, 0] * 0 - (np.abs(x).min(axis=1))


def corner_dist(x):
    return H - np.abs(x).min(axis=1)


def sources(J, q, angles=(0.0,), n_base=10, depth0=0.2, d_first=None):
    h = H - depth0
    s = (np.arange(n_base) + 0.5) / n_base * 2 * h - h
    base = []
    for ax, sg in itertools.product((0, 1), (-1, 1)):
        y = np.zeros((n_base, 2))
        y[:, ax], y[:, 1 - ax] = sg * h, s
        base.append(y)
    pts = [np.concatenate(base)]
    d0 = d_first or depth0 * np.sqrt(2)
    for j in range(J + 1):
        d = d0 * q ** j
        for s1, s2 in itertools.product((-1, 1), repeat=2):
            for a in angles:
                ang = np.pi / 4 + np.deg2rad(a)
                pts.append(np.array([[s1 * (H - d * np.cos(ang)), s2 * (H - d * np.sin(ang))]]))
    return np.concatenate(pts)


def kernels(x, y, nrm=None, dip=False):
    """Monopole (and optionally dipole) fields; with nrm, their normal derivatives."""
    dv = x[:, None, :] - y[None, :, :]
    r = np.linalg.norm(dv, axis=2)
    h0, h1 = hankel1(0, K * r), hankel1(1, K * r)
    cols = []
    if nrm is None:
        cols.append(h0)
        if dip:
            for c in range(2):
                cols.append(h1 * dv[..., c] / r)               # -grad_y-type dipoles (up to sign/k)
    else:
        rn = np.einsum("mnc,mc->mn", dv, nrm) / r
        cols.append(-K * h1 * rn)
        if dip:
            h2 = hankel1(2, K * r)
            for c in range(2):
                e = dv[..., c] / r
                # d/dn of (h1(kr) e_c) = k h1'(kr) rn e_c + h1 (n_c - rn e_c)/r
                dh1 = h0 - h1 / (K * r)
                cols.append(K * dh1 * rn * e + h1 * (nrm[:, None, c] - rn * e) / r)
    return np.concatenate(cols, axis=1)


def incident(x, nrm=None):
    u = np.exp(1j * K * (x @ DIR))
    return u if nrm is None else 1j * K * (nrm @ DIR) * u


def far_field(y, c, dip, n=360):
    th = np.arange(n) / n * 2 * np.pi
    xh = np.column_stack((np.cos(th), np.sin(th)))
    ph = np.exp(-1j * K * (xh @ y.T))
    cols = [ph]
    if dip:
        for cc in range(2):
            cols.append(-1j * ph * xh[:, cc:cc + 1])            # far field of h1 * e_c
    return np.concatenate(cols, axis=1) @ c


def solve(pol, J, q, beta, angles=(0.0,), dip=False, n_base=10, col_q=0.5, p=4, cond=1e-13):
    y = sources(J, q, angles, n_base)
    d_min = 0.2 * np.sqrt(2) * q ** max(J, 0)
    lv = max(3, int(np.ceil(np.log(d_min / H) / np.log(col_q))) + 2)
    x, w, nrm, _ = boundary(lv, p, col_q)
    rw = np.sqrt(w) * corner_dist(x) ** beta
    if pol == "TM":
        A, b = kernels(x, y, None, dip), -incident(x)
    else:
        A, b = kernels(x, y, nrm, dip), -incident(x, nrm)
    A, b = A * rw[:, None], b * rw
    sc = np.linalg.norm(A, axis=0)
    c, _, rank, _ = lstsq(A / sc, b, lapack_driver="gelsy", cond=cond)
    res = np.linalg.norm(A / sc @ c - b) / np.linalg.norm(b)
    return y, c / sc, {"N": len(y), "unk": A.shape[1], "rows": A.shape[0], "rank": int(rank), "res": res}


def study(pol, configs, ref):
    yr, cr, _ = solve(pol, **ref)
    fr = far_field(yr, cr, ref.get("dip", False))
    for cfg in configs:
        y, c, info = solve(pol, **cfg)
        f = far_field(y, c, cfg.get("dip", False))
        print(f"{pol} {cfg}: unk={info['unk']} rows={info['rows']} res={info['res']:.1e} "
              f"far-field diff vs ref={np.linalg.norm(f - fr) / np.linalg.norm(fr):.2e}", flush=True)


if __name__ == "__main__":
    pol = sys.argv[1]
    ref = dict(J=60, q=0.6, beta=0.0 if pol == "TM" else 0.5, dip=True, angles=(0.0,))
    cfgs = []
    for J, q in ((-1, 0.5), (3, 0.25), (6, 0.25), (6, 0.5), (12, 0.5), (24, 0.5), (20, 0.7), (40, 0.7)):
        for dip in (False, True):
            cfgs.append(dict(J=J, q=q, beta=ref["beta"], dip=dip))
    study(pol, cfgs, ref)
