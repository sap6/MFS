"""2D prototype: singular "tails" for the corner sources of the unit square.

The free corner sources of proto2d_square.py stop at depth d_J. Below that
depth the density on the bisector follows the corner expansion of the field
(powers rho^(2n/3) and integers), so the sources at depths d_j < d_J are not
given free coefficients: they are tied together in a few fixed combinations

    T_e = sum_{j > J} (d_j / d_J)^e  * source(d_j),     e = 2/3, 1, 4/3, ...

one unknown per exponent e and per source type (monopole / normal dipole).

Usage: python3 proto2d_square_tail.py TE|TM
"""
import itertools
import sys

import numpy as np
from scipy.linalg import lstsq

import proto2d_square as p2

H, K = p2.H, p2.K
D0 = 0.2 * np.sqrt(2)


def tail_columns(x, nrm, J, q, n_tail, exps_m, exps_d, qt=None):
    """Columns for the tails of the four corners. Returns matrix, plus a
    function giving their far field."""
    qt = qt or q
    dJ = D0 * q ** J
    d = dJ * qt ** np.arange(1, n_tail + 1)
    cols, ffs = [], []
    th = np.arange(360) / 360 * 2 * np.pi
    xh = np.column_stack((np.cos(th), np.sin(th)))
    for s1, s2 in itertools.product((-1, 1), repeat=2):
        y = np.column_stack((s1 * (H - d / np.sqrt(2)), s2 * (H - d / np.sqrt(2))))
        Kf = p2.kernels(x, y, nrm, dip=True)                  # [mono | dip_x | dip_y]
        n = len(d)
        mono, dx, dy = Kf[:, :n], Kf[:, n:2 * n], Kf[:, 2 * n:]
        ph = np.exp(-1j * K * (xh @ y.T))
        fm, fx, fy = ph, -1j * ph * xh[:, :1], -1j * ph * xh[:, 1:]
        # unit normal to the bisector (s1, s2)/sqrt2 is (s1, -s2)/sqrt2
        dn, fn = (s1 * dx - s2 * dy) / np.sqrt(2), (s1 * fx - s2 * fy) / np.sqrt(2)
        for e in exps_m:
            w = (d / dJ) ** e
            cols.append(mono @ w)
            ffs.append(fm @ w)
        for e in exps_d:
            w = (d / dJ) ** e / dJ
            cols.append(dn @ w * dJ)
            ffs.append(fn @ w * dJ)
    if not cols:
        return np.zeros((len(x), 0)), np.zeros((360, 0))
    return np.column_stack(cols), np.column_stack(ffs)


def solve(pol, J, q=0.5, beta=None, n_tail=40, exps_m=(), exps_d=(), col_extra=3, p=4, n_base=10,
          cond=1e-13, qt=None):
    beta = (0.0 if pol == "TM" else 0.5) if beta is None else beta
    y = p2.sources(J, q, (0.0,), n_base)
    d_min = D0 * q ** J
    lv = max(3, int(np.ceil(np.log(d_min / H) / np.log(0.5))) + col_extra)
    x, w, nrm, _ = p2.boundary(lv, p, 0.5)
    rw = np.sqrt(w) * p2.corner_dist(x) ** beta
    nn = nrm if pol == "TE" else None
    A = p2.kernels(x, y, nn, True)
    b = -p2.incident(x, nn)
    At, Ft = tail_columns(x, nn, J, q, n_tail, exps_m, exps_d, qt)
    A = np.column_stack((A, At)) * rw[:, None]
    b = b * rw
    sc = np.linalg.norm(A, axis=0)
    c, _, rank, _ = lstsq(A / sc, b, lapack_driver="gelsy", cond=cond)
    res = np.linalg.norm(A / sc @ c - b) / np.linalg.norm(b)
    c = c / sc
    nf = A.shape[1] - At.shape[1]
    ff = p2.far_field(y, c[:nf], True) + Ft @ c[nf:]
    return ff, dict(unk=A.shape[1], rows=A.shape[0], res=res)


EXP_M = (2 / 3, 1, 4 / 3, 5 / 3, 2, 7 / 3, 8 / 3, 3)
EXP_D = (5 / 3, 2, 7 / 3, 8 / 3, 3, 10 / 3)

if __name__ == "__main__":
    pol = sys.argv[1]
    yr, cr, _ = p2.solve(pol, J=40, q=0.6, beta=0.0 if pol == "TM" else 0.5, dip=True)   # 70 depths (down to 1e-16) break down with numpy 2.5 / scipy 1.18; 30, 40 and 50 agree to 3e-7
    ref = p2.far_field(yr, cr, True)
    print("far-field error against the J=40, q=0.6 free solve")
    print(" J  d_J      free     " + "  ".join(f"tail{n}m{m}d" for n, m in ((1, 0), (2, 1), (3, 1), (4, 2), (6, 4), (8, 6))))
    for J in (2, 4, 6, 8, 10, 12):
        line = f"{J:2d}  {D0 * 0.5 ** J:.1e} "
        for nm, nd in ((0, 0), (1, 0), (2, 1), (3, 1), (4, 2), (6, 4), (8, 6)):
            ff, info = solve(pol, J, exps_m=EXP_M[:nm], exps_d=EXP_D[:nd])
            line += f" {np.linalg.norm(ff - ref) / np.linalg.norm(ref):.1e}"
        print(line, flush=True)
