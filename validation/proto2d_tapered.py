"""2D prototype (unit square, TE or TM): which distribution of source depths
on the corner bisectors gives the best far field for a given number of
depths? Free solves (no tails), deepest depth 1e-9.

  geometric   d_j = d0 q^j                       (fixed ratio: a fixed attainable accuracy, NOTES item 53)
  tapered     d_j = d0 exp(-sigma (sqrt(n) - sqrt(j)))   (lightning solvers: dense near d0, sparse near the corner)

    python3 proto2d_tapered.py TE|TM
"""
import itertools
import sys

import numpy as np
from scipy.linalg import lstsq

import proto2d_square as p2

H = p2.H
D0 = 0.2 * np.sqrt(2)
DEEP = 1e-9


def solve(pol, depths, beta, n_base=10, p=4, col_q=0.5, per_octave=1.0):
    h = H - 0.2
    s = (np.arange(n_base) + 0.5) / n_base * 2 * h - h
    pts = []
    for ax, sg in itertools.product((0, 1), (-1, 1)):
        y = np.zeros((n_base, 2))
        y[:, ax], y[:, 1 - ax] = sg * h, s
        pts.append(y)
    for d in depths:
        for s1, s2 in itertools.product((-1, 1), repeat=2):
            pts.append(np.array([[s1 * (H - d / np.sqrt(2)), s2 * (H - d / np.sqrt(2))]]))
    y = np.concatenate(pts)
    lv = max(3, int(np.ceil(np.log(min(depths) / H) / np.log(col_q))) + 2)
    x, w, nrm, _ = p2.boundary(lv, max(p, int(np.ceil(p * per_octave))), col_q)
    rw = np.sqrt(w) * p2.corner_dist(x) ** beta
    nn = nrm if pol == "TE" else None
    A, b = p2.kernels(x, y, nn, True) * rw[:, None], -p2.incident(x, nn) * rw
    sc = np.linalg.norm(A, axis=0)
    c, _, rank, _ = lstsq(A / sc, b, lapack_driver="gelsy", cond=1e-13)
    res = np.linalg.norm(A / sc @ c - b) / np.linalg.norm(b)
    return p2.far_field(y, c / sc, True), dict(unk=A.shape[1], rows=A.shape[0], res=res, cnorm=np.linalg.norm(c))


def geometric(q):
    n = int(np.ceil(np.log(DEEP / D0) / np.log(q)))
    return D0 * q ** np.arange(n + 1)


def tapered(n):
    sigma = np.log(D0 / DEEP) / (np.sqrt(n) - 1)
    return D0 * np.exp(-sigma * (np.sqrt(n) - np.sqrt(np.arange(1, n + 1))))


if __name__ == "__main__":
    pol = sys.argv[1]
    beta = 0.5 if pol == "TE" else 0.0
    cases = [("geometric, ratio 1/2", geometric(0.5), 1.0), ("geometric, ratio 0.6", geometric(0.6), 1.0),
             ("geometric, ratio 2^-1/2", geometric(2 ** -0.5), 2.0), ("geometric, ratio 0.8", geometric(0.8), 3.0)]
    cases += [(f"tapered, n = {n}", tapered(n), 2.0 if n <= 60 else 3.0) for n in (25, 30, 40, 50, 60, 80, 100)]
    out = {}
    for name, depths, po in cases:
        out[name] = solve(pol, depths, beta, per_octave=po)
    ref = out["tapered, n = 100"][0]
    alt = out["geometric, ratio 0.8"][0]
    print(f"{pol}: far field against the tapered solve with 100 depths (and against geometric 0.8: {np.linalg.norm(alt - ref) / np.linalg.norm(ref):.1e})")
    for name, depths, po in cases:
        ff, info = out[name]
        print(f"  {name:26s} depths={len(depths):4d} unknowns={info['unk']:5d} rows={info['rows']:5d} residual={info['res']:.1e} "
              f"|c|={info['cnorm']:.1e} error={np.linalg.norm(ff - ref) / np.linalg.norm(ref):.1e} (vs geometric 0.8: {np.linalg.norm(ff - alt) / np.linalg.norm(alt):.1e})")
