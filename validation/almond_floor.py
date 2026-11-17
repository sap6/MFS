"""Almond, nose-on, 1.19 GHz: is the self-convergence of the skeleton layouts
real? Fields of several layouts at 400 points, their mutual differences and
optical mismatches: more rows towards the focal curve (J), two rows per
octave (q = 2^-1/2), closer spacing along the rows (s), finer collocation.

    python3 almond_floor.py
"""
import itertools
import time

import numpy as np

import pec_almond_skeleton as a
from pec_sphere_convergence import fibonacci_sphere

S = 2 ** -0.5
CASES = {                     # name: (n_mid, J, q, s, s_col)
    "6:3": (6, 3, 0.5, 0.7, 0.35),
    "8:5": (8, 5, 0.5, 0.7, 0.35),
    "10:7": (10, 7, 0.5, 0.7, 0.35),
    "8:5 s0.5": (8, 5, 0.5, 0.5, 0.35),
    "8:5 col0.25": (8, 5, 0.5, 0.7, 0.25),
    "8:10 q.71": (8, 10, S, 0.7, 0.35),
    "8:10 q.71 s0.5": (8, 10, S, 0.5, 0.35),
    "10:14 q.71 s0.5": (10, 14, S, 0.5, 0.35),
}

if __name__ == "__main__":
    k = 2 * np.pi * 1.19e9 / 299792458.0
    pts = 0.3 * fibonacci_sphere(400)
    a.solve(k, a.source_rows(2, 0, 0.5))
    E = {}
    for name, (n_mid, J, q, s, s_col) in CASES.items():
        t0 = time.time()
        y = a.source_rows(n_mid, J, q, s=s)
        c, info = a.solve(k, y, s_col=s_col)
        E[name] = (a.matrix(pts, y, k) @ c).reshape(-1, 3)
        ext, sca = a.cross_sections(k, y, c)
        print(f"{name:16s} sources={len(y):5d} unknowns={info['unknowns']:5d} rows={info['rows']:6d} {info['solver']:8s} t={time.time() - t0:5.1f}s "
              f"|c|={info['coef_norm']:.1e} optical={abs(ext - sca) / sca:.1e} sigma_ext={ext:.8e}", flush=True)
    names = list(E)
    print("differences:")
    for n1, n2 in itertools.combinations(names, 2):
        print(f"  {n1:16s} - {n2:16s}: {np.linalg.norm(E[n1] - E[n2]) / np.linalg.norm(E[n2]):.2e}")
    np.savez("almond_floor_t.npz", **{n.replace(" ", "_"): v for n, v in E.items()})
