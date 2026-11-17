"""Cube MFS, balanced family: time split and error against the stored J = 9
reference field, for the assembly selected by MFS_LINES (tree or old).

    python3 cube_speed.py 4,5,6 [PM ALPHA_M NG ALPHA_G] [col_p]
"""
import sys
import time

import numpy as np

import fig_cost_cube as f
import line_tree as lt
import pec_cube_lightning3d as m

if __name__ == "__main__":
    Js = [int(v) for v in sys.argv[1].split(",")]
    if len(sys.argv) > 5:
        lt.PM, lt.ALPHA_M, lt.NG, lt.ALPHA_G = int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
    col_p = int(sys.argv[6]) if len(sys.argv) > 6 else f.COL["col_p"]
    d = np.load("fig_cost_cube.npz")
    pts, ref = d["pts"], d["tail_E9"]
    f.run_mfs(2, pts)
    for J in Js:
        t0 = time.time()
        kw, hmax = f.balanced(J)
        src = m.SymSources(J, tail_n=9, **kw)
        lv = max(2, int(np.ceil(np.log(src.depths[-1] / (2 * m.H)) / np.log(f.COL["col_q"]))) + 1)
        c, info = m.solve(src, lv, col_p, f.COL["col_q"], f.COL["beta"], hmax=hmax)
        t2 = time.time()
        E = (src.matrix(pts) @ c).reshape(-1, 3)
        t3 = time.time()
        ext, sca = m.cross_sections(src, c)
        old = d[f"tail_E{J}"] if f"tail_E{J}" in d else ref
        print(f"J={J} lines={m.TREE and 'tree' or 'old'} unknowns={info['unknowns']} rows={info['rows']} assemble={info['t_assemble']:.1f} "
              f"qr={info['t_solve']:.1f} eval={t3 - t2:.1f} total={t3 - t0:.1f}s  error vs J9: {f.rel(E, ref):.3e}  "
              f"change from stored field: {f.rel(E, old):.1e}  optical {abs(ext - sca) / sca:.1e} |c|={info['coef_norm']:.1e}", flush=True)
