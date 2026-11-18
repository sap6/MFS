"""Check of the shared-node line assembly (line_tree.py) against the
target-adapted quadrature it replaces (line_kernel.py), on the collocation
points of the cube solver and on far points.

    python3 test_line_tree.py [J]

Prints, for every kind of line of the balanced family at depth J, the largest
difference between the two blocks relative to the largest entry, with the old
rule at its production setting (12 Gauss nodes per panel) and at a tighter
one (20 nodes, panels half as wide), which serves as the reference.
"""
import sys
import time

import numpy as np

import pec_cube_lightning3d as m
from pec_sphere_convergence import fibonacci_sphere
from fig_cost_cube import COL, balanced


def old_block(line, x, T, ng=12, cap=0.17):
    """Block from line_kernel with a chosen rule."""
    g, gw = np.polynomial.legendre.leggauss(ng)
    b = line.base if isinstance(line, m.TailLine) else line
    saved = b.g, b.gw, b.cap
    b.g, b.gw, b.cap = g, gw, cap
    tree = m.TREE
    m.TREE = False
    try:
        return line.block(x, T)
    finally:
        m.TREE = tree
        b.g, b.gw, b.cap = saved


if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    kw, hmax = balanced(J)
    src = m.Sources(J, tail_n=9, **kw)
    src.lines = [ln for ln in src.lines if ln.yb > 0 and ln.yc > 0]      # the lines of one octant, as in cube_mfs.py
    lv = max(2, int(np.ceil(np.log(src.depths[-1] / (2 * m.H)) / np.log(COL["col_q"]))) + 1)
    x, w, T, cells, hmin = m.cube_cells(lv, COL["col_p"], COL["col_q"], hmax)
    rng = np.random.default_rng(0)
    sel = rng.choice(len(x), 3000, replace=False)
    near = np.argsort(m.edge_dist(x))[:1500]                 # the points closest to the edges
    idx = np.unique(np.concatenate((sel, near)))
    x, T = x[idx], T[idx]
    print(f"J = {J}: {len(src.lines)} lines, {len(x)} targets, smallest distance to an edge {m.edge_dist(x).min():.2e}")
    worst = 0.0
    for line in src.lines:
        if isinstance(line, m.VertexTail):
            continue
        t0 = time.time()
        new = line.block(x, T)
        t1 = time.time()
        old = old_block(line, x, T)
        t2 = time.time()
        ref = old_block(line, x, T, 20, 0.085)
        sc = np.abs(ref).max(axis=(0, 1, 3), keepdims=True)   # per mode
        e_new, e_old = np.abs((new - ref) / sc).max(), np.abs((old - ref) / sc).max()
        worst = max(worst, e_new)
        kind = ("tail " if isinstance(line, m.TailLine) else "line ") + ("m" if line.kernel is m.magnetic_pairs else "e")
        b = line.base if isinstance(line, m.TailLine) else line
        print(f"{kind} ax={b.ax} d={b.d:.2e} n={b.n:2d}: tree-ref {e_new:.1e}  old-ref {e_old:.1e}   "
              f"time tree {t1 - t0:.2f}s old {t2 - t1:.2f}s", flush=True)
    far = 2.0 * fibonacci_sphere(400)
    for line in src.lines[:4]:
        new, ref = line.block(far), old_block(line, far, None, 20, 0.085)
        worst = max(worst, np.abs(new - ref).max() / np.abs(ref).max())
    print(f"worst relative difference tree - reference: {worst:.1e}")
