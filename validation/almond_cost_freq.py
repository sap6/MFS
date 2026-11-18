"""Almond, nose-on incidence: error against wall time at a frequency other
than 1.19 GHz (default 3.5 GHz, where the almond is three wavelengths long),
MFS layouts against Bempp with the two mirror planes; the reference is the
finest MFS layout (8:5). Sources are spaced at 0.5 times their distance to
the surface, as in the benchmark runs. NOTES.md item 85.

    python3 almond_cost_freq.py [GHz]     # writes almond_cost_freq_<GHz>GHz.npz (unknowns, time, error)
"""
import sys
import time

import numpy as np

import pec_almond_skeleton as a
from pec_sphere_convergence import fibonacci_sphere

if __name__ == "__main__":
    freq = float(sys.argv[1]) if len(sys.argv) > 1 else 3.5
    k = 2 * np.pi * freq * 1e9 / 299792458.0
    pts = 0.3 * fibonacci_sphere(400)
    rel = lambda x, y: np.linalg.norm(x - y) / np.linalg.norm(y)
    a.solve(k, a.source_rows(2, 0, 0.5))
    out = []
    for n_mid, J in ((3, 0), (4, 1), (5, 2), (6, 3), (7, 4), (8, 5)):
        t0 = time.time()
        y = a.source_rows(n_mid, J, 0.5, s=0.5)
        c, info = a.solve(k, y)
        E = (a.matrix(pts, y, k) @ c).reshape(-1, 3)
        out.append((f"{n_mid}:{J}", info["unknowns"], time.time() - t0, E))
        print(f"MFS {n_mid}:{J} unknowns {info['unknowns']} rows {info['rows']} {out[-1][2]:.1f} s", flush=True)
    ref = out[-1][3]
    for name, unk, t, E in out[:-1]:
        print(f"MFS {name}: {unk} unknowns, {t:.1f} s, difference from 8:5: {rel(E, ref):.2e}")
    a.bem_reference(k, pts, ((12, 12),))
    bem = []
    for ndof, E, t in a.bem_reference(k, pts, ((28, 28), (40, 36), (60, 56), (84, 76), (100, 92))):
        bem.append((ndof, t, rel(E, ref)))
        print(f"Bempp {ndof} unknowns, {t:.1f} s, difference from the MFS 8:5: {bem[-1][2]:.2e}", flush=True)
    np.savez(f"almond_cost_freq_{freq:g}GHz.npz", mfs=np.array([(unk, t, rel(E, ref)) for name, unk, t, E in out[:-1]]), bem=np.array(bem))
