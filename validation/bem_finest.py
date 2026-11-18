"""Bempp on the finest cube mesh that fits into memory with the symmetry
(n = 48 intervals per edge, 41472 unknowns: the matrix of one octant of test
functions is 3.4 GB, the OpenCL driver allows 4 GB): the gradings and
quadrature orders near the best of the scan, each timed twice, appended to
the scan and merged into the timed runs of the case (NOTES.md item 87).

    python3 bem_finest.py k1|k1o|k2|k3 [repetitions]
"""
import sys

import numpy as np

import pec_cube_lightning3d as m

N_FINE = 48
GS = (1.5, 1.6, 1.75, 1.9)
case = sys.argv[1]
reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
m.K = 2 * np.pi * {"k2": 2, "k3": 3}.get(case, 1)
import cube_mfs as c  # noqa: E402
import bem_grading_scan as scan  # noqa: E402
import bem_cube_sym as b  # noqa: E402
import bempp_cl.api as bempp  # noqa: E402

b.K = m.K
pts = np.load("bem_cube_h0.05.npz")["pts"]
if case == "k1":
    ref, args = np.load("fig_cost_cube.npz")["ref_E"], (pts,)
else:
    ref = np.load(f"cube_cases_{case}_ref.npz")["E"]
    args = (pts,) + (tuple(c.wave(60.0, 25.0, 35.0)) if case == "k1o" else ())
fname = "bem_grading_scan_k1_final.npz" if case == "k1" else f"cube_cases_{case}.npz"
tkey, ekey = ("all", "rows") if case == "k1" else ("bem", "env")
d = dict(np.load(fname))
done = {(round(r[0], 3), int(r[1]), int(r[3])) for r in d[tkey]}
new = []
for reg in (4, 6, 8):
    bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
    bempp.GLOBAL_PARAMETERS.quadrature.singular = reg + 2
    b.solve(4, 1.0, *args)
    for g in GS:
        if (round(g, 3), reg, N_FINE) in done:
            continue
        ts = []
        for rep in range(reps):
            E, N, sizes, t4, grid = b.solve(N_FINE, g, *args)
            ts.append(sum(t4))
        row = (g, reg, reg + 2, N_FINE, N, min(ts), np.linalg.norm(E - ref) / np.linalg.norm(ref))
        new.append(row)
        print(f"{case} g = {g:4.2f} quadrature {reg}/{reg + 2} n = {N_FINE} dofs {N} " + " ".join(f"{t:6.1f}" for t in ts) + f" s  {row[6]:.2e}", flush=True)
        allrows = [tuple(r) for r in d[tkey]] + new
        out = dict(d)
        out[tkey], out[ekey] = np.array(allrows), np.array(scan.envelope(allrows))
        np.savez(fname, **out)
        s = np.load(f"bem_grading_scan_{case}.npz")["rows"]
        np.savez(f"bem_grading_scan_{case}.npz", rows=np.vstack((s, np.array(row)[None, :])))
print("envelope:")
for r in np.load(fname)[ekey][-6:]:
    print(f"   {r[5]:6.1f} s  {r[6]:.2e}   g = {r[0]:4.2f}  {int(r[1])}/{int(r[2])}  dofs {int(r[4])}")
