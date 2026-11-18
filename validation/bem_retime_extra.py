"""Bempp runs of the grading scan that were left out of the timed comparison
because the scan had been timed with a slower version of the symmetry fold
(NOTES.md item 72): time every scanned run that could now lie on the lower
envelope and merge it into the data of the case.

    python3 bem_retime_extra.py k1|k1o|k2|k3 [repetitions] [list]

A scanned run is timed if its scan time, scaled by the smallest new/old
ratio of the timed runs with the same quadrature and mesh, is not beaten by
a timed run that is at least as accurate and 1.25 times faster; if no run
with its quadrature and mesh was timed, it is timed in any case. Writes
bem_grading_scan_k1_final.npz (k1) or cube_cases_<case>.npz (k1o, k2).
"""
import sys

import numpy as np

import pec_cube_lightning3d as m

case = sys.argv[1]
reps = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 2
S = {"k2": 2, "k3": 3}.get(case, 1)
m.K = 2 * np.pi * S
import cube_mfs as c  # noqa: E402
import bem_grading_scan as scan  # noqa: E402

fname = "bem_grading_scan_k1_final.npz" if case == "k1" else f"cube_cases_{case}.npz"
d = dict(np.load(fname))
tkey = "all" if case == "k1" else "bem"
timed = [tuple(r) for r in d[tkey]]
rows = [tuple(r) for r in np.load(f"bem_grading_scan_{case}.npz")["rows"]]
key = lambda r: (round(r[0], 3), int(r[1]), int(r[3]))
tk = {key(r): r for r in timed}
ratio = {}
for r in rows:
    if key(r) in tk:
        g = (int(r[1]), int(r[3]))
        ratio[g] = min(ratio.get(g, 9.9), tk[key(r)][5] / r[5])
for reg in (4, 6, 8):                                   # groups without a timed run: time them all (ratio 0)
    for n in sorted({int(r[3]) for r in rows}):
        ratio.setdefault((reg, n), 0.0)
todo = []
for r in rows:
    if key(r) in tk:
        continue
    est = r[5] * ratio[(int(r[1]), int(r[3]))]
    if not any(q[5] <= est / 1.25 and q[6] <= r[6] for q in timed):
        todo.append(r)
todo.sort(key=lambda r: (r[1], r[5]))
print(f"{case}: {len(timed)} timed runs, {len(todo)} more to time ({sum(r[5] * ratio[(int(r[1]), int(r[3]))] for r in todo):.0f} s per repetition)", flush=True)
if "list" in sys.argv:
    for r in todo:
        print(f"   g = {r[0]:4.2f} {int(r[1])}/{int(r[2])} n = {int(r[3])} dofs {int(r[4])} scan {r[5]:.1f} s {r[6]:.2e}")
    sys.exit()
import bem_cube_sym as b  # noqa: E402
import bempp_cl.api as bempp  # noqa: E402
b.K = m.K
pts = np.load("bem_cube_h0.05.npz")["pts"]
if case == "k1":
    ref, args = np.load("fig_cost_cube.npz")["ref_E"], (pts,)
else:
    ref = np.load(f"cube_cases_{case}_ref.npz")["E"]
    args = (pts,) + (tuple(c.wave(60.0, 25.0, 35.0)) if case == "k1o" else ())
t, e = {r: [] for r in todo}, {}
for rep in range(reps):
    reg = None
    for r in todo:
        if int(r[1]) != reg:
            reg = int(r[1])
            bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
            bempp.GLOBAL_PARAMETERS.quadrature.singular = reg + 2
            b.solve(4, 1.0, *args)
        E, N, sizes, ts, grid = b.solve(int(r[3]), r[0], *args)
        t[r].append(sum(ts))
        e[r] = np.linalg.norm(E - ref) / np.linalg.norm(ref)
    new = timed + [(r[0], r[1], r[2], r[3], r[4], min(t[r]), e[r]) for r in todo]
    d[tkey] = np.array(new)
    d["rows" if case == "k1" else "env"] = np.array(scan.envelope(new))
    np.savez(fname, **d)
    print(f"repetition {rep + 1} done", flush=True)
for r in todo:
    print(f"   g = {r[0]:4.2f} {int(r[1])}/{int(r[2])} dofs {int(r[4])}: {min(t[r]):6.1f} s  {e[r]:.2e}   (scan: {r[5]:.1f} s)")
print("envelope:")
for r in d["rows" if case == "k1" else "env"]:
    print(f"   {r[5]:6.1f} s  {r[6]:.2e}   g = {r[0]:4.2f}  {int(r[1])}/{int(r[2])}  dofs {int(r[4])}")
