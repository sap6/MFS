"""Final timing of the cube comparison on one machine: the MFS levels of the
table and the best Bempp runs (with the mirror symmetry) are run alternately,
several times, and each keeps its smallest wall time. The cloud containers
differ in speed and drift during a session, so times measured at different
moments are not comparable (NOTES.md item 60).

    python3 cube_duel.py [repetitions]      # default 3; about 10 minutes per repetition

Updates the time columns of fig_cost_cube.npz ("tail", "dense") and writes
bem_grading_scan_k1_final.npz (lower envelope with the new times; the Bempp
errors are computed again against the reference "ref_E"). The MFS fields
are checked against the stored ones and not changed.
"""
import sys
import time

import numpy as np

import bem_env  # noqa: F401
import bempp_cl.api as bempp

import bem_cube_sym as b
import bem_grading_scan as scan
import fig_cost_cube as fc

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    d = dict(np.load("fig_cost_cube.npz"))
    pts = d["pts"]
    ref = d["ref_E"]                                    # reference field: J = 12 with the finer interior grid
    mfs = [("tail", int(r[0]), "balanced") for r in d["tail"]] + [("dense", int(r[0]), "dense") for r in d["dense"]]
    fin = np.load("bem_grading_scan_k1_final.npz")
    rows = [tuple(r) for r in fin["all"]]
    cand = [r for r in rows if r[6] < 2.5e-3 and not any(q[5] <= 0.8 * r[5] and q[6] <= r[6] for q in rows)]
    for N in sorted({r[4] for r in rows}):              # and the most accurate run on each mesh (fig_cube_convergence.py)
        best = min((r for r in rows if r[4] == N), key=lambda r: r[6])
        if best not in cand:
            cand.append(best)
    cand.sort(key=lambda r: (r[1], r[5]))
    e_bem = {}
    print(f"{len(mfs)} MFS runs and {len(cand)} Bempp runs, {reps} repetitions", flush=True)
    t_mfs = {m: [] for m in mfs}
    t_bem = {r: [] for r in cand}
    fc.run_mfs(2, pts)
    fc.run_mfs(2, pts, family="dense")
    for rep in range(reps):
        t0 = time.time()
        for m in mfs:
            E, t, unk, nrow, opt = fc.run_mfs(m[1], pts, family=m[2])
            t_mfs[m].append(t)
            assert np.linalg.norm(E - d[f"{m[0]}_E{m[1]}"]) <= 1e-9 * np.linalg.norm(E), m
        reg = None
        for r in cand:
            if int(r[1]) != reg:
                reg = int(r[1])
                bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
                bempp.GLOBAL_PARAMETERS.quadrature.singular = reg + 2
                b.solve(4, 1.0, pts)
            E, N, sizes, ts, grid = b.solve(int(r[3]), r[0], pts)
            t_bem[r].append(sum(ts))
            e_bem[r] = np.linalg.norm(E - ref) / np.linalg.norm(ref)     # against the present reference
        print(f"repetition {rep + 1}: {time.time() - t0:.0f} s", flush=True)
        for key in ("tail", "dense"):
            a = d[key].copy()
            for i, row in enumerate(a):
                a[i, 1] = min(t_mfs[(key, int(row[0]), "balanced" if key == "tail" else "dense")])
            d[key] = a
        np.savez("fig_cost_cube.npz", **d)
        out = [(r[0], r[1], r[2], r[3], r[4], min(t_bem[r]), e_bem[r]) for r in cand]
        np.savez("bem_grading_scan_k1_final.npz", rows=np.array(scan.envelope(out)), all=np.array(out))
    print("MFS: family, J, times of the repetitions, smallest, error")
    for m in mfs:
        e = np.linalg.norm(d[f"{m[0]}_E{m[1]}"] - ref) / np.linalg.norm(ref)
        print(f"  {m[2]:8s} {m[1]:2d}  " + " ".join(f"{t:6.1f}" for t in t_mfs[m]) + f"   {min(t_mfs[m]):6.1f} s  {e:.2e}")
    print("Bempp with symmetry: g, quadrature, dofs, times, smallest, error")
    for r in cand:
        print(f"  g = {r[0]:4.2f} {int(r[1])}/{int(r[2])} {int(r[4]):5d}  " + " ".join(f"{t:6.1f}" for t in t_bem[r]) + f"   {min(t_bem[r]):6.1f} s  {e_bem[r]:.2e}")
    print("final Bempp envelope:")
    for r in scan.envelope([(r[0], r[1], r[2], r[3], r[4], min(t_bem[r]), e_bem[r]) for r in cand]):
        print(f"   {r[5]:6.1f} s  {r[6]:.2e}   g = {r[0]:4.2f}  {int(r[1])}/{int(r[2])}  dofs {int(r[4])}")
