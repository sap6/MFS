"""Adds one more MFS level to a finished case of cube_cases.py (timed twice,
smaller time kept) without repeating the Bempp runs.

    python3 cube_cases_add_level.py k2 9
"""
import sys

import numpy as np

import pec_cube_lightning3d as m
import cube_cases as cc

case, J = sys.argv[1], int(sys.argv[2])
S = cc.SIDE[case]
m.K = 2 * np.pi * S
import cube_mfs as c  # noqa: E402

pts = np.load("bem_cube_h0.05.npz")["pts"]
dp = c.wave(*cc.OBLIQUE) if case.endswith("o") else (m.D, m.P)
ref = cc.reference(case)
c.run(2, dp[0], dp[1], pts, family="dense", **cc.layout(c, S, 2))
ts = []
for rep in range(2):
    E, t, info = c.run(J, dp[0], dp[1], pts, family="dense", **cc.layout(c, S, J))
    ts.append(t)
row = (1, J, max(info["per_class"]), info["unknowns"], info["rows"], min(ts), np.linalg.norm(E - ref) / np.linalg.norm(ref), info["optical"])
print(f"{case} J = {J}: unknowns/class {row[2]} rows {row[4]} times {ts[0]:.1f} {ts[1]:.1f} s  error {row[6]:.2e} optical {row[7]:.1e} |c| {info['coef_norm']:.1f}")
d = dict(np.load(f"cube_cases_{case}.npz"))
d["mfs"] = np.array(sorted([tuple(r) for r in d["mfs"] if int(r[1]) != J] + [row], key=lambda r: r[1]))
np.savez(f"cube_cases_{case}.npz", **d)
