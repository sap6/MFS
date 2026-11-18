"""Three wavelengths with the layout of the paper as it stands (first depth
0.2, eight base points per side, interior panels 0.1, two depths per
octave), and 4 (S - 1) = 8 more modes per line: the rule of cube_cases.layout.
Stores the fields of the levels and prints their successive differences
(NOTES.md item 79).

    python3 cube_k3_levels.py J1,J2,...  [dn=8] [hmax=0.1] [nb=8] [gb=6]     (gb: memory for the matrices of one group of classes)
"""
import os
import sys

import numpy as np

import pec_cube_lightning3d as m

m.K = 6 * np.pi
import cube_mfs as c  # noqa: E402

if __name__ == "__main__":
    opt = dict(a.split("=") for a in sys.argv[2:])
    dn, hmax, nb = int(opt.get("dn", 8)), float(opt.get("hmax", 0.1)), int(opt.get("nb", 8))
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    bem = np.load("cube_k3_t_bem_n40.npy")
    store = dict(np.load("cube_k3_t_levels.npz")) if os.path.exists("cube_k3_t_levels.npz") else {}
    rel = lambda a, b: np.linalg.norm(a - b) / np.linalg.norm(b)
    c.run(2, pts=pts, family="dense")
    for J in (int(v) for v in sys.argv[1].split(",")):
        kw = dict(c.dense(J), hmax=hmax)
        kw["kw"] = dict(kw["kw"], n0=kw["kw"]["n0"] + dn, n_base=nb)
        E, t, info = c.run(J, pts=pts, family="dense", max_gb=float(opt.get("gb", 6.0)), **kw)
        name = f"J{J},dn={dn},h={hmax:g},nb={nb}"
        store[name] = E
        np.savez("cube_k3_t_levels.npz", **store)
        print(f"{name}: unknowns/class {max(info['per_class'])} rows {info['rows']} t = {t:5.1f} s  from Bempp {rel(E, bem):.2e}  "
              f"optical {info['optical']:.1e}  |c| {info['coef_norm']:.1f}", flush=True)
    names = sorted((k for k in store if k.endswith(f",dn={dn},h={hmax:g},nb={nb}")), key=lambda k: int(k.split(",")[0][1:]))
    for a, b in zip(names[:-1], names[1:]):
        print(f"   {a.split(',')[0]} -> {b.split(',')[0]}: {rel(store[a], store[b]):.2e}")
