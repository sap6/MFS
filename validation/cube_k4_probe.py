"""Cube of four wavelengths (k = 8 pi): does the layout of the paper carry
over, and how many more modes per line does it need? (NOTES.md item 86.)
Self-convergence only: no independent reference exists at this size.

    python3 cube_k4_probe.py J1,J2,... dn [nb=8] [hmax=0.1] [gb=2.5] [colp=4]
"""
import os
import sys

import numpy as np

import pec_cube_lightning3d as m

m.K = 8 * np.pi
import cube_mfs as c  # noqa: E402

if __name__ == "__main__":
    dn = int(sys.argv[2])
    opt = dict(a.split("=") for a in sys.argv[3:])
    nb, hmax, gb = int(opt.get("nb", 8)), float(opt.get("hmax", 0.1)), float(opt.get("gb", 2.5))
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    store = dict(np.load("cube_k4_t_probe.npz")) if os.path.exists("cube_k4_t_probe.npz") else {}
    c.run(2, pts=pts, family="dense")
    for J in (int(v) for v in sys.argv[1].split(",")):
        kw = dict(c.dense(J), hmax=hmax)
        kw["kw"] = dict(kw["kw"], n0=kw["kw"]["n0"] + dn, n_base=nb)
        if "colp" in opt:                                     # Gauss nodes per boundary panel (default 4)
            kw["col"] = dict(kw["col"], col_p=int(opt["colp"]))
        E, t, info = c.run(J, pts=pts, family="dense", max_gb=gb, **kw)
        name = f"J{J},dn={dn},nb={nb},h={hmax:g}" + (f",p={opt['colp']}" if "colp" in opt else "")
        store[name] = E
        np.savez("cube_k4_t_probe.npz", **store)
        print(f"{name}: unknowns/class {max(info['per_class'])} rows {info['rows']} t = {t:5.1f} s  optical {info['optical']:.1e}  |c| {info['coef_norm']:.1f}  "
              f"sigma_ext {info['ext']:.6f}", flush=True)
