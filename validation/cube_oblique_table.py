"""Cube, oblique wave (from theta 60, phi 25, polarisation angle 35 degrees):
errors of the MFS levels and of the graded Bempp meshes against the MFS run
with J = 9, and the times.

    python3 cube_oblique_table.py run    # solves J = 4..9 (about 4 min) and writes cube_oblique.npz
    python3 cube_oblique_table.py        # table from cube_oblique.npz and bem_cube_graded_*_obl.npz
"""
import glob
import sys

import numpy as np

import cube_mfs as cm

ANG = (60.0, 25.0, 35.0)
JS = (4, 5, 6, 7, 8, 9)


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


if __name__ == "__main__":
    d, p = cm.wave(*ANG)
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        pts = np.load("bem_cube_h0.05.npz")["pts"]
        cm.run(2, d, p, pts)
        out = {"pts": pts, "d": d, "p": p}
        for J in JS:
            E, t, info = cm.run(J, d, p, pts)
            out[f"E{J}"] = E
            out[f"info{J}"] = np.array([t, info["unknowns"], info["rows"], info["t_assemble"], info["t_factor"],
                                        info["optical"], info["ext"], info["residual"]])
            print(f"J={J} time={t:.1f}s optical={info['optical']:.1e}", flush=True)
            np.savez("cube_oblique.npz", **out)
    r = np.load("cube_oblique.npz")
    ref = r[f"E{JS[-1]}"]
    print("MFS, eight classes:   J  unknowns(all classes)  rows/class  time[s]  (assembly, factorisations)  error vs J=9  optical  sigma_ext")
    for J in JS:
        t, unk, rows, ta, tf, opt, ext, res = r[f"info{J}"]
        print(f"                     {J:2d} {int(unk):10d} {int(rows):12d} {t:9.1f}   ({ta:.1f}, {tf:.1f}) {rel(r[f'E{J}'], ref):16.2e} {opt:9.1e} {ext:.6f}")
    print("Bempp, graded meshes: n  quadrature  unknowns  time[s]  error vs MFS J=9  sigma_ext")
    for f in sorted(glob.glob("bem_cube_graded_*_obl.npz"), key=lambda f: ("_q" in f, int(np.load(f)["ndof"]))):
        b = np.load(f)
        if np.allclose(b["d"], d) and np.allclose(b["p"], p):
            q = f.split("_q")[1].split("_obl")[0].replace("_", "/") if "_q" in f else "4/6"
            print(f"                     {int(b['n']):2d} {q:>8s} {int(b['ndof']):10d} {float(b['t_solve']):8.1f} {rel(b['E'], ref):14.2e}      {float(b['ext']):.6f}")
