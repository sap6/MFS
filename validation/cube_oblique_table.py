"""Cube, oblique wave (from theta 60, phi 25, polarisation angle 35 degrees):
errors of the MFS levels and of the graded Bempp meshes, and the times.

The reference is the family with two source depths per octave at J = 10,
whose eight classes are solved in groups (memory); it is not timed. The
timed rows are the family with one depth per octave, J = 4..9, and the dense
family up to J = 8 (the largest whose eight classes fit in memory together).

    python3 cube_oblique_table.py run    # solves everything (about 12 min) and writes cube_oblique.npz
    python3 cube_oblique_table.py        # table from cube_oblique.npz and bem_cube_graded_*_obl.npz
"""
import glob
import sys

import numpy as np

import cube_mfs as cm

ANG = (60.0, 25.0, 35.0)
JS = (4, 5, 6, 7, 8, 9)
JD = (4, 5, 6, 7, 8)
J_REF = 10


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def pack(t, info):
    return np.array([t, info["unknowns"], info["rows"], info["t_assemble"], info["t_factor"], info["optical"],
                     info["ext"], info["residual"]])


if __name__ == "__main__":
    d, p = cm.wave(*ANG)
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        pts = np.load("bem_cube_h0.05.npz")["pts"]
        cm.run(2, d, p, pts)
        cm.run(2, d, p, pts, family="dense")
        out = {"pts": pts, "d": d, "p": p}
        if "resume" in sys.argv:
            out = dict(np.load("cube_oblique.npz"))
        for fam, Js in (("", JS), ("dense", JD)):
            for J in Js:
                if f"E{fam}{J}" in out:
                    continue
                E, t, info = cm.run(J, d, p, pts, family=fam or "balanced")
                out[f"E{fam}{J}"], out[f"info{fam}{J}"] = E, pack(t, info)
                print(f"{fam or 'fast':5s} J={J} time={t:.1f}s optical={info['optical']:.1e}", flush=True)
                np.savez("cube_oblique.npz", **out)
        E, t, info = cm.run(J_REF, d, p, pts, family="dense")
        out["Eref"], out["inforef"] = E, pack(t, info)
        np.savez("cube_oblique.npz", **out)
    r = np.load("cube_oblique.npz")
    ref = r["Eref"]
    print(f"reference: two depths per octave, J = {J_REF}: optical mismatch {r['inforef'][5]:.1e}, sigma_ext {r['inforef'][6]:.6f}")
    for fam, Js, name in (("", JS, "one depth per octave"), ("dense", JD, "two depths per octave")):
        print(f"MFS, eight classes, {name}:   J  unknowns(all classes)  rows/class  time[s]  (assembly, factorisations)  error  optical  sigma_ext")
        for J in Js:
            if f"E{fam}{J}" not in r:
                continue
            t, unk, rows, ta, tf, opt, ext, res = r[f"info{fam}{J}"]
            print(f"                     {J:2d} {int(unk):10d} {int(rows):12d} {t:9.1f}   ({ta:.1f}, {tf:.1f}) {rel(r[f'E{fam}{J}'], ref):12.2e} {opt:9.1e} {ext:.6f}")
    print("Bempp, graded meshes: n  quadrature  unknowns  time[s]  error  sigma_ext")
    for f in sorted(glob.glob("bem_cube_graded_*_obl.npz"), key=lambda f: ("_q" in f, int(np.load(f)["ndof"]))):
        b = np.load(f)
        if np.allclose(b["d"], d) and np.allclose(b["p"], p):
            q = f.split("_q")[1].split("_obl")[0].replace("_", "/") if "_q" in f else "4/6"
            print(f"                     {int(b['n']):2d} {q:>8s} {int(b['ndof']):10d} {float(b['t_solve']):8.1f} {rel(b['E'], ref):12.2e}   {float(b['ext']):.6f}")
