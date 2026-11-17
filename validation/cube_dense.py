"""Cube, standard wave: the family with two source depths per octave
(cube_mfs.dense), J = 3..10, timed; writes cube_dense.npz. Then errors of
both families and of the graded Bempp meshes against its finest member.

    python3 cube_dense.py run [Jmax]    # solve (about 5 min) and report
    python3 cube_dense.py               # report from the saved data
"""
import glob
import sys

import numpy as np

import cube_mfs as cm


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


if __name__ == "__main__":
    d0 = np.load("fig_cost_cube.npz")
    pts = d0["pts"]
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        jmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10
        cm.run(2, pts=pts, family="dense")
        out = {}
        for J in range(3, jmax + 1):
            E, t, info = cm.run(J, pts=pts, family="dense")
            out[f"E{J}"] = E
            out[f"info{J}"] = np.array([t, info["unknowns"], info["rows"], info["t_assemble"], info["t_factor"],
                                        info["optical"], info["ext"], info["residual"], info["coef_norm"]])
            print(f"dense J={J} unknowns={info['unknowns']} rows={info['rows']} time={t:.1f}s residual={info['residual']:.2e} "
                  f"optical={info['optical']:.1e} sigma_ext={info['ext']:.6f}", flush=True)
            np.savez("cube_dense.npz", **out)
    r = np.load("cube_dense.npz")
    Js = sorted(int(k[1:]) for k in r.files if k.startswith("E"))
    ref = r[f"E{Js[-1]}"]
    print(f"reference: dense family, J = {Js[-1]}")
    print("dense family:     J  unknowns/class  rows   time[s]  error      change to next   optical   sigma_ext   residual")
    for i, J in enumerate(Js):
        t, unk, rows, ta, tf, opt, ext, res, cn = r[f"info{J}"]
        nxt = f"{rel(r[f'E{J}'], r[f'E{Js[i + 1]}']):.2e}" if i + 1 < len(Js) else "   -    "
        print(f"                 {J:2d} {int(unk) // 2:10d} {int(rows):8d} {t:8.1f}  {rel(r[f'E{J}'], ref):.2e}   {nxt}       {opt:.1e}  {ext:.6f}  {res:.2e}")
    print("one depth per octave (fig_cost_cube.npz):")
    for J, t, unk, nrow, opt in d0["tail"]:
        print(f"                 {int(J):2d} {int(unk) // 2:10d} {int(nrow):8d} {t:8.1f}  {rel(d0[f'tail_E{int(J)}'], ref):.2e}                  {opt:.1e}")
    print("Bempp, graded meshes:  quadrature  unknowns  time[s]  error    sigma_ext")
    for q, tag in (("4/6", ""), ("6/8", "_q6_8"), ("8/10", "_q8_10")):
        for f in sorted(glob.glob(f"bem_cube_graded_n*_g2{tag}.npz"), key=lambda f: int(np.load(f)["ndof"])):
            b = np.load(f)
            print(f"                 {q:>9s} {int(b['ndof']):10d} {float(b['t_solve']):8.1f}  {rel(b['E'], ref):.2e}  {float(b['ext']):.6f}")
    fine = sorted((np.load(f) for f in glob.glob("bem_cube_graded_n*_g2_q8_10.npz")), key=lambda b: int(b["ndof"]))
    if len(fine) >= 3:
        from scipy.optimize import brentq
        (n1, E1), (n2, E2), (n3, E3) = ((float(b["n"]), b["E"]) for b in fine[-3:])
        q_ = np.linalg.norm(E1 - E2) / np.linalg.norm(E2 - E3)
        p = brentq(lambda p: (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p) - q_, 0.2, 12)
        Ex = E3 + (E3 - E2) * n3 ** -p / (n2 ** -p - n3 ** -p)
        s2, s3 = (float(b["ext"]) for b in fine[-2:])
        print(f"graded Bempp 8/10 extrapolated from its last three meshes (rate n^-{p:.2f}): error {rel(Ex, ref):.2e}, "
              f"sigma_ext {s3 + (s3 - s2) * n3 ** -p / (n2 ** -p - n3 ** -p):.6f}")
