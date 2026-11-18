"""The cube beyond the standard case: a general incident wave, and a side of
two wavelengths. For each case the MFS levels (two depths per octave) and
the best Bempp runs (with the mirror symmetry; candidates from
bem_grading_scan.py) are run alternately, several times, and each keeps its
smallest wall time, as cube_duel.py does for the standard case (NOTES.md
items 60, 66, 68).

    python3 cube_cases.py k1o ref             # reference field (J = 11; J = 10 for k2) -> cube_cases_k1o_ref.npz; once
    python3 cube_cases.py k1o [repetitions]   # k = 2 pi, oblique wave (eight parity classes)
    python3 cube_cases.py k2 ref; python3 cube_cases.py k2 [repetitions]      # k = 4 pi, standard wave
    python3 cube_cases.py k3 ref; python3 cube_cases.py k3 [repetitions]      # k = 6 pi, standard wave
    python3 cube_cases.py report              # tables of all finished cases

Errors of both codes are differences from the reference field of the case.
Writes cube_cases_<case>.npz: mfs (1; J; unknowns per class; unknowns of all
classes; equations per class; time; error; optical mismatch), bem (all timed
Bempp runs) and env (lower envelope), both with the columns of
bem_grading_scan.py.
"""
import os
import sys

import numpy as np

import pec_cube_lightning3d as m

SIDE = {"k1o": 1, "k2": 2, "k3": 3}
OBLIQUE = (60.0, 25.0, 35.0)
LEVELS = {"k1o": range(4, 11), "k2": range(4, 9), "k3": range(4, 9)}     # two depths per octave
J_REFS = {"k1o": 11, "k2": 10, "k3": 10}                                 # level of the reference field
MODES = {1: 0, 2: 4, 3: 12}                                              # additional modes per line at S wavelengths per side


def layout(c, S, J, family="dense"):
    """The layout of cube_mfs.py at S wavelengths per side: the same rule as at
    one wavelength (first depth 0.2, eight base points per side), with MODES[S]
    more modes per line (NOTES.md items 57, 68, 79). The line quadrature
    takes more nodes by itself (line_tree.nodes). At three wavelengths the
    matrices of the two classes are built one after the other from J = 9 on,
    and for the general wave at one wavelength up to 8.5 GB of matrices are
    held at once (peak memory 12 GB), so that J = 9 needs one assembly."""
    if S == 1:
        return dict(max_gb=8.5)                 # the four matrices of a general wave at J = 9 (8 GB) are kept together
    kw = dict(c.dense(J))
    kw["kw"] = dict(kw["kw"], n0=kw["kw"]["n0"] + MODES[S])
    if S == 3 and J >= 9:
        kw["max_gb"] = 2.5
    return kw


def reference(case):
    return np.load(f"cube_cases_{case}_ref.npz")["E"]


def report():
    import bem_grading_scan as scan
    for case in SIDE:
        f = f"cube_cases_{case}.npz"
        if not os.path.exists(f):
            continue
        d = np.load(f)
        print(f"{case}: MFS  J  unknowns/class  all  equations  time  error  optical")
        for r in d["mfs"]:
            print(f"   {int(r[1]):2d} {int(r[2]):6d} {int(r[3]):6d} {int(r[4]):6d} {r[5]:7.1f} s  {r[6]:.2e}  {r[7]:.1e}")
        if "env" in d:
            print("   Bempp with symmetry, lower envelope:")
            for r in d["env"]:
                print(f"   {r[5]:7.1f} s  {r[6]:.2e}   g = {r[0]:4.2f}  {int(r[1])}/{int(r[2])}  dofs {int(r[4])}")


if __name__ == "__main__":
    case = sys.argv[1]
    if case == "report":
        report()
        sys.exit()
    S = SIDE[case]
    J_REF = J_REFS[case]
    m.K = 2 * np.pi * S                                     # before any solver object is built
    import cube_mfs as c
    import bem_grading_scan as scan
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    oblique = case.endswith("o")
    dp = c.wave(*OBLIQUE) if oblique else (m.D, m.P)
    c.run(2, dp[0], dp[1], pts, family="dense", **layout(c, S, 2))          # compilation
    if len(sys.argv) > 2 and sys.argv[2] == "ref":
        E, t, info = c.run(J_REF, dp[0], dp[1], pts, family="dense", **layout(c, S, J_REF))
        print(f"{case} reference, two depths per octave, J = {J_REF}: unknowns/class {max(info['per_class'])} rows {info['rows']} t = {t:.0f} s "
              f"optical {info['optical']:.1e} sigma_ext {info['ext']:.6f} |c| {info['coef_norm']:.1f}")
        np.savez(f"cube_cases_{case}_ref.npz", E=E, pts=pts, info=np.array([J_REF, t, max(info["per_class"]), info["rows"], info["optical"], info["ext"]]))
        if case != "k3":
            old = scan.reference(case)[1]
            print(f"   differs from the reference of the earlier family by {np.linalg.norm(E - old) / np.linalg.norm(old):.2e}")
        sys.exit()
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    ref = reference(case)
    runs = list(LEVELS[case])
    import bem_cube_sym as b
    import bempp_cl.api as bempp
    b.K = m.K
    rows = [tuple(r) for r in np.load(f"bem_grading_scan_{case}.npz")["rows"]]
    cand = [r for r in rows if not any(q[5] <= 0.5 * r[5] and q[6] <= r[6] for q in rows)]      # generous: the scan was timed before the fold was compiled
    cand.sort(key=lambda r: (r[1], r[5]))
    if "CUBE_CASES_TEST" in os.environ:                    # smoke test: two runs of each code
        runs, cand = runs[:2], cand[:2]
    print(f"{case}: {len(runs)} MFS runs and {len(cand)} Bempp runs, {reps} repetitions", flush=True)
    t_mfs, t_bem, res, e_bem = {J: [] for J in runs}, {r: [] for r in cand}, {}, {}
    args = (pts, dp[0], dp[1]) if oblique else (pts,)
    for rep in range(reps):
        for J in runs:
            E, t, info = c.run(J, dp[0], dp[1], pts, family="dense", **layout(c, S, J))
            t_mfs[J].append(t)
            res[J] = (max(info["per_class"]), info["unknowns"], info["rows"], np.linalg.norm(E - ref) / np.linalg.norm(ref), info["optical"])
            if rep == 0:
                print(f"  J = {J:2d}: unknowns/class {res[J][0]} rows {info['rows']} t = {t:6.1f} s  error {res[J][3]:.2e} "
                      f"optical {info['optical']:.1e} |c| {info['coef_norm']:.1f}", flush=True)
        reg = None
        for r in cand:
            if int(r[1]) != reg:
                reg = int(r[1])
                bempp.GLOBAL_PARAMETERS.quadrature.regular = reg
                bempp.GLOBAL_PARAMETERS.quadrature.singular = reg + 2
                b.solve(4, 1.0, *args)
            E, N, sizes, ts, grid = b.solve(int(r[3]), r[0], *args)
            t_bem[r].append(sum(ts))
            e_bem[r] = np.linalg.norm(E - ref) / np.linalg.norm(ref)       # against the present reference
        bem = [(r[0], r[1], r[2], r[3], r[4], min(t_bem[r]), e_bem[r]) for r in cand]
        out = dict(mfs=np.array([(1, J, res[J][0], res[J][1], res[J][2], min(t_mfs[J]), res[J][3], res[J][4]) for J in runs], float),
                   bem=np.array(bem), env=np.array(scan.envelope(bem)))
        np.savez(f"cube_cases_{case}{'_t_test' if 'CUBE_CASES_TEST' in os.environ else ''}.npz", **out)
        print(f"repetition {rep + 1} done", flush=True)
    if "CUBE_CASES_TEST" not in os.environ:
        report()
