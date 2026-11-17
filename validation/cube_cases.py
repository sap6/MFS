"""The cube beyond the standard case: a general incident wave, and a side of
two and three wavelengths. For each case the MFS levels and the best Bempp
runs (with the mirror symmetry; candidates from bem_grading_scan.py) are run
alternately, several times, and each keeps its smallest wall time, as
cube_duel.py does for the standard case (NOTES.md items 60, 64).

    python3 cube_cases.py k1o [repetitions]   # k = 2 pi, oblique wave (eight parity classes)
    python3 cube_cases.py k2  [repetitions]   # k = 4 pi, standard wave
    python3 cube_cases.py k2o [repetitions]   # k = 4 pi, oblique wave
    python3 cube_cases.py k3  [repetitions]   # k = 6 pi, standard wave
    python3 cube_cases.py k3 mfs              # MFS only, once (no Bempp scan needed)
    python3 cube_cases.py report              # tables of all finished cases

Errors are differences from the reference of bem_grading_scan.reference(case).
Writes cube_cases_<case>.npz: mfs (family 0 = one depth per octave, 1 = two;
J; unknowns per class; unknowns of all classes; equations per class; time;
error; optical mismatch), bem (all timed Bempp runs) and env (lower envelope),
both with the columns of bem_grading_scan.py.
"""
import os
import sys

import numpy as np

import pec_cube_lightning3d as m

SIDE = {"k1o": 1, "k2": 2, "k2o": 2, "k3": 3}
OBLIQUE = (60.0, 25.0, 35.0)
LEVELS = {                                                  # (levels with one depth per octave, with two); the two-depth
    "k1o": (range(4, 10), range(4, 9)),                     # levels stop two below the reference level, whose own error
    "k2": (range(4, 9), range(5, 9)),                       # would otherwise bias the last entry (NOTES.md item 64)
    "k2o": (range(5, 9), (7, 8)),
    "k3": (range(5, 9), range(6, 9)),
}


def layout(c, S, J, family):
    """How the layout of cube_mfs.py follows the wavenumber (NOTES.md items 57, 64).
    S = 1: the family as it is. S = 2: interior collocation panels scaled with
    the wavelength, 4 more modes per line. S = 3: in addition the first depth
    is 0.2 / S and the base grid has 16 points per face."""
    if S == 1:
        return {}
    kw = c.dense(J) if family == "dense" else dict(zip(("kw", "hmax"), c.balanced(J)))
    kw = dict(kw)
    n_base, hmax = (4, 0.25 / S) if J <= 5 else (6, 0.15 / S)
    depth0 = 0.2
    if S > 2:
        n_base, hmax, depth0 = BASE3, HMAX3, 0.2 / S
    kw["kw"] = dict(kw["kw"], n_base=n_base, n0=kw["kw"]["n0"] + int(round(4 * (S - 1))), depth0=depth0)
    kw["hmax"] = hmax
    return kw


BASE3, HMAX3 = 18, 0.05                                     # three wavelengths: base points per face, interior panel length


def report():
    import bem_grading_scan as scan
    for case in SIDE:
        f = f"cube_cases_{case}.npz"
        if not os.path.exists(f):
            continue
        d = np.load(f)
        print(f"{case}: MFS  family  J  unknowns/class  all  equations  time  error  optical")
        for r in d["mfs"]:
            print(f"   {'two' if r[0] else 'one'} {int(r[1]):2d} {int(r[2]):6d} {int(r[3]):6d} {int(r[4]):6d} {r[5]:7.1f} s  {r[6]:.2e}  {r[7]:.1e}")
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
    m.K = 2 * np.pi * S                                     # before any solver object is built
    import cube_mfs as c
    import bem_grading_scan as scan
    mfs_only = len(sys.argv) > 2 and sys.argv[2] == "mfs"
    reps = 1 if mfs_only else (int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    k, ref, wave = scan.reference(case)
    pts = np.load("bem_cube_h0.05.npz")["pts"]
    dp = c.wave(*OBLIQUE) if wave is not None else (m.D, m.P)
    runs = [("balanced", J) for J in LEVELS[case][0]] + [("dense", J) for J in LEVELS[case][1]]
    for fam in ("balanced", "dense"):                       # compilation
        c.run(2, dp[0], dp[1], pts, family=fam, **layout(c, S, 2, fam))
    cand = []
    if not mfs_only:
        import bem_cube_sym as b
        import bempp_cl.api as bempp
        b.K = k
        rows = [tuple(r) for r in np.load(f"bem_grading_scan_{case}.npz")["rows"]]
        cand = [r for r in rows if not any(q[5] <= 0.75 * r[5] and q[6] <= r[6] for q in rows)]
        cand.sort(key=lambda r: (r[1], r[5]))
    print(f"{case}: {len(runs)} MFS runs and {len(cand)} Bempp runs, {reps} repetitions", flush=True)
    t_mfs, t_bem, res = {r: [] for r in runs}, {r: [] for r in cand}, {}
    args = (pts,) if wave is None else (pts, dp[0], dp[1])
    for rep in range(reps):
        for fam, J in runs:
            E, t, info = c.run(J, dp[0], dp[1], pts, family=fam, **layout(c, S, J, fam))
            t_mfs[(fam, J)].append(t)
            res[(fam, J)] = (max(info["per_class"]), info["unknowns"], info["rows"], np.linalg.norm(E - ref) / np.linalg.norm(ref), info["optical"])
            if rep == 0:
                print(f"  {fam:8s} J = {J:2d}: unknowns/class {res[(fam, J)][0]} rows {info['rows']} t = {t:6.1f} s  error {res[(fam, J)][3]:.2e} "
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
        out = dict(mfs=np.array([(fam == "dense", J, res[(fam, J)][0], res[(fam, J)][1], res[(fam, J)][2], min(t_mfs[(fam, J)]),
                                  res[(fam, J)][3], res[(fam, J)][4]) for fam, J in runs], float))
        if cand:
            bem = [(r[0], r[1], r[2], r[3], r[4], min(t_bem[r]), r[6]) for r in cand]
            out.update(bem=np.array(bem), env=np.array(scan.envelope(bem)))
        np.savez(f"cube_cases_{case}{'_mfs' if mfs_only else ''}.npz", **out)
        print(f"repetition {rep + 1} done", flush=True)
    if not mfs_only:
        report()
