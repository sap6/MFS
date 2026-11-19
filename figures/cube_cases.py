"""
Figure 7 and Table 5: the unit cube, error against wall time for the MFS (the layout of
the paper, two depths per octave, and one depth per octave) and for Bempp with the
same symmetry, in four cases: the standard wave at one wavelength
(k1, panel (d) of Figure 10), a general wave (k1o), and the standard wave at two
and three wavelengths (k2, k3). Errors are differences from the MFS reference
of each case.

    python figures/cube_cases.py                 # the figure, from data/cube_cases_*.npz
    python figures/cube_cases.py --scan CASE     # Bempp: scan of grading, quadrature and mesh size (hours)
    python figures/cube_cases.py --run CASE      # reference, then both codes timed alternately (hours)

For a fair comparison both codes run on the same machine with nothing else,
each run twice, alternately, and each keeps its smallest time.
"""
import gc
import os
import pickle
import sys
import time

import numpy as np
from style import BLUE, NAVY, RUST, data, envelope, finish, plt, rel

from mfs import core, cube

CASES = {"k1": (1, None), "k1o": (1, (60, 25, 35)), "k2": (2, None), "k3": (3, None)}   # wavelengths per side, wave
LEVELS = {"k1": range(4, 12), "k1o": range(4, 11), "k2": range(4, 10), "k3": range(4, 10)}   # two depths per octave
LEVELS_ONE = range(4, 10)                                # one depth per octave
J_REF = {"k1": 12, "k1o": 11, "k2": 10, "k3": 10}
EXTRA = {1: 0, 2: 4, 3: 12}                              # modes added to every line at two and three wavelengths
GS = (1.25, 1.5, 1.6, 1.75, 1.9, 2.0, 2.25)              # Bempp: grading exponents, quadrature orders and meshes
QUAD = ((4, 6, (12, 16, 20, 24, 28, 32, 40)), (6, 8, (12, 16, 20, 24, 28, 32, 40, 48)), (8, 10, (16, 20, 24, 28, 32, 40, 48)))


def setup(case):
    """Wavenumber, wave, test points, and the layout keywords of a level"""
    S, angles = CASES[case]
    d, p = core.wave(*angles) if angles else core.wave()
    layout = lambda J: dict(extra_modes=EXTRA[S])
    return 2 * np.pi * S, d, p, 2.0 * core.fibonacci(400), layout


def bempp_run(case, g, reg, sing, n, ref):
    """Bempp with the symmetry of the cube: one row (g, regular and singular order, n, unknowns, time, error)"""
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    bempp = bem.bempp
    k, d, p, pts, _ = setup(case)
    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = int(reg), int(sing)
    gc.collect()
    t0 = time.time()
    E, N = bem.cube(int(n), g, k, d, p, pts)
    return (g, reg, sing, n, N, time.time() - t0, rel(E, ref))


if "--scan" in sys.argv or "--run" in sys.argv:
    case = sys.argv[-1]
    k, d, p, pts, layout = setup(case)
    if "--scan" in sys.argv:
        # Every grading with every quadrature and mesh size, against the reference of the case
        ref = np.load(data(f"cube_cases_{case}_ref.npz"))["E"]
        rows = [bempp_run(case, g, reg, sing, n, ref) for reg, sing, ns in QUAD for g in GS for n in ns]
        np.savez(data(f"cube_scan_{case}.npz"), rows=np.array(rows))
    else:
        # The reference, then the MFS levels and the Bempp runs near the lower envelope of the scan, alternately
        t_start, max_s = time.time(), float(os.environ.get("MFS_MAX_SECONDS", 0))
        ck_path = data(f"cube_cases_{case}.ckpt.pkl")
        runs = [(2, J) for J in LEVELS[case]] + [(1, J) for J in LEVELS_ONE]   # (depths per octave, level)
        rows = [tuple(r) for r in np.load(data(f"cube_scan_{case}.npz"))["rows"]]
        cand = [r for r in rows if not any(q[5] <= 0.5 * r[5] and q[6] <= r[6] for q in rows)]
        st = pickle.load(open(ck_path, "rb")) if os.path.exists(ck_path) else {
            "E": None, "t_mfs": {r: [] for r in runs}, "t_bem": {r[:4]: [] for r in cand}, "res": {}}
        did, warm_m, warm_b = False, False, False
        if st["E"] is None:
            cube.solve(2, d, p, pts, k)                                     # warm-up: compilation is not timed
            warm_m = True
            gc.collect()
            E, t, info = cube.solve(J_REF[case], d, p, pts, k, **layout(J_REF[case]))
            print(f"{case} ref J={J_REF[case]}: {t:.1f} s, opt={info['optical']:.1e}", flush=True)
            info = [J_REF[case], t, max(info["unknowns"]), info["rows"], info["optical"], info["ext"]]
            np.savez(data(f"cube_cases_{case}_ref.npz"), E=E, pts=pts, info=info)
            st["E"] = E
            pickle.dump(st, open(ck_path, "wb"))
            did = True
        E, t_mfs, t_bem, res = st["E"], st["t_mfs"], st["t_bem"], st["res"]
        for rep in range(2):
            for nd, J in runs:
                if len(t_mfs[nd, J]) > rep:
                    continue
                if did and max_s and time.time() - t_start > (20 if (nd == 2 and J >= 9) else max_s):
                    pickle.dump(st, open(ck_path, "wb"))
                    sys.exit(2)
                if not warm_m:
                    cube.solve(2, d, p, pts, k)
                    warm_m = True
                gc.collect()
                EJ, t, info = cube.solve(J, d, p, pts, k, q=2 ** -0.5 if nd == 2 else 0.5, **layout(J))
                t_mfs[nd, J].append(t)
                res[nd, J] = (max(info["unknowns"]), sum(info["unknowns"]), info["rows"], rel(EJ, E), info["optical"])
                pickle.dump(st, open(ck_path, "wb"))
                did = True
                print(f"{case} pass {rep + 1} MFS ({nd}, {J}): {t:.2f} s, err={res[nd, J][3]:.2e}", flush=True)
            for r in cand:
                if len(t_bem[r[:4]]) > rep:
                    continue
                if did and max_s and time.time() - t_start > (20 if r[3] >= 48 else max_s):
                    pickle.dump(st, open(ck_path, "wb"))
                    sys.exit(2)
                if not warm_b:
                    bempp_run(case, 1.0, 4, 6, 4, E)
                    warm_b = True
                b = bempp_run(case, *r[:4], E)
                t_bem[r[:4]].append(b[5])
                res[r[:4]] = b
                pickle.dump(st, open(ck_path, "wb"))
                did = True
                print(f"{case} pass {rep + 1} Bempp {r[:4]}: {b[5]:.2f} s, err={b[6]:.2e}", flush=True)
        mfs = np.array([(nd, J, *res[nd, J][:3], min(t_mfs[nd, J]), *res[nd, J][3:]) for nd, J in runs])
        bem = np.array([(*res[key][:5], min(t_bem[key]), res[key][6]) for key in t_bem])
        best = envelope([(r[5], r[6]) for r in bem])                     # (time, error) of the lower envelope
        env = np.array([r for r in bem if (r[5], r[6]) in set(map(tuple, best))])
        np.savez(data(f"cube_cases_{case}.npz"), mfs=mfs, bem=bem, env=env[np.argsort(env[:, 5])])
        os.remove(ck_path)

# The figure: three panels (the standard wave is panel (d) of Figure 10)
plt.rcParams["font.size"] = 8.5
fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9), sharey=True)
panels = (("k1o", "(a) one wavelength, general wave"), ("k2", "(b) two wavelengths"), ("k3", "(c) three wavelengths"))
for ax, (case, title) in zip(axs, panels):
    res = np.load(data(f"cube_cases_{case}.npz"))
    for nd, col, mk, ls, label in ((2, BLUE, "o", "-", "MFS"), (1, NAVY, "h", (0, (6, 1.5)), "MFS, one depth per octave")):
        m = res["mfs"][res["mfs"][:, 0] == nd]
        if len(m):
            ax.loglog(m[:, 5], m[:, 6], linestyle=ls, color=col, marker=mk, ms=5, lw=1.5, mec="white", mew=0.7, label=label)
    env = res["env"][res["env"][:, 6] < 2e-2]
    ax.loglog(env[:, 5], env[:, 6], linestyle=(0, (3, 1)), color=RUST, marker="X", ms=5, lw=1.5, mec="white", mew=0.7,
              label="Bempp, graded, with symmetry")
    ax.set_title(title, loc="left", fontsize=8.5)
    ax.set_xlabel("wall time [s]")
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
axs[0].set_ylabel("relative error of the scattered field")
axs[0].legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=2.6)
finish(fig, axs, "fig_cube_cases")
