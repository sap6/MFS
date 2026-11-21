"""
Figure 7 and Table 5: the unit cube, error against wall time on two cores for the MFS
(the layout of the paper, two depths per octave, and one depth per octave) and for Bempp
with the same symmetry, in four cases: the standard wave at one wavelength
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
import sys
import time

import numpy as np
from style import BLUE, NAVY, RUST, curve, data, envelope, finish, plt, rel

from mfs import core, cube

CASES = {"k1": (1, None), "k1o": (1, (60, 25, 35)), "k2": (2, None), "k3": (3, None)}   # wavelengths per side, wave
LEVELS = {"k1": range(4, 12), "k1o": range(4, 11), "k2": range(4, 10), "k3": range(4, 10)}   # two depths per octave
LEVELS_ONE = range(4, 10)                                # one depth per octave
J_REF = {"k1": 12, "k1o": 11, "k2": 10, "k3": 10}
GS = (1.25, 1.5, 1.6, 1.75, 1.9, 2.0, 2.25)              # Bempp: grading exponents, quadrature orders and meshes
QUAD = ((4, 6, (12, 16, 20, 24, 28, 32, 40, 48)), (6, 8, (12, 16, 20, 24, 28, 32, 40, 48)),
        (8, 10, (16, 20, 24, 28, 32, 40, 48)))


def setup(case):
    """Wavenumber, wave, and test points"""
    S, angles = CASES[case]
    d, p = core.wave(*angles) if angles else core.wave()
    return 2 * np.pi * S, d, p, 2.0 * core.fibonacci(400)


def bempp_run(case, g, reg, sing, n, ref):
    """Bempp with the symmetry of the cube: one row (g, regular and singular order, n, unknowns, time, error)"""
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    bempp = bem.bempp
    k, d, p, pts = setup(case)
    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = int(reg), int(sing)
    gc.collect()
    t0 = time.time()
    E, N = bem.cube(int(n), g, k, d, p, pts)
    return (g, reg, sing, n, N, time.time() - t0, rel(E, ref))


if "--scan" in sys.argv or "--run" in sys.argv:
    case = sys.argv[-1]
    if case not in CASES:
        sys.exit(f"usage: python figures/cube_cases.py [--scan|--run] {'|'.join(CASES)}")
    k, d, p, pts = setup(case)
    if "--scan" in sys.argv:
        # Every grading with every quadrature and mesh size (1.5 to 1.9 on the finest mesh), against the reference
        ref = np.load(data(f"cube_cases_{case}_ref.npz"))["E"]
        rows = [bempp_run(case, g, reg, sing, n, ref) for reg, sing, ns in QUAD for n in ns for g in (GS[1:5] if n == 48 else GS)]
        np.savez(data(f"cube_scan_{case}.npz"), rows=np.array(rows))
    else:
        # The reference, then the MFS levels and the best-graded Bempp run of each (quadrature, mesh), alternately
        cube.solve(2, d, p, pts, k)                                     # warm-up: compilation is not timed
        E, t, info = cube.solve(J_REF[case], d, p, pts, k)
        print(f"{case} ref J={J_REF[case]}: {t:.1f} s, opt={info['optical']:.1e}", flush=True)
        info = [J_REF[case], t, max(info["unknowns"]), info["rows"], info["optical"], info["ext"]]
        np.savez(data(f"cube_cases_{case}_ref.npz"), E=E, pts=pts, info=info)
        bempp_run(case, 1.0, 4, 6, 4, E)
        runs = [(2, J) for J in LEVELS[case]] + [(1, J) for J in LEVELS_ONE]   # (depths per octave, level)
        rows = [tuple(r) for r in np.load(data(f"cube_scan_{case}.npz"))["rows"]]
        cand = [min((r for r in rows if r[1:4] == qm), key=lambda r: r[6]) for qm in dict.fromkeys(r[1:4] for r in rows)]
        t_mfs, t_bem, res = {r: [] for r in runs}, {r[:4]: [] for r in cand}, {}
        for rep in range(2):
            for nd, J in runs:
                gc.collect()
                EJ, t, info = cube.solve(J, d, p, pts, k, q=2 ** -0.5 if nd == 2 else 0.5)
                t_mfs[nd, J].append(t)
                res[nd, J] = (max(info["unknowns"]), sum(info["unknowns"]), info["rows"], rel(EJ, E), info["optical"])
                print(f"{case} pass {rep + 1} MFS ({nd}, {J}): {t:.2f} s, err={res[nd, J][3]:.2e}", flush=True)
            for r in cand:
                b = bempp_run(case, *r[:4], E)
                t_bem[r[:4]].append(b[5])
                res[r[:4]] = b
                print(f"{case} pass {rep + 1} Bempp {r[:4]}: {b[5]:.2f} s, err={b[6]:.2e}", flush=True)
        mfs = np.array([(nd, J, *res[nd, J][:3], min(t_mfs[nd, J]), *res[nd, J][3:]) for nd, J in runs])
        bem = np.array([(*res[key][:5], min(t_bem[key]), res[key][6]) for key in t_bem])
        best = envelope([(r[5], r[6]) for r in bem])                     # (time, error) of the lower envelope
        env = np.array([r for r in bem if (r[5], r[6]) in set(map(tuple, best))])
        np.savez(data(f"cube_cases_{case}.npz"), mfs=mfs, bem=bem, env=env[np.argsort(env[:, 5])])

# The figure: three panels (the standard wave is panel (d) of Figure 10)
fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.9), sharex=True, sharey=True)
panels = (("k1o", "(a) $k = 2\\pi$, general wave"), ("k2", "(b) $k = 4\\pi$, two wavelengths"), ("k3", "(c) $k = 6\\pi$, three wavelengths"))
for ax, (case, title) in zip(axs, panels):
    res = np.load(data(f"cube_cases_{case}.npz"))
    mfs = res["mfs"]
    m2 = mfs[mfs[:, 0] == 2]
    for pts, col, mk, ls, zo, label in (
        (m2, BLUE, "o", "-", 4, "MFS (two depths per octave)"),
        (mfs[mfs[:, 0] == 1], NAVY, "h", (0, (6, 1.5)), 3, "MFS (one depth per octave)"),
        (res["env"][res["env"][:, 6] < 2.5e-2], RUST, "X", (0, (3, 1)), 2, "Bempp, graded, with symmetry"),
    ):
        curve(ax, pts[:, 5], pts[:, 6], col, mk, ls, label, 4.8, 1.4, 0.6, zorder=zo)
    ax.text(m2[0, 5] * 1.18, m2[0, 6] * 1.12, f"$J={int(m2[0, 1])}$", fontsize=7.0, color=BLUE)
    ax.text(m2[-1, 5] * 1.18, m2[-1, 6] * 0.92, f"$J={int(m2[-1, 1])}$", fontsize=7.0, color=BLUE, va="center")
    ax.set_title(title, loc="left", fontsize=8.5)
    ax.set_xlabel("wall time [s]  (2 cores)")
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
axs[0].set(xlim=(0.8, 400), ylim=(1e-5, 3.5e-2), ylabel="relative error of the scattered field")
h, lab = axs[0].get_legend_handles_labels()
fig.legend(h, lab, frameon=False, fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.005))
finish(fig, axs, "fig_cube_cases", rect=(0, 0.09, 1, 1))
