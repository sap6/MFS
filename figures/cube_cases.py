"""
Figure 7 and Table 5: the unit cube, error against wall time for the MFS and for
Bempp with the same symmetry, in four cases: the standard wave at one wavelength
(k1, panel (d) of Figure 10), a general wave (k1o), and the standard wave at two
and three wavelengths (k2, k3). Errors are differences from the MFS reference
of each case.

    python figures/cube_cases.py                 # the figure, from data/cube_cases_*.npz
    python figures/cube_cases.py --scan CASE     # Bempp: scan of grading, quadrature and mesh size (hours)
    python figures/cube_cases.py --run CASE      # reference, then both codes timed alternately (hours)

For a fair comparison both codes run on the same machine with nothing else,
each run twice, alternately, and each keeps its smallest time.
"""
import sys
import time

import numpy as np
from style import BLUE, RUST, data, envelope, finish, plt, rel

from mfs import core, cube

CASES = {"k1": (1, None), "k1o": (1, (60, 25, 35)), "k2": (2, None), "k3": (3, None)}   # wavelengths per side, wave
LEVELS = {"k1": range(4, 12), "k1o": range(4, 11), "k2": range(4, 10), "k3": range(4, 10)}
J_REF = {"k1": 12, "k1o": 11, "k2": 10, "k3": 10}
EXTRA = {1: 0, 2: 4, 3: 12}                              # modes added to every line at two and three wavelengths
GS = (1.25, 1.5, 1.6, 1.75, 1.9, 2.0, 2.25)              # Bempp: grading exponents, quadrature orders and meshes
QUAD = ((4, 6, (12, 16, 20, 24, 28, 32, 40)), (6, 8, (12, 16, 20, 24, 28, 32, 40, 48)), (8, 10, (16, 20, 24, 28, 32, 40, 48)))


def setup(case):
    """Wavenumber, wave, test points, and the layout keywords of a level"""
    S, angles = CASES[case]
    d, p = core.wave(*angles) if angles else core.wave()
    layout = lambda J: dict(extra_modes=EXTRA[S], max_gb=2.5 if (S == 3 and J >= 9) else (8.5 if case == "k1o" else 6.0))
    return 2 * np.pi * S, d, p, 2.0 * core.fibonacci(400), layout


def bempp_run(case, g, reg, sing, n, ref):
    import bempp_cl.api as bempp
    from mfs import bem
    k, d, p, pts, _ = setup(case)
    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = reg, sing
    t0 = time.time()
    E, N = bem.cube(n, g, k, d, p, pts)
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
        E, t, info = cube.solve(J_REF[case], d, p, pts, k, **layout(J_REF[case]))
        np.savez(data(f"cube_cases_{case}_ref.npz"), E=E, pts=pts, info=[J_REF[case], t, max(info["unknowns"]), info["rows"], info["optical"], info["ext"]])
        rows = [tuple(r) for r in np.load(data(f"cube_scan_{case}.npz"))["rows"]]
        cand = [r for r in rows if not any(q[5] <= 0.5 * r[5] and q[6] <= r[6] for q in rows)]
        t_mfs, t_bem, res = {J: [] for J in LEVELS[case]}, {r[:4]: [] for r in cand}, {}
        for rep in range(2):
            for J in LEVELS[case]:
                EJ, t, info = cube.solve(J, d, p, pts, k, **layout(J))
                t_mfs[J].append(t)
                res[J] = (max(info["unknowns"]), sum(info["unknowns"]), info["rows"], rel(EJ, E), info["optical"])
            for r in cand:
                b = bempp_run(case, *r[:4], E)
                t_bem[r[:4]].append(b[5])
                res[r[:4]] = b
        mfs = np.array([(1, J, *res[J][:3], min(t_mfs[J]), *res[J][3:]) for J in LEVELS[case]])
        bem = np.array([(*res[key][:5], min(t_bem[key]), res[key][6]) for key in t_bem])
        best = envelope([(r[5], r[6]) for r in bem])                     # (time, error) of the lower envelope
        env = np.array([r for r in bem if (r[5], r[6]) in set(map(tuple, best))])
        np.savez(data(f"cube_cases_{case}.npz"), mfs=mfs, bem=bem, env=env[np.argsort(env[:, 5])])

# The figure: three panels (the standard wave is panel (d) of Figure 10)
plt.rcParams["font.size"] = 8.5
fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9), sharey=True)
for ax, (case, title) in zip(axs, (("k1o", "(a) one wavelength, general wave"), ("k2", "(b) two wavelengths"), ("k3", "(c) three wavelengths"))):
    d = np.load(data(f"cube_cases_{case}.npz"))
    ax.loglog(d["mfs"][:, 5], d["mfs"][:, 6], "-", color=BLUE, marker="o", ms=5, lw=1.5, mec="white", mew=0.7, label="MFS")
    env = d["env"][d["env"][:, 6] < 2e-2]
    ax.loglog(env[:, 5], env[:, 6], linestyle=(0, (3, 1)), color=RUST, marker="X", ms=5, lw=1.5, mec="white", mew=0.7,
              label="Bempp, graded, with symmetry")
    ax.set_title(title, loc="left", fontsize=8.5)
    ax.set_xlabel("wall time [s]")
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
axs[0].set_ylabel("relative error of the scattered field")
axs[0].legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=2.6)
finish(fig, axs, "fig_cube_cases")
