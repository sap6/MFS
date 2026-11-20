"""
Figure 10 and Table 4: error of the scattered field against wall time on two
cores, for the MFS and for Bempp: (a) the unit sphere, against the Mie series;
(b) the torus and (c) the almond, against the finest MFS run; (d) the cube
(data of cube_convergence.py and cube_cases.py k1). Both codes use the same
symmetry: none on the sphere, the rotation on the torus, the two mirror planes
on the almond, and the mirrors and the rotation on the cube.

    python figures/cost.py                  # the figure, from data/
    python figures/cost.py --run SHAPE      # time both codes again: sphere, torus or almond (needs bempp-cl)

Run alone on the machine: times include assembly, solution and the field at the test points.
"""
import gc
import sys
import time

import numpy as np
from style import BLUE, NAVY, ORANGE, RUST, curve, data, finish, plt, rel

from mfs import almond, core, mie, sphere, torus

K = 2 * np.pi


def timed(f, *args):
    ts = []
    for _ in range(2):
        gc.collect()
        t0 = time.time()
        out = f(*args)
        ts.append(time.time() - t0)
    return out, min(ts)


if "--run" in sys.argv:
    shape = sys.argv[-1]
    if shape not in ("sphere", "torus", "almond"):
        sys.exit("usage: python figures/cost.py --run sphere|torus|almond")
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    bempp = bem.bempp
    d, p = core.wave()
    if shape == "sphere":
        # 400 points on the sphere of radius 2, against the Mie series
        pts = 2.0 * core.fibonacci(400)
        exact = mie.scattered_field(pts, K, d, p)
        mfs = lambda N: core.field(pts, *sphere.solve(N), K)
        mfs(50), bem.efie(bempp.shapes.sphere(h=0.4), K, d, p, pts)          # warm-up: compilation is not timed
        rows = [(3 * N, t, rel(E, exact)) for N in (50, 100, 150, 200, 300, 400, 600, 800, 1200) for E, t in [timed(mfs, N)]]
        bems = [(N, t, rel(E, exact)) for h in (0.2, 0.14, 0.1, 0.07, 0.056)
                for (E, N), t in [timed(bem.efie, bempp.shapes.sphere(h=h), K, d, p, pts)]]
        np.savez(data("cost_sphere.npz"), mfs=np.array(rows), bem=np.array(bems))
    elif shape == "torus":
        # 200 points on the sphere of radius 2.5, against the finest MFS run (n_s = 56 rings of 240 dipoles)
        pts = 2.5 * core.sphere_rule(10, 20)[0]
        ref = torus.scatter(K, d, p, 56, 240)
        bem.torus(6, K, d, p, pts)                                          # warm-up
        Eref, rows = ref.field(pts), []
        def run(ns, nphi):
            sol = torus.scatter(K, d, p, ns, nphi)
            return sol, sol.field(pts)
        for ns, nphi in ((6, 20), (8, 28), (12, 40), (16, 52), (20, 64), (24, 80), (28, 96), (36, 128), (44, 160)):
            (sol, E), t = timed(run, ns, nphi)                              # the checks below are not timed
            checks = (torus.residual(sol, d, p), torus.optical(sol, d, p)[1], np.linalg.norm(sol.alpha))
            rows.append((3 * ns * nphi, t, rel(E, Eref), *checks))
        bems = [(N, t, rel(E, Eref)) for n in (12, 20, 28, 40, 56, 80, 112, 160, 224)
                for (E, N), t in [timed(bem.torus, n, K, d, p, pts)]]
        ext, opt = torus.optical(ref, d, p)
        np.savez(data("cost_torus.npz"), ref=[3 * 56 * 240, torus.residual(ref, d, p), opt, ext],
                 mfs=np.array(rows), bem=np.array(bems))
    elif shape == "almond":
        # Nose-on at 1.19 GHz, 400 points on the sphere of radius 0.3 m; the layout 8:5 is the reference
        k = 2 * np.pi * 1.19e9 / almond.C0
        pts = 0.3 * core.fibonacci(400)
        da, pa, out, rows, brows = np.array([-1.0, 0, 0]), np.array([0, 1.0, 0]), {}, [], []
        almond.nose_on(k, 2, 0), bem.almond(12, 12, k, da, pa, pts)        # warm-up
        for i, (n_mid, n_F) in enumerate(((3, 0), (4, 1), (5, 2), (6, 3), (8, 5))):
            (src, c, E, shape_), t = timed(almond.nose_on, k, n_mid, n_F, pts)
            ext, sca = almond.optical(k, src, c)                            # the check is not timed, as for Bempp
            rows.append((len(src[0]), shape_[1], t, abs(ext - sca) / sca))
            out[f"mfs_E{i}"] = E
        for i, (n_t, n_s) in enumerate(((28, 28), (40, 36), (60, 56), (84, 76), (100, 92))):
            (E, N), t = timed(bem.almond, n_t, n_s, k, da, pa, pts)
            brows.append((N, t))
            out[f"bem_E{i}"] = E
        np.savez(data("cost_almond.npz"), mfs=np.array(rows), bem=np.array(brows), **out)

# The figure
fig, axs = plt.subplots(2, 2, figsize=(6.5, 4.8))
ax = axs[0, 0]
res = np.load(data("cost_sphere.npz"))
curve(ax, res["mfs"][:, 1], res["mfs"][:, 2], BLUE, "o", "-", "MFS")
curve(ax, res["bem"][:, 1], res["bem"][:, 2], ORANGE, "D", "-.", "Bempp")
ax.set_title(r"(a) unit sphere, $k = 2\pi$ (against Mie)", loc="left", fontsize=8.5)

ax = axs[0, 1]
res = np.load(data("cost_torus.npz"))
curve(ax, res["mfs"][:, 1], res["mfs"][:, 2], BLUE, "o", "-", "MFS")
curve(ax, res["bem"][:, 1], res["bem"][:, 2], ORANGE, "D", "-.", "Bempp, with symmetry")
ax.set_title(r"(b) torus, $R = 1$, $a = 0.5$, $k = 2\pi$", loc="left", fontsize=8.5)
ax.set_xlim(right=250)

ax = axs[1, 0]
res = np.load(data("cost_almond.npz"))
ref = res[f"mfs_E{len(res['mfs']) - 1}"]
curve(ax, res["mfs"][:-1, 2], [rel(res[f"mfs_E{i}"], ref) for i in range(len(res["mfs"]) - 1)], BLUE, "o", "-", "MFS")
curve(ax, res["bem"][:, 1], [rel(res[f"bem_E{i}"], ref) for i in range(len(res["bem"]))], ORANGE, "D", "-.", "Bempp, with symmetry")
ax.set_title(r"(c) almond, $1.19$ GHz, nose-on", loc="left", fontsize=8.5)
ax.set_ylim(3e-5, 1.2e-1)
ax.set_xticks([1, 2, 5, 10, 20, 50])
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:g}"))

ax = axs[1, 1]
res = np.load(data("cube_convergence.npz"))
ref, jref = res["ref_E"], int(res["ref_info"][0])
for key, col, mk, ls, lab in (("dense", BLUE, "o", "-", "MFS (two depths per octave)"),
                              ("tail", NAVY, "h", (0, (6, 1.5)), "MFS (one depth per octave)")):
    rows = [r for r in res[key] if 4 <= r[0] < jref]
    curve(ax, [r[1] for r in rows], [rel(res[f"{key}_E{int(r[0])}"], ref) for r in rows], col, mk, ls, lab)
u = np.load(data("cube_bempp_uniform.npz"))
keep = u["ndof"] > 1000
curve(ax, u["t"][keep], [rel(E, ref) for E in u["E"][keep]], ORANGE, "D", "-.", "Bempp, quasi-uniform")
env = np.load(data("cube_cases_k1.npz"))["env"]
curve(ax, env[:, 5], env[:, 6], RUST, "X", (0, (3, 1)), "Bempp, graded, with symmetry")
ax.set_title(r"(d) unit cube, $k = 2\pi$, standard wave", loc="left", fontsize=8.5)
ax.set_ylim(1.5e-7, 3e-2)

locs = {(0, 0): ("lower left", None), (0, 1): ("lower left", None),
        (1, 0): ("upper right", None), (1, 1): ("lower left", (0.0, -0.01))}
for idx, ax in np.ndenumerate(axs):
    loc, anchor = locs[idx]
    ax.legend(frameon=False, fontsize=7.5, loc=loc, bbox_to_anchor=anchor,
              handlelength=1.9, handletextpad=0.45, labelspacing=0.35)
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
for ax in axs[1]:
    ax.set_xlabel("wall time [s]  (2 cores)")
for ax in axs[:, 0]:
    ax.set_ylabel("relative error of the scattered field")
finish(fig, axs, "fig_cost_all")
