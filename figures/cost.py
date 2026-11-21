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


def alternate(mfs_jobs, bem_jobs):
    """Run MFS and Bempp jobs twice in alternating passes, keeping the minimum wall time of each job"""
    tm, tb = [[] for _ in mfs_jobs], [[] for _ in bem_jobs]
    om, ob = [None] * len(mfs_jobs), [None] * len(bem_jobs)
    for rep in range(2):
        for i, fn in enumerate(mfs_jobs):
            gc.collect()
            t0 = time.time()
            om[i] = fn()
            tm[i].append(time.time() - t0)
        for i, fn in enumerate(bem_jobs):
            gc.collect()
            t0 = time.time()
            ob[i] = fn()
            tb[i].append(time.time() - t0)
        print(f"completed alternating pass {rep + 1}/2", flush=True)
    return [(o, min(t)) for o, t in zip(om, tm)], [(o, min(t)) for o, t in zip(ob, tb)]


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
        Ns = (50, 100, 150, 200, 300, 400, 600, 800, 1200)
        grids = [bempp.shapes.sphere(h=h) for h in (0.2, 0.14, 0.1, 0.07, 0.056)]
        m_res, b_res = alternate([(lambda N=N: mfs(N)) for N in Ns],
                                 [(lambda g=g: bem.efie(g, K, d, p, pts)) for g in grids])
        rows = [(3 * N, t, rel(E, exact)) for N, (E, t) in zip(Ns, m_res)]
        bems = [(N, t, rel(E, exact)) for (E, N), t in b_res]
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
        m_specs = ((6, 20), (8, 28), (12, 40), (16, 52), (20, 64), (24, 80), (28, 96), (36, 128), (44, 160))
        b_specs = (12, 20, 28, 40, 56, 80, 112, 160, 224)
        m_res, b_res = alternate([(lambda s=s: run(*s)) for s in m_specs],
                                 [(lambda n=n: bem.torus(n, K, d, p, pts)) for n in b_specs])
        for (ns, nphi), ((sol, E), t) in zip(m_specs, m_res):               # the checks below are not timed
            checks = (torus.residual(sol, d, p), torus.optical(sol, d, p)[1], np.linalg.norm(sol.alpha))
            rows.append((3 * ns * nphi, t, rel(E, Eref), *checks))
        bems = [(N, t, rel(E, Eref)) for (E, N), t in b_res]
        ext, opt = torus.optical(ref, d, p)
        np.savez(data("cost_torus.npz"), ref=[3 * 56 * 240, torus.residual(ref, d, p), opt, ext],
                 mfs=np.array(rows), bem=np.array(bems))
    elif shape == "almond":
        # Nose-on at 1.19 GHz, 400 points on the sphere of radius 0.3 m; the layout 8:5 is the reference
        k = 2 * np.pi * 1.19e9 / almond.C0
        pts = 0.3 * core.fibonacci(400)
        da, pa, out, rows, brows = np.array([-1.0, 0, 0]), np.array([0, 1.0, 0]), {}, [], []
        almond.nose_on(k, 2, 0, pts), bem.almond(12, 12, k, da, pa, pts)    # warm-up
        m_specs = ((3, 0), (4, 1), (5, 2), (6, 3), (8, 5))
        b_specs = ((28, 28), (40, 36), (60, 56), (84, 76), (100, 92))
        m_res, b_res = alternate([(lambda s=s: almond.nose_on(k, s[0], s[1], pts)) for s in m_specs],
                                 [(lambda s=s: bem.almond(s[0], s[1], k, da, pa, pts)) for s in b_specs])
        for i, ((src, c, E, shape_), t) in enumerate(m_res):
            ext, sca = almond.optical(k, src, c)                            # the check is not timed, as for Bempp
            rows.append((len(src[0]), shape_[1], t, abs(ext - sca) / sca))
            out[f"mfs_E{i}"] = E
        for i, ((E, N), t) in enumerate(b_res):
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
