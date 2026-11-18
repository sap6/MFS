"""
Figure 10 and Table 4: error of the scattered field against wall time on four
cores, for the MFS and for Bempp: (a) the unit sphere, against the Mie series;
(b) the torus and (c) the almond, against the finest MFS run; (d) the cube
(data of cube_convergence.py and cube_cases.py k1). Both codes use the same
symmetry: none on the sphere, the rotation on the torus, the two mirror planes
on the almond, and the mirrors and the rotation on the cube.

    python figures/cost.py                  # the figure, from data/
    python figures/cost.py --run SHAPE      # time both codes again: sphere, torus or almond (needs bempp-cl)

Run alone on the machine: times include assembly, solution and the field at the test points.
"""
import sys
import time

import matplotlib.ticker
import numpy as np
from style import BLUE, NAVY, ORANGE, RUST, curve, data, finish, plt, rel

from mfs import almond, core, mie, sphere, torus

K = 2 * np.pi


def timed(f, *args):
    t0 = time.time()
    out = f(*args)
    return out, time.time() - t0


if "--run" in sys.argv:
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    bempp = bem.bempp
    shape, (d, p) = sys.argv[-1], core.wave()
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
        Eref = ref.field(pts)
        rows = []
        for ns, nphi in ((6, 20), (8, 28), (12, 40), (16, 52), (20, 64), (24, 80), (28, 96), (36, 128), (44, 160)):
            t0 = time.time()
            sol = torus.scatter(K, d, p, ns, nphi)
            E, t = sol.field(pts), time.time() - t0                  # the checks below are not timed
            rows.append((3 * ns * nphi, t, rel(E, Eref), torus.residual(sol, d, p), torus.optical(sol, d, p)[1], np.linalg.norm(sol.alpha)))
        bems = [(N, t, rel(E, Eref)) for n in (12, 20, 28, 40, 56, 80, 112, 160, 224)
                for (E, N, _), t in [timed(bem.torus, n, K, d, p, pts)]]
        ext, opt = torus.optical(ref, d, p)
        np.savez(data("cost_torus.npz"), ref=[3 * 56 * 240, torus.residual(ref, d, p), opt, ext], mfs=np.array(rows), bem=np.array(bems))
    elif shape == "almond":
        # Nose-on at 1.19 GHz, 400 points on the sphere of radius 0.3 m; the layout 8:5 is the reference
        k = 2 * np.pi * 1.19e9 / almond.C0
        pts = 0.3 * core.fibonacci(400)
        out, rows = {}, []
        almond.nose_on(k, 2, 0), bem.almond(12, 12, k, np.array([-1.0, 0, 0]), np.array([0, 1.0, 0]), pts)   # warm-up
        for i, (n_mid, n_F) in enumerate(((3, 0), (4, 1), (5, 2), (6, 3), (8, 5))):
            (y, c, E, (ext, sca), shape_), t = timed(almond.nose_on, k, n_mid, n_F, pts)
            rows.append((len(y), shape_[1], t, abs(ext - sca) / sca))
            out[f"mfs_E{i}"] = E
        brows = []
        for i, (n_t, n_s) in enumerate(((28, 28), (40, 36), (60, 56), (84, 76), (100, 92))):
            (E, N), t = timed(bem.almond, n_t, n_s, k, np.array([-1.0, 0, 0]), np.array([0, 1.0, 0]), pts)
            brows.append((N, t))
            out[f"bem_E{i}"] = E
        np.savez(data("cost_almond.npz"), mfs=np.array(rows), bem=np.array(brows), **out)

# The figure
plt.rcParams["font.size"] = 8.5
fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.4))
ax = axs[0, 0]
d = np.load(data("cost_sphere.npz"))
curve(ax, d["mfs"][:, 1], d["mfs"][:, 2], BLUE, "o", "-", "MFS")
curve(ax, d["bem"][:, 1], d["bem"][:, 2], ORANGE, "D", "-.", "Bempp")
ax.set_title("(a) sphere, two wavelengths across (against Mie)", loc="left", fontsize=8.5)

ax = axs[0, 1]
d = np.load(data("cost_torus.npz"))
curve(ax, d["mfs"][:, 1], d["mfs"][:, 2], BLUE, "o", "-", "MFS")
curve(ax, d["bem"][:, 1], d["bem"][:, 2], ORANGE, "D", "-.", "Bempp, with symmetry")
ax.set_title("(b) torus", loc="left", fontsize=8.5)

ax = axs[1, 0]
d = np.load(data("cost_almond.npz"))
ref = d[f"mfs_E{len(d['mfs']) - 1}"]
curve(ax, d["mfs"][:-1, 2], [rel(d[f"mfs_E{i}"], ref) for i in range(len(d["mfs"]) - 1)], BLUE, "o", "-", "MFS")
curve(ax, d["bem"][:, 1], [rel(d[f"bem_E{i}"], ref) for i in range(len(d["bem"]))], ORANGE, "D", "-.", "Bempp, with symmetry")
ax.set_title("(c) almond, 1.19 GHz, nose-on", loc="left", fontsize=8.5)
ax.set_xticks([2, 5, 10, 20])
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, pos: f"{v:g}"))

ax = axs[1, 1]
d = np.load(data("cube_convergence.npz"))
ref, jref = d["ref_E"], int(d["ref_info"][0])
for key, col, mk, ls, lab in (("dense", BLUE, "o", "-", "MFS"), ("tail", NAVY, "h", (0, (6, 1.5)), "MFS, one depth per octave")):
    rows = [r for r in d[key] if r[0] >= 4 and not (key == "dense" and int(r[0]) == jref)]
    curve(ax, [r[1] for r in rows], [rel(d[f"{key}_E{int(r[0])}"], ref) for r in rows], col, mk, ls, lab)
u = np.load(data("cube_bempp_uniform.npz"))
keep = u["ndof"] > 1000
curve(ax, u["t"][keep], [rel(E, ref) for E in u["E"][keep]], ORANGE, "D", "-.", "Bempp, quasi-uniform")
env = np.load(data("cube_cases_k1.npz"))["env"]
curve(ax, env[:, 5], env[:, 6], RUST, "X", (0, (3, 1)), "Bempp, graded, with symmetry")
ax.set_title("(d) cube, side one wavelength", loc="left", fontsize=8.5)
ax.set_ylim(8e-7, 3e-2)

for ax in axs.ravel():
    ax.legend(frameon=False, fontsize=6.8, loc="upper right" if ax is axs[1, 0] else "lower left", handlelength=2.6)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for ax in axs[1]:
    ax.set_xlabel("wall time [s]  (4 cores)")
for ax in axs[:, 0]:
    ax.set_ylabel("relative error of the scattered field")
finish(fig, axs, "fig_cost_all")
