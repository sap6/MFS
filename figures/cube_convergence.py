"""
Figure 6 and Table 2: the unit cube at k = 2π, error of the scattered field
against √N_u, N_u the number of unknowns of the whole cube. The errors are
differences from the finest layout, J = 12.

    python figures/cube_convergence.py          # the figure, from data/
    python figures/cube_convergence.py --run    # solve again (about 12 minutes; Bempp: 8 minutes more)
"""
import gc
import sys
import time

import numpy as np
from style import BLUE, GREEN, MUTED, NAVY, ORANGE, PURPLE, RUST, curve, data, finish, plt, rel

from mfs import core, cube

if "--run" in sys.argv:
    # The paper's layout (two depths per octave), one depth per octave, no tails, and the usual layout (an inner cube)
    pts = 2.0 * core.fibonacci(400)
    cube.solve(2, pts=pts)                                              # warm-up: compilation is not timed
    out = {"pts": pts}
    for name, Js, kw in (("dense", range(3, 13), {}), ("tail", range(3, 10), dict(q=0.5)),
                         ("notail", (4, 6, 8), dict(tails=False)), ("plain", (4, 6, 8, 10), None)):
        reps = 2 if name in ("dense", "tail") else 1
        rows = []
        for J in Js:
            runs_J = []
            for _ in range(reps):
                gc.collect()
                runs_J.append(cube.solve(None, pts=pts, n_base=J) if kw is None else cube.solve(J, pts=pts, **kw))
            E, t, info = min(runs_J, key=lambda r: r[1])
            rows.append((J, t, sum(info["unknowns"]), info["rows"], info["optical"]))
            out[f"{name}_E{J}"] = E
            print(name, J, f"{t:.1f} s", flush=True)
        out[name] = np.array(rows)
    out["ref_E"], out["ref_info"] = out["dense_E12"], np.array([12, *out["dense"][-1][1:], 0.0])
    np.savez(data("cube_convergence.npz"), **out)

    # Bempp as distributed: quasi-uniform meshes, no symmetry, quadrature orders 4 and 6
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    bempp = bem.bempp
    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = 4, 6
    bem.efie(bempp.shapes.cube(length=1, origin=(-0.5, -0.5, -0.5), h=0.4), 2 * np.pi, *core.wave(), pts)
    runs = []
    for h in (0.3, 0.1, 0.07, 0.05, 0.04):
        grid = bempp.shapes.cube(length=1, origin=(-0.5, -0.5, -0.5), h=h)
        gc.collect()
        t0 = time.time()
        E, N = bem.efie(grid, 2 * np.pi, *core.wave(), pts)
        runs.append((N, time.time() - t0, E, h))
        print("bempp uniform", h, N, f"{runs[-1][1]:.1f} s", flush=True)
    ndof, t, E, h = zip(*runs)
    np.savez(data("cube_bempp_uniform.npz"), ndof=ndof, t=t, E=E, h=h)

# Errors against the reference, as a function of the unknowns of the whole cube
res = np.load(data("cube_convergence.npz"))
ref, jref = res["ref_E"], int(res["ref_info"][0])
curves = {name: np.array([(4 * r[2], rel(res[f"{name}_E{int(r[0])}"], ref)) for r in res[name] if 4 <= r[0] < jref])
          for name in ("dense", "tail", "notail", "plain")}
u = np.load(data("cube_bempp_uniform.npz"))
curves["uniform"] = np.array([(n, rel(E, ref)) for n, E in zip(u["ndof"], u["E"]) if n > 1000])
best = {}                                                        # Bempp on graded meshes: the smallest error on each mesh
for *_, N, _, err in np.load(data("cube_cases_k1.npz"))["bem"]:
    best[int(N)] = min(best.get(int(N), 1.0), err)
curves["graded"] = np.array(sorted(best.items()))

# Root-exponential fit to the levels 4 to 10 of the paper's layout (J = 11 is too close to the reference)
c = curves["dense"][:-1]
p = np.polyfit(np.sqrt(c[:, 0]), np.log(c[:, 1]), 1)
print(f"error = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N_u)), largest deviation a factor "
      f"{np.exp(np.abs(np.log(c[:, 1]) - np.polyval(p, np.sqrt(c[:, 0]))).max()):.2f}")

# The figure
fig, ax = plt.subplots(figsize=(5.2, 3.85))
xs = np.array([np.sqrt(curves["dense"][0, 0]) - 4, np.sqrt(curves["dense"][-1, 0]) + 8])
ax.semilogy(xs, np.exp(np.polyval(p, xs)), color=MUTED, lw=0.8, zorder=1)
ax.annotate(rf"${np.exp(p[1]):.2f}\,\exp({p[0]:.4f}\sqrt{{N_u}})$", (153, 5.5e-4), fontsize=8.5, color=MUTED)
for key, (col, mk, ls, label) in {"plain": (PURPLE, "^", ":", "plain MFS (inner cube)"),
                                  "notail": (GREEN, "s", "--", "MFS without the tails"),
                                  "tail": (NAVY, "h", (0, (6, 1.5)), "MFS (one depth per octave)"),
                                  "dense": (BLUE, "o", "-", "MFS (two depths per octave)"),
                                  "uniform": (ORANGE, "D", "-.", "Bempp, quasi-uniform meshes"),
                                  "graded": (RUST, "X", (0, (3, 1)), "Bempp, graded meshes (best setting per mesh)")}.items():
    v = curves[key]
    curve(ax, np.sqrt(v[:, 0]), v[:, 1], col, mk, ls, label, lw=1.4, logx=False, zorder=3)

# Level callouts linking markers to Table 2
d0_xy = (np.sqrt(curves["dense"][0, 0]), curves["dense"][0, 1])
d1_xy = (np.sqrt(curves["dense"][-1, 0]), curves["dense"][-1, 1])
t1_xy = (np.sqrt(curves["tail"][-1, 0]), curves["tail"][-1, 1])
ax.text(d0_xy[0] + 2.5, d0_xy[1] * 1.18, "$J=4$", fontsize=7.5, color=BLUE)
ax.text(d1_xy[0] + 3.5, d1_xy[1] * 1.25, "$J=11$", fontsize=7.5, color=BLUE)
ax.text(t1_xy[0] - 2.0, t1_xy[1] * 1.45, "$J=9$", fontsize=7.4, color=NAVY, ha="right")

ax.set(xlabel=r"$\sqrt{N_u}$,  $N_u$ = unknowns of the whole cube", ylabel="relative error of the scattered field")
ax.set_title(r"unit cube, $k = 2\pi$: error against $\sqrt{N_u}$", loc="left", fontsize=8.5)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, columnspacing=1.0)
finish(fig, ax, "fig_cube_convergence", tight="tight")
