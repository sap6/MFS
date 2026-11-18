"""
Figure 6 and Table 2: the unit cube at k = 2π, error of the scattered field
against √N_u, N_u the number of unknowns of the whole cube. The errors are
differences from the finest layout, J = 12.

    python figures/cube_convergence.py          # the figure, from data/
    python figures/cube_convergence.py --run    # solve again (about 15 minutes; Bempp: half an hour more)
"""
import sys

import numpy as np
from style import BLUE, GREEN, INK, MUTED, NAVY, ORANGE, PURPLE, RUST, data, plt, rel

from mfs import core, cube

if "--run" in sys.argv:
    # The paper's layout (two depths per octave), one depth per octave, no tails, and the usual layout (an inner cube)
    pts = 2.0 * core.fibonacci(400)
    out = {"pts": pts}
    for name, Js, kw in (("dense", range(3, 13), {}), ("tail", range(3, 10), dict(q=0.5)),
                         ("notail", (4, 6, 8), dict(tails=False)), ("plain", (4, 6, 8, 10), None)):
        rows = []
        for J in Js:
            E, t, info = cube.solve(None, pts=pts, n_base=J) if kw is None else cube.solve(J, pts=pts, **kw)
            rows.append((J, t, sum(info["unknowns"]), info["rows"], info["optical"]))
            out[f"{name}_E{J}"] = E
            print(name, J, f"{t:.1f} s", flush=True)
        out[name] = np.array(rows)
    out["ref_E"], out["ref_info"] = out["dense_E12"], np.array([12, *out["dense"][-1][1:], 0.0])
    np.savez(data("cube_convergence.npz"), **out)

    # Bempp as distributed: quasi-uniform meshes, no symmetry, quadrature orders 4 and 6
    import bempp_cl.api as bempp
    from mfs import bem
    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = 4, 6
    rows = []
    for h in (0.3, 0.1, 0.07, 0.05, 0.04):
        t0 = __import__("time").time()
        E, N = bem.efie(bempp.shapes.cube(length=1, origin=(-0.5, -0.5, -0.5), h=h), 2 * np.pi, *core.wave(), pts)
        rows.append((N, __import__("time").time() - t0, E, h))
    np.savez(data("cube_bempp_uniform.npz"), ndof=[r[0] for r in rows], t=[r[1] for r in rows], E=[r[2] for r in rows], h=[r[3] for r in rows])

# Errors against the reference, as a function of the unknowns of the whole cube
d = np.load(data("cube_convergence.npz"))
ref, jref = d["ref_E"], int(d["ref_info"][0])
curves = {name: np.array([(4 * r[2], rel(d[f"{name}_E{int(r[0])}"], ref)) for r in d[name]
                          if (name not in ("tail", "dense") or r[0] >= 4) and not (name == "dense" and int(r[0]) == jref)])
          for name in ("dense", "tail", "notail", "plain")}
u = np.load(data("cube_bempp_uniform.npz"))
curves["uniform"] = np.array([(n, rel(E, ref)) for n, E in zip(u["ndof"], u["E"]) if n > 1000])
best = {}                                                        # Bempp on graded meshes: the smallest error on each mesh
for g, reg, sing, n, N, t, err in np.load(data("cube_cases_k1.npz"))["bem"]:
    best[int(N)] = min(best.get(int(N), 1.0), err)
curves["graded"] = np.array(sorted(best.items()))

# Root-exponential fit to the levels 4 to 10 of the paper's layout (J = 11 is too close to the reference)
c = curves["dense"][:-1]
p = np.polyfit(np.sqrt(c[:, 0]), np.log(c[:, 1]), 1)
print(f"error = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N_u)), largest deviation a factor "
      f"{np.exp(np.abs(np.log(c[:, 1]) - np.polyval(p, np.sqrt(c[:, 0]))).max()):.2f}")

# The figure
fig, ax = plt.subplots(figsize=(6.0, 4.4))
xs = np.array([np.sqrt(curves["dense"][0, 0]) - 4, np.sqrt(curves["dense"][-1, 0]) + 8])
ax.semilogy(xs, np.exp(np.polyval(p, xs)), color=MUTED, lw=0.8, zorder=1)
ax.annotate(rf"${np.exp(p[1]):.2f}\,\exp({p[0]:.3f}\sqrt{{N_u}})$", (xs[1] - 42, np.exp(np.polyval(p, xs[1] - 42)) * 3.2), fontsize=8, color=MUTED)
for key, (col, mk, ls, label) in {"plain": (PURPLE, "^", ":", "plain MFS (inner cube)"),
                                  "notail": (GREEN, "s", "--", "MFS without the tails"),
                                  "tail": (NAVY, "h", (0, (6, 1.5)), "MFS, one depth per octave"),
                                  "dense": (BLUE, "o", "-", "MFS (two depths per octave)"),
                                  "uniform": (ORANGE, "D", "-.", "Bempp, quasi-uniform meshes"),
                                  "graded": (RUST, "X", (0, (3, 1)), "Bempp, graded meshes (best setting per mesh)")}.items():
    v = curves[key]
    ax.semilogy(np.sqrt(v[:, 0]), v[:, 1], linestyle=ls, color=col, marker=mk, ms=5, lw=1.4, mec="white", mew=0.7, label=label, zorder=3)
ax.set_xlabel(r"$\sqrt{N_u}$,  $N_u$ = number of unknowns")
ax.set_ylabel("relative error of the scattered field")
ax.grid(True, which="major", color="#e4e3df", lw=0.6)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, columnspacing=1.2)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{__import__('style').HERE}/fig_cube_convergence.{ext}", dpi=200, bbox_inches="tight")
