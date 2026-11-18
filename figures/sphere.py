"""
Figure 1: the unit sphere at k = 2π, error against the Mie series.
Left: against √N for four source radii. Right: against the source radius for N = 400.

    python figures/sphere.py          # the figure, from data/sphere_convergence.json
    python figures/sphere.py --run    # solve again (about two minutes)
"""
import json
import sys

import numpy as np
from style import BLUE, GREEN, MUTED, ORANGE, PURPLE, data, finish, plt

from mfs import core, mie, sphere

K = 2 * np.pi

if "--run" in sys.argv:
    # Errors on the sphere of radius 1.5, the boundary residual, and the coefficient norm
    d, p = core.wave()
    u = core.fibonacci(1500)
    exact, T = mie.scattered_field(1.5 * u, K, d, p), core.tangent_basis(u)
    res = []
    Ns = (50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200)
    runs = [("N", R, N) for R in (0.3, 0.5, 0.7, 0.9) for N in Ns] + [("R", R, 400) for R in np.round(np.linspace(0.1, 0.95, 18), 3)]
    for sweep, R, N in runs:
        y, alpha = sphere.solve(N, K, R)
        bc = np.einsum("mab,mb->ma", T, core.field(u, y, alpha, K) + core.plane_wave(u, K, d, p))
        res.append(dict(k=K, sweep=sweep, R=float(R), N=N, mid=core.rel(core.field(1.5 * u, y, alpha, K), exact),
                        bc=float(np.linalg.norm(bc) / np.linalg.norm(np.einsum("mab,mb->ma", T, core.plane_wave(u, K, d, p)))),
                        alpha_norm=float(np.linalg.norm(alpha))))
        print(sweep, R, N, f"{res[-1]['mid']:.2e}", flush=True)
    json.dump(res, open(data("sphere_convergence.json"), "w"), indent=1)

# The figure, from the k = 2π runs
res = [r for r in json.load(open(data("sphere_convergence.json"))) if abs(r["k"] - K) < 1e-9]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 2.9))
for (col, mk, ls), R in zip(((BLUE, "o", "-"), (ORANGE, "D", "-."), (GREEN, "s", "--"), (PURPLE, "^", ":")), (0.3, 0.5, 0.7, 0.9)):
    s = sorted((r for r in res if r["sweep"] == "N" and r["R"] == R), key=lambda r: r["N"])
    ax1.semilogy([np.sqrt(r["N"]) for r in s], [r["mid"] for r in s], ls, color=col, marker=mk, ms=4.5, lw=1.4,
                 mec="white", mew=0.6, label=f"$R_s = {R}$")
ax1.set(xlabel=r"$\sqrt{N}$", ylabel="relative error against Mie")
ax1.set_title(r"$k = 2\pi$: error against the number of sources", loc="left", fontsize=8.5)
ax1.legend(frameon=False, fontsize=7.5, loc="lower left")

s = sorted((r for r in res if r["sweep"] == "R"), key=lambda r: r["R"])
R = [r["R"] for r in s]
ax2.semilogy(R, [r["mid"] for r in s], "-", color=BLUE, marker="o", ms=4.5, lw=1.4, mec="white", mew=0.6, label="error of the field")
ax2.semilogy(R, [r["bc"] for r in s], "-.", color=ORANGE, marker="D", ms=4.5, lw=1.4, mec="white", mew=0.6, label="boundary residual")
ax2.semilogy(R, [r["alpha_norm"] for r in s], ":", color=MUTED, marker="s", ms=4, lw=1.2, mec="white", mew=0.5,
             label="coefficient norm")
ax2.set(xlabel="source radius $R_s$")
ax2.set_title(rf"$k = 2\pi$, $N = {s[0]['N']}$: dependence on the source radius", loc="left", fontsize=8.5)
ax2.legend(frameon=False, fontsize=7.5, loc="upper center")
finish(fig, (ax1, ax2), "fig_sphere")
