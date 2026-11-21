"""
Figure 1: the unit sphere at k = 2π, error against the Mie series.
(a) Equatorial cutaway showing outer boundary collocation points M = 2N and
inner Fibonacci source dipoles N at radius R_s; (b) error against √N; (c)
dependence on R_s for N = 400.

    python figures/sphere.py          # the figure, from data/sphere_convergence.json
    python figures/sphere.py --run    # solve again (about 20 seconds)
"""
import json
import os
import sys

import numpy as np
from matplotlib.patches import FancyArrowPatch
from style import BLUE, BODY, GREEN, INK, MUTED, ORANGE, PURPLE, RUST, curve, data, finish, plt, rel

from mfs import core, mie, sphere

K = 2 * np.pi

if "--run" in sys.argv:
    d, p = core.wave()
    u = core.fibonacci(1500)
    exact, T = mie.scattered_field(1.5 * u, K, d, p), core.tangent_basis(u)
    inc = np.einsum("mab,mb->ma", T, core.plane_wave(u, K, d, p))
    res = []
    Ns = (50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200)
    runs = [("N", R, N) for R in (0.3, 0.5, 0.7, 0.9) for N in Ns] + [("R", R, 400) for R in np.round(np.linspace(0.1, 0.95, 18), 3)]
    for sweep, R, N in runs:
        y, alpha = sphere.solve(N, K, R)
        res.append(dict(k=K, sweep=sweep, R=float(R), N=N, mid=rel(core.field(1.5 * u, y, alpha, K), exact),
                        bc=rel(np.einsum("mab,mb->ma", T, core.field(u, y, alpha, K)), -inc),
                        alpha_norm=float(np.linalg.norm(alpha))))
        print(sweep, R, N, f"{res[-1]['mid']:.2e}", flush=True)
    path = data("sphere_convergence.json")
    old = [r for r in json.load(open(path)) if abs(r["k"] - K) > 1e-9] if os.path.exists(path) else []
    json.dump(res + old, open(path, "w"), indent=1)

res = [r for r in json.load(open(data("sphere_convergence.json"))) if abs(r["k"] - K) < 1e-9]
fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.75), gridspec_kw=dict(width_ratios=[0.94, 1.16, 1.16]))

# (a) Equatorial cutaway showing boundary collocation points (M = 2N) and inner dipoles (N at R_s)
ax = axs[0]
th = np.linspace(0, 2 * np.pi, 240)
ax.fill(np.cos(th), np.sin(th), color=BODY, lw=0, zorder=0)
ax.plot(np.cos(th), np.sin(th), color=INK, lw=1.2, zorder=2)

# Illustrative ring slice: 24 boundary points on r = 1, 12 source dipoles on r = R_s = 0.5
Rs = 0.5
ang_b = 2 * np.pi * np.arange(24) / 24
ang_s = 2 * np.pi * (np.arange(12) + 0.5) / 12
ax.plot(np.cos(ang_b), np.sin(ang_b), ".", color=INK, ms=3.4, zorder=4, label=r"boundary $\Gamma$ ($M=2N$)")
ax.plot(Rs * np.cos(th), Rs * np.sin(th), "-.", color=ORANGE, lw=0.85, alpha=0.75, zorder=2)
ax.plot(Rs * np.cos(ang_s), Rs * np.sin(ang_s), "o", color=ORANGE, ms=3.3, mec="white", mew=0.4,
        zorder=4, label=r"dipoles at $R_s$ ($N$)")
ax.plot(0.0, 0.0, "*", color=RUST, ms=6.5, zorder=5, label="singularity ($r=0$)")

# Show one tangent vector at a boundary point and R_s radius arrow
ang0 = np.deg2rad(45.0)
bx0, by0 = np.cos(ang0), np.sin(ang0)
ax.add_patch(FancyArrowPatch((bx0, by0), (bx0 - 0.34 * by0, by0 + 0.34 * bx0),
                             arrowstyle="-|>", mutation_scale=7, color=INK, lw=0.85, zorder=5))
ax.text(bx0 - 0.36 * by0 + 0.04, by0 + 0.36 * bx0, r"$\mathbf{t}_i$", fontsize=7.5, color=INK)
ax.annotate("", xy=(Rs * np.cos(-0.35), Rs * np.sin(-0.35)), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.8))
ax.text(0.24, -0.22, "$R_s$", fontsize=7.5, color=ORANGE)

ax.set(xlim=(-1.22, 1.22), ylim=(-1.18, 1.22), aspect="equal")
ax.axis("off")
ax.legend(frameon=False, fontsize=6.9, loc="upper center", bbox_to_anchor=(0.50, -0.06), ncol=1,
          handletextpad=0.3, labelspacing=0.25)
ax.set_title("(a) placement ($R_s < 1$)", loc="left", fontsize=8.4)

# (b) Convergence vs sqrt(N)
ax1 = axs[1]
for (col, mk, ls), R in zip(((BLUE, "o", "-"), (ORANGE, "D", "-."), (GREEN, "s", "--"), (PURPLE, "^", ":")), (0.3, 0.5, 0.7, 0.9)):
    s = sorted((r for r in res if r["sweep"] == "N" and r["R"] == R), key=lambda r: r["N"])
    curve(ax1, [np.sqrt(r["N"]) for r in s], [r["mid"] for r in s], col, mk, ls, f"$R_s = {R}$", 4.1, 1.35, 0.55, logx=False)
ax1.set(xlabel=r"$\sqrt{N}$ ($3N$ unknowns)", ylabel="relative error against Mie")
ax1.set_title(r"(b) error against $\sqrt{N}$", loc="left", fontsize=8.4)
ax1.legend(frameon=False, fontsize=7.5, loc="lower left")

# (c) Dependence on R_s at N = 400
ax2 = axs[2]
s = sorted((r for r in res if r["sweep"] == "R"), key=lambda r: r["R"])
R = [r["R"] for r in s]
curve(ax2, R, [r["mid"] for r in s], BLUE, "o", "-", "field error", 4.1, 1.35, 0.55, logx=False)
curve(ax2, R, [r["bc"] for r in s], ORANGE, "D", "-.", "boundary residual", 4.1, 1.35, 0.55, logx=False)
curve(ax2, R, [r["alpha_norm"] for r in s], MUTED, "s", ":", r"norm $\|\boldsymbol{\alpha}\|$", 3.8, 1.2, 0.5, logx=False)
ax2.set(xlabel="source radius $R_s$", ylabel=r"error, residual, $\|\boldsymbol{\alpha}\|$")
ax2.set_title(rf"(c) $N = {s[0]['N']}$: role of $R_s$", loc="left", fontsize=8.4)
ax2.legend(frameon=False, fontsize=7.3, loc="upper right")

for ax in (ax1, ax2):
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.02, right=0.99, top=0.90, bottom=0.25, wspace=0.33)
finish(fig, None, "fig_sphere", tight="tight")
