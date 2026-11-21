"""
Figure 2: the torus R = 1 with three minor radii a, k = 2π, 19200 unknowns.
(a) Meridian cross-section showing the boundary collocation ring, the inner
source ring a_s, and the toroidal focal ring sqrt(R² - a²); (b) boundary
residual and (c) coefficient norm against a_s / a.

    python figures/torus_sources.py          # the figure, from data/torus_sources.npz
    python figures/torus_sources.py --run    # solve again (about 15 seconds)
"""
import sys

import numpy as np
from style import BLUE, BODY, GREEN, INK, MUTED, ORANGE, RUST, curve, data, finish, plt

from mfs import core, torus

K = 2 * np.pi
TUBES = {0.25: (BLUE, "o", "-"), 0.5: (ORANGE, "D", "-."), 0.8: (GREEN, "s", "--")}
FRACS = (0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05)

if "--run" in sys.argv:
    d, p = core.wave()
    out = {}
    for a in TUBES:
        rows = []
        for f in FRACS:
            s = torus.scatter(K, d, p, 40, 160, a=a, a_s=f * a, rings=3)
            rows.append((f, torus.residual(s, d, p, a=a), np.linalg.norm(s.alpha), torus.optical(s, d, p)[1]))
            print(f"a = {a}, a_s/a = {f}: residual {rows[-1][1]:.1e}, |alpha| {rows[-1][2]:.1e}", flush=True)
        out[f"a{a:g}"] = np.array(rows)
    np.savez(data("torus_sources.npz"), **out)

# Three-panel figure: (a) meridian geometry & focal ring, (b) residual, (c) coefficient norm
res = np.load(data("torus_sources.npz"))
fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.75), gridspec_kw=dict(width_ratios=[0.96, 1.15, 1.15]))

# Panel (a): meridian cross-section (rho, z) for R = 1, a = 0.8 (focal ring at rho = 0.6, bound a_s = 0.4)
ax = axs[0]
R, a, a_good, a_bad = 1.0, 0.8, 0.52, 0.26
th = np.linspace(0, 2 * np.pi, 240)
ax.fill(R + a * np.cos(th), a * np.sin(th), color=BODY, lw=0, zorder=0)
ax.plot(R + a * np.cos(th), a * np.sin(th), color=INK, lw=1.2, zorder=2)

# Boundary collocation points (2 n_s = 28 illustrative) and inner source dipoles (n_s = 14)
rb, zb, _ = torus.circle(28, R, a)
rs, zs, _ = torus.circle(14, R, a_good)
ax.plot(rb, zb, ".", color=INK, ms=3.2, zorder=4, label=r"boundary $\Gamma$ ($2n_s$)")
ax.plot(R + a_good * np.cos(th), a_good * np.sin(th), "-", color=GREEN, lw=0.85, alpha=0.7, zorder=2)
ax.plot(rs, zs, "o", color=GREEN, ms=3.2, mec="white", mew=0.35, zorder=4, label=r"valid $a_s$ (encloses)")
ax.plot(R + a_bad * np.cos(th), a_bad * np.sin(th), ":", color=RUST, lw=1.0, zorder=2, label=r"too thin ($a_s < \text{bound}$)")

# Symmetry axis rho = 0, tube centre rho = R, and focal ring rho = sqrt(R^2 - a^2) = 0.6
rho_f = np.sqrt(R * R - a * a)
ax.axvline(0.0, color=MUTED, lw=0.8, ls=(0, (4, 2)), zorder=1)
ax.plot(R, 0.0, "+", color=MUTED, ms=5.0, mew=0.9, zorder=3)
ax.plot(rho_f, 0.0, "*", color=RUST, ms=6.8, zorder=6)
ax.text(0.03, 0.78, "axis\n$\\rho=0$", fontsize=7.0, color=MUTED, va="top")
ax.annotate("focal ring\n$\\sqrt{R^2-a^2}$", xy=(rho_f, -0.05), xytext=(0.64, -0.38),
            fontsize=6.6, color=RUST, ha="center", va="top",
            bbox=dict(boxstyle="square,pad=0.08", fc="#eeede8", ec="none"),
            arrowprops=dict(arrowstyle="-", color=RUST, lw=0.55, shrinkA=0, shrinkB=2))
ax.text(R + 0.03, 0.03, "$R$", fontsize=7.3, color=MUTED)
ax.set(xlim=(-0.06, 1.88), ylim=(-0.92, 0.92), aspect="equal", xlabel=r"radial $\rho$", ylabel="$z$",
       xticks=[0, 0.6, 1.0, 1.8])
ax.legend(frameon=False, fontsize=6.9, loc="upper center", bbox_to_anchor=(0.48, -0.25), ncol=1,
          handletextpad=0.3, labelspacing=0.25)
ax.set_title("(a) meridian cut ($a/R = 0.8$)", loc="left", fontsize=8.4)
ax.grid(False)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

for a_val, (col, mk, ls) in TUBES.items():
    r = res[f"a{a_val:g}"]
    for ax, j in zip(axs[1:], (1, 2)):
        curve(ax, r[:, 0], r[:, j], col, mk, ls, rf"$a/R = {a_val:g}$", 4.2, 1.35, 0.55, logx=False)
        ax.axvline((1 - np.sqrt(1 - a_val * a_val)) / a_val, color=col, lw=0.9, ls=":")

axs[1].set(xlabel=r"source torus radius $a_s/a$", ylabel="boundary residual")
axs[2].set(xlabel=r"source torus radius $a_s/a$", ylabel=r"coefficient norm $\|\boldsymbol{\alpha}\|$")
axs[2].legend(frameon=False, fontsize=7.6, loc="upper right")
axs[1].set_title(r"(b) residual ($19200$ unknowns)", loc="left", fontsize=8.4)
axs[2].set_title(r"(c) $\|\boldsymbol{\alpha}\|$ (dotted: focal bound)", loc="left", fontsize=8.4)
for ax in axs[1:]:
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.06, right=0.99, top=0.90, bottom=0.28, wspace=0.33)
finish(fig, None, "fig_torus_sources", tight="tight")
