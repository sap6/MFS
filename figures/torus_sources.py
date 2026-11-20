"""
Figure 2: the torus R = 1 with three minor radii a, k = 2π, 19200 unknowns.
Boundary residual and coefficient norm against the radius of the source torus;
the dotted lines mark the conjectured bound a_s = R - sqrt(R² - a²).

    python figures/torus_sources.py          # the figure, from data/torus_sources.npz
    python figures/torus_sources.py --run    # solve again (about a minute)
"""
import sys

import numpy as np
from style import BLUE, GREEN, ORANGE, curve, data, finish, plt

from mfs import core, torus

K = 2 * np.pi
TUBES = {0.25: (BLUE, "o", "-"), 0.5: (ORANGE, "D", "-."), 0.8: (GREEN, "s", "--")}
FRACS = (0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05)

if "--run" in sys.argv:
    # 40 rings of 160 dipoles, 120 rings of boundary points; columns a_s/a, residual, coefficient norm, optical mismatch
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

# The figure
res = np.load(data("torus_sources.npz"))
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.7))
for a, (col, mk, ls) in TUBES.items():
    r = res[f"a{a:g}"]
    for ax, j in zip(axes, (1, 2)):
        curve(ax, r[:, 0], r[:, j], col, mk, ls, rf"$a/R = {a:g}$", 4.5, 1.4, 0.6, logx=False)
        ax.axvline((1 - np.sqrt(1 - a * a)) / a, color=col, lw=0.9, ls=":")
axes[0].set_ylabel("boundary residual")
axes[1].set_ylabel(r"coefficient norm $\|\boldsymbol{\alpha}\|$")
for ax in axes:
    ax.set_xlabel(r"source torus radius $a_s/a$")
axes[1].legend(frameon=False, fontsize=8, loc="upper right")
axes[0].set_title(r"(a) boundary residual ($k = 2\pi$, $19200$ unknowns)", loc="left", fontsize=8.5)
axes[1].set_title(r"(b) coefficient norm (dotted: focal ring)", loc="left", fontsize=8.5)
finish(fig, axes, "fig_torus_sources")
