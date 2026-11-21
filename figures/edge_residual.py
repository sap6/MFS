"""
Figure 5: the unit cube at J = 6, tangential residual on a face as a function of
the distance to an edge, with and without the singular tails, showing the
discrete free and tail source depths along the top axis.

    python figures/edge_residual.py          # the figure, from data/edge_residual.npz
    python figures/edge_residual.py --run    # solve again (about a minute)
"""
import sys

import numpy as np
from style import BLUE, GREEN, MUTED, ORANGE, data, finish, plt

from mfs import core, cube

J = 6

if "--run" in sys.argv:
    # Root mean square over five positions along the edge, on the face y = 1/2, near the edge x = y = 1/2
    d, p = core.wave()
    rho = np.logspace(-6, np.log10(0.3), 120)
    x = np.array([(0.5 - r, 0.5, z) for z in (0.03, 0.08, 0.13, 0.18, 0.23) for r in rho])
    inc = core.plane_wave(x, 2 * np.pi, d, p)
    out = {"rho": rho, "J": J}
    for name, tails in (("tail", True), ("notail", False)):
        c = cube.Cube(J, tails=tails, classes=core.classes_of(2 * np.pi, d, p, cube.REFL))
        tot = c.field(c.solve(d, p)[0], x) + inc
        out[name] = np.sqrt(np.mean(np.abs(tot[:, [0, 2]].reshape(5, len(rho), 2)) ** 2, axis=(0, 2)))
    out["dJ"] = c.depths[-1]
    np.savez(data("edge_residual.npz"), **out)

# The figure
res = np.load(data("edge_residual.npz"))
fig, ax = plt.subplots(figsize=(4.6, 3.15))
ax.loglog(res["rho"], res["notail"], "--", color=GREEN, lw=1.6, label="edge sources, no tail")
ax.loglog(res["rho"], res["tail"], "-", color=BLUE, lw=1.6, label="edge sources and tail")

dJ = float(res["dJ"])
d_free = 0.2 * (2.0 ** -0.5) ** np.arange(0, 2 * int(res["J"]) + 1)
d_tail = dJ * (0.5 ** np.arange(1, 10))
d_end = d_tail[-1]

# Source-depth rug at the top of the plot showing exact free (2/octave) and tail (1/octave) depths
y_rug = 38.0
ax.plot(d_free, np.full_like(d_free, y_rug), "o", color=BLUE, ms=3.4, mec="white", mew=0.4, zorder=4)
ax.plot(d_tail, np.full_like(d_tail, y_rug), "o", mfc="white", mec=ORANGE, ms=3.2, mew=0.8, zorder=4)
ax.text(np.sqrt(d_free[0] * dJ), y_rug * 1.55, "free depths (2/oct)", color=BLUE, fontsize=7.2, ha="center")
ax.text(np.sqrt(d_tail[0] * d_end), y_rug * 1.55, "9-line tail (1/oct)", color=ORANGE, fontsize=7.2, ha="center")

ax.axvline(dJ, color=MUTED, lw=0.75, ls=":", zorder=1)
ax.axvline(d_end, color=MUTED, lw=0.75, ls=":", zorder=1)
ax.text(dJ * 1.15, 1.8e-5, r"$d_{\min}$", color=MUTED, fontsize=8)
ax.text(d_end * 1.15, 2.5e-3, r"$2^{-9}d_{\min}$", color=MUTED, fontsize=8)
ax.set(ylim=(3e-6, 110.0), xlabel=r"distance $\rho$ to the edge",
       ylabel=r"tangential residual  $|\mathbf{t}\cdot(\mathbf{E}_N+\mathbf{E}^{\mathrm{inc}})|$")
ax.legend(frameon=False, fontsize=8.0, loc="lower left", bbox_to_anchor=(0.12, 0.01))
ax.set_title(rf"unit cube, $k = 2\pi$, $J = {int(res['J'])}$: residual near an edge", loc="left", fontsize=8.5, pad=10)
finish(fig, ax, "fig_edge_residual")
