"""
Figure 5: the unit cube at J = 6, tangential residual on a face as a function of
the distance to an edge, with and without the singular tails.

    python figures/edge_residual.py          # the figure, from data/edge_residual.npz
    python figures/edge_residual.py --run    # solve again (about a minute)
"""
import sys

import numpy as np
from style import BLUE, GREEN, MUTED, data, finish, plt

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
fig, ax = plt.subplots(figsize=(4.5, 3.0))
ax.loglog(res["rho"], res["notail"], "--", color=GREEN, lw=1.6, label="edge sources, no tail")
ax.loglog(res["rho"], res["tail"], "-", color=BLUE, lw=1.6, label="edge sources and tail")
ax.axvline(res["dJ"], color=MUTED, lw=0.7, ls=":")
ax.text(res["dJ"] * 1.15, 1.6e-5, r"$d_{\min}$", color=MUTED, fontsize=8.5)
ax.set(xlabel=r"distance $\rho$ to the edge", ylabel=r"tangential residual  $|\mathbf{t}\cdot(\mathbf{E}_N+\mathbf{E}^{\mathrm{inc}})|$")
ax.legend(frameon=False, fontsize=8.5, loc="lower left")
ax.set_title(f"Unit cube, $J = {int(res['J'])}$: residual on a face near an edge", loc="left", fontsize=9)
finish(fig, ax, "fig_edge_residual")
