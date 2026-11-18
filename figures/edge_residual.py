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
    T = np.tile([[1.0, 0, 0], [0, 0, 1.0]], (len(x), 1, 1))
    out = {"rho": rho, "J": J}
    for name, tails in (("tail", True), ("notail", False)):
        c = cube.Cube(J, tails=tails, classes=cube.classes_of(2 * np.pi, d, p))
        coef, _, _ = c.solve(d, p)
        A = c.matrix(x, T, list(coef))
        E = sum(A[i][:, :c.unknowns[cl]] @ coef[cl] for i, cl in enumerate(coef)).reshape(-1, 2)
        r = np.abs(E + np.einsum("mab,mb->ma", T, core.plane_wave(x, 2 * np.pi, d, p))) ** 2
        out[name] = np.sqrt(np.mean(r.reshape(5, len(rho), 2), axis=(0, 2)))
        out["dJ"] = c.depths[-1]
    np.savez(data("edge_residual.npz"), **out)

# The figure
d = np.load(data("edge_residual.npz"))
fig, ax = plt.subplots(figsize=(5.2, 3.4))
ax.loglog(d["rho"], d["notail"], "--", color=GREEN, lw=1.6, label="edge sources, no tail")
ax.loglog(d["rho"], d["tail"], "-", color=BLUE, lw=1.6, label="edge sources and tail")
ax.axvline(d["dJ"], color=MUTED, lw=0.7, ls=":")
ax.text(d["dJ"] * 1.15, ax.get_ylim()[0] * 2.5, r"$d_{\min}$", color=MUTED, fontsize=8)
ax.set(xlabel="distance to the edge", ylabel="tangential residual  $|\\mathbf{t}\\cdot(\\mathbf{E}_N+\\mathbf{E}^{\\rm inc})|$")
ax.legend(frameon=False, fontsize=8, loc="lower left")
ax.set_title(f"Unit cube, $J = {int(d['J'])}$: residual on a face near an edge", loc="left", fontsize=9)
finish(fig, ax, "fig_edge_residual")
