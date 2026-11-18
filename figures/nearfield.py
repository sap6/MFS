"""
Figure 11: the unit sphere at k = 2π, error of the scattered field against the
Mie series as a function of the distance from the surface, for the MFS and Bempp.

    python figures/nearfield.py          # the figure, from data/sphere_nearfield.npz
    python figures/nearfield.py --run    # solve again (needs bempp-cl; a few minutes)
"""
import sys

import numpy as np
from style import BLUE, NAVY, ORANGE, RUST, data, finish, plt

from mfs import core, mie, sphere

K = 2 * np.pi

if "--run" in sys.argv:
    # 200 points on each of 17 spheres at distances 1e-4 to 1 from the surface
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    bempp = bem.bempp
    d, p = core.wave()
    delta = np.logspace(-4, 0, 17)
    pts = np.concatenate([(1 + dl) * core.fibonacci(200) for dl in delta])
    exact = mie.scattered_field(pts, K, d, p)
    shell = lambda E: np.array([core.rel(E[i * 200:(i + 1) * 200], exact[i * 200:(i + 1) * 200]) for i in range(17)])
    out = {"delta": delta}
    for N in (400, 1200):
        y, alpha = sphere.solve(N)
        E = core.field(pts, y, alpha, K)
        out[f"mfs{N}"] = shell(E)
    for h in (0.1, 0.07):
        E, N = bem.efie(bempp.shapes.sphere(h=h), K, d, p, pts)
        out[f"bem{N}"] = shell(E)
    np.savez(data("sphere_nearfield.npz"), **out)

# The figure
res = np.load(data("sphere_nearfield.npz"))
fig, ax = plt.subplots(figsize=(5.0, 3.3))
styles = {"mfs": [(BLUE, "o", "-"), (NAVY, "h", (0, (6, 1.5)))], "bem": [(ORANGE, "D", "-."), (RUST, "X", (0, (3, 1)))]}
for key in [k for k in res if k != "delta"]:
    col, mk, ls = styles[key[:3]].pop(0)
    label = f"MFS, $N = {int(key[3:])}$" if key[:3] == "mfs" else f"Bempp, {int(key[3:])} unknowns"
    ax.loglog(res["delta"], res[key], linestyle=ls, color=col, marker=mk, ms=4.5, lw=1.4, mec="white", mew=0.6, label=label)
ax.set_xlabel("distance from the surface")
ax.set_ylabel("relative error against Mie")
ax.legend(frameon=False, fontsize=8, loc="center", bbox_to_anchor=(0.5, 0.24), ncol=2)
finish(fig, ax, "fig_sphere_nearfield")
