"""
Figure 11: the unit sphere at k = 2π, error of the scattered field against the
Mie series as a function of the distance from the surface, for the MFS and Bempp.

    python figures/nearfield.py          # the figure, from data/sphere_nearfield.npz
    python figures/nearfield.py --run    # solve again (needs bempp-cl; a few minutes)
"""
import sys

import numpy as np
from style import BLUE, NAVY, ORANGE, RUST, curve, data, finish, plt

from mfs import core, mie, sphere

K = 2 * np.pi

if "--run" in sys.argv:
    # 200 points on each of 17 spheres at distances 1e-4 to 1 from the surface
    from mfs import bem                         # before bempp: it lets the OpenCL driver see all of the memory
    d, p = core.wave()
    delta = np.logspace(-4, 0, 17)
    pts = np.concatenate([(1 + dl) * core.fibonacci(200) for dl in delta])
    exact = mie.scattered_field(pts, K, d, p)
    shell = lambda E: np.array([core.rel(E[i * 200:(i + 1) * 200], exact[i * 200:(i + 1) * 200]) for i in range(17)])
    out = {"delta": delta}
    for N in (400, 1200):
        out[f"mfs{N}"] = shell(core.field(pts, *sphere.solve(N), K))
    for h in (0.1, 0.07):
        E, N = bem.efie(bem.bempp.shapes.sphere(h=h), K, d, p, pts)
        out[f"bem{N}"] = shell(E)
    np.savez(data("sphere_nearfield.npz"), **out)

# The figure
res = np.load(data("sphere_nearfield.npz"))
fig, ax = plt.subplots(figsize=(4.4, 2.9))
for key, col, mk, ls, label in (
    ("mfs400", BLUE, "o", "-", "MFS, 1200 unknowns"),
    ("mfs1200", NAVY, "h", (0, (6, 1.5)), "MFS, 3600 unknowns"),
    ("bem4824", ORANGE, "D", "-.", "Bempp, 4824 unknowns"),
    ("bem9750", RUST, "X", (0, (3, 1)), "Bempp, 9750 unknowns"),
):
    curve(ax, res["delta"], res[key], col, mk, ls, label, 4.5, 1.4, 0.6)
ax.set(xlabel="distance from the surface", ylabel="relative error against Mie")
ax.set_title(r"unit sphere, $k = 2\pi$: error near the surface", loc="left", fontsize=8.5)
ax.legend(frameon=False, fontsize=8, loc="center", bbox_to_anchor=(0.5, 0.24), ncol=2, columnspacing=0.8, handlelength=2.2)
finish(fig, ax, "fig_sphere_nearfield")
