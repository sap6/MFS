"""
Figure 4: cross-section normal to an edge of the cube (a sketch). The field is
continued into the conductor through each of the two faces, and the two
continuations meet on the wall of the skeleton.

    python figures/wedge.py
"""
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch
from style import BLUE, BODY, INK, MUTED, ORANGE, finish, plt

fig, ax = plt.subplots(figsize=(3.3, 3.1))
L = 1.0
ax.fill([-L, 0, 0, -L], [-L, -L, 0, 0], color=BODY, lw=0)                       # the conductor: the quadrant x, y < 0
ax.plot([0, 0], [-L, 0], color=INK, lw=1.4)                                     # face phi = 0
ax.plot([-L, 0], [0, 0], color=INK, lw=1.4)                                     # face phi = 3 pi / 2
ax.plot([0, 0.78], [0, 0.78], color=MUTED, lw=0.8, ls=(0, (4, 3)))              # exterior bisector, psi = 0
ax.plot([0, -0.82], [0, -0.82], color=BLUE, lw=2.0)                             # wall of the skeleton
ax.plot(0, 0, "o", color=INK, ms=4)
# A point of the exterior, its radius and the two angles
ang = np.deg2rad(118.0)
P = 0.74 * np.array([np.cos(ang), np.sin(ang)])
ax.plot([0, P[0]], [0, P[1]], color=INK, lw=0.8)
ax.plot(*P, "o", color=INK, ms=3)
ax.text(0.5 * P[0] - 0.07, 0.5 * P[1] + 0.02, r"$\rho$", ha="center", va="center")
ax.add_patch(Arc((0, 0), 0.44, 0.44, theta1=-90, theta2=118, color=ORANGE, lw=1.0))
ax.text(0.27, -0.13, r"$\phi$", color=ORANGE, ha="center", va="center")
ax.add_patch(Arc((0, 0), 0.92, 0.92, theta1=45, theta2=118, color=MUTED, lw=1.0))
ax.text(0.19, 0.53, r"$\psi$", color=MUTED, ha="center", va="center")
ax.text(0.80, 0.80, r"$\psi = 0$", color=MUTED, ha="left", va="bottom")
ax.text(0.04, -0.93, r"$\phi = 0$", ha="left", va="center")
ax.text(-0.96, 0.05, r"$\phi = 3\pi/2$", ha="left", va="bottom")
ax.text(-0.56, -0.81, "wall of the\nskeleton,\n" + r"$\psi = \pm\pi$", color=BLUE, ha="left", va="center", fontsize=8.5)
ax.text(0.45, 0.15, "exterior", color=MUTED, ha="center", va="center", fontsize=8.5)
ax.text(-0.74, -0.42, "conductor", color=MUTED, ha="center", va="center", fontsize=8.5)
# The two continuations into the conductor
for (x0, y0), (x1, y1), rad in (((0.15, -0.32), (-0.24, -0.33), 0.3), ((-0.32, 0.15), (-0.33, -0.24), -0.3)):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>", mutation_scale=8,
                                 color=INK, lw=0.8))
ax.set(xlim=(-1.02, 1.12), ylim=(-1.02, 1.0), aspect="equal")
ax.axis("off")
finish(fig, None, "fig_wedge", tight="tight")
