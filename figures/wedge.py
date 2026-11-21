"""
Figure 4: cross-section normal to an edge of the cube. Shows how analytic
continuations across the two faces meet on the skeleton wall to create jump
currents, discretised by geometric free + tail line sources matched by graded
boundary collocation panels on the two faces.

    python figures/wedge.py
"""
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch
from style import BLUE, BODY, INK, MUTED, ORANGE, finish, plt

fig, ax = plt.subplots(figsize=(3.55, 3.25))
L = 1.02
ax.fill([-L, 0, 0, -L], [-L, -L, 0, 0], color=BODY, lw=0, zorder=0)
ax.plot([0, 0], [-L, 0], color=INK, lw=1.45, zorder=3)                         # face phi = 0
ax.plot([-L, 0], [0, 0], color=INK, lw=1.45, zorder=3)                         # face phi = 3 pi / 2
ax.plot([0, 0.76], [0, 0.76], color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)   # exterior bisector, psi = 0
ax.plot([0, -0.94], [0, -0.94], color=BLUE, lw=1.35, alpha=0.55, zorder=2)     # wall of the skeleton

# Dyadically graded boundary panels + 3 Gauss nodes per panel on both faces
g3 = 0.5 * (1.0 + np.polynomial.legendre.leggauss(3)[0])
edges_b = np.array([0.0, 0.045, 0.09, 0.18, 0.36, 0.72, 1.00])
for a, b in zip(edges_b[:-1], edges_b[1:]):
    if b >= 0.09:
        ax.plot([0, 0.032], [-b, -b], color=MUTED, lw=0.65, zorder=3)
        ax.plot([-b, -b], [0, 0.032], color=MUTED, lw=0.65, zorder=3)
    pts_g = a + (b - a) * g3
    if b >= 0.09:
        ax.plot(np.full_like(pts_g, 0.015), -pts_g, ".", color=INK, ms=3.1, zorder=4)
        ax.plot(-pts_g, np.full_like(pts_g, 0.015), ".", color=INK, ms=3.1, zorder=4)

# Geometric free depths (2 per octave) and singular tail depths (1 per octave) on skeleton wall
c = 1.0 / np.sqrt(2.0)
d_free = 0.82 * (2.0 ** -0.5) ** np.arange(0, 6)
d_tail = d_free[-1] * (0.5 ** np.arange(1, 5))
# Strip dividers showing strip width proportional to depth
div_d = np.sqrt(d_free[:-1] * d_free[1:])
for dv in div_d[:3]:
    px, py = -dv * c, -dv * c
    ax.plot([px - 0.032, px + 0.032], [py + 0.032, py - 0.032], color=BLUE, lw=0.6, alpha=0.65, zorder=2)

ax.plot(-d_free * c, -d_free * c, "o", color=BLUE, ms=4.8, mec="white", mew=0.55, zorder=5)
ax.plot(-d_tail * c, -d_tail * c, "o", mfc="white", mec=ORANGE, ms=4.0, mew=0.95, zorder=5)
ax.plot(0, 0, "o", color=INK, ms=4.0, zorder=6)

# Exterior point (rho, phi, psi)
ang = np.deg2rad(118.0)
P = 0.74 * np.array([np.cos(ang), np.sin(ang)])
ax.plot([0, P[0]], [0, P[1]], color=INK, lw=0.8)
ax.plot(*P, "o", color=INK, ms=3.0)
ax.text(0.5 * P[0] - 0.07, 0.5 * P[1] + 0.02, r"$\rho$", ha="center", va="center", fontsize=8.5)
ax.add_patch(Arc((0, 0), 0.44, 0.44, theta1=-90, theta2=118, color=ORANGE, lw=1.0))
ax.text(0.27, -0.13, r"$\phi$", color=ORANGE, ha="center", va="center", fontsize=8.5)
ax.add_patch(Arc((0, 0), 0.90, 0.90, theta1=45, theta2=118, color=MUTED, lw=1.0))
ax.text(0.19, 0.52, r"$\psi$", color=MUTED, ha="center", va="center", fontsize=8.5)
ax.text(0.78, 0.78, r"$\psi = 0$", color=MUTED, ha="left", va="bottom", fontsize=8.2)
ax.text(0.06, -0.95, r"$\phi = 0$  (graded $\Gamma$)", ha="left", va="center", fontsize=7.8, color=INK)
ax.text(-1.00, 0.07, r"$\phi = 3\pi/2$", ha="left", va="bottom", fontsize=8.0, color=INK)

ax.text(-0.08, -0.84, "skeleton $\\mathcal{M}$ ($\\psi = \\pm\\pi$)\nfree ($\\bullet$) + tail ($\\circ$)",
        color=BLUE, ha="right", va="center", fontsize=7.3, zorder=6,
        bbox=dict(boxstyle="round,pad=0.08", fc=BODY, ec="none"))
ax.text(-0.99, -0.45, "conductor", color=MUTED, ha="left", va="center", fontsize=8.0)
ax.text(0.46, 0.16, "exterior", color=MUTED, ha="center", va="center", fontsize=8.2)

# Analytic continuations across both faces meeting on M
for (x0, y0), (x1, y1), rad in (((0.16, -0.34), (-0.23, -0.35), 0.28), ((-0.34, 0.16), (-0.35, -0.23), -0.28)):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                                 mutation_scale=8, color=INK, lw=0.8))

ax.set(xlim=(-1.05, 1.12), ylim=(-1.06, 0.98), aspect="equal")
ax.axis("off")
finish(fig, None, "fig_wedge", tight="tight")
