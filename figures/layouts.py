"""
Figures 3 and 8: the source layouts of the cube and of the almond.

    python figures/layouts.py
"""
import numpy as np
from style import BLUE, BODY, GREEN, INK, MUTED, ORANGE, finish, plt

from mfs import almond, cube


def cube_layout(J=6, fname="fig_cube_layout"):
    """Figure 3: the source layout of the cube at level J"""
    H, layout = cube.H, cube.Cube(J, classes=[])                   # the layout of the paper, not assembled
    depths = layout.depths
    tail = depths[-1] * 0.5 ** np.arange(1, 10)
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.9), gridspec_kw=dict(width_ratios=[1, 1, 1.05]))
    # (a) Cross-section of the cube through the middle of an edge, corner region
    ax = axs[0]
    ax.fill([0, H, H, 0], [0, 0, H, H], color=BODY, lw=0)
    ax.plot([0, H, H], [H, H, 0], color=INK, lw=1.2)
    ax.plot([0, H], [0, H], color=MUTED, lw=0.7, ls=(0, (4, 3)))
    c = np.sqrt(0.5)
    ax.plot(H - depths * c, H - depths * c, "o", color=BLUE, ms=4, mec="white", mew=0.5, label="free lines")
    ax.plot(H - tail * c, H - tail * c, "o", mfc="white", mec=ORANGE, ms=3, mew=0.8, label="tail lines")
    s = (np.arange(8) + 0.5) / 8 * 2 * (H - 0.2) - (H - 0.2)        # n_base = 8
    s = s[s > 0]
    ax.plot(s, np.full_like(s, H - 0.2), "s", color=GREEN, ms=3.5, mec="white", mew=0.4, label="base points")
    ax.plot(np.full_like(s, H - 0.2), s, "s", color=GREEN, ms=3.5, mew=0.4, mec="white")
    nodes = np.unique(layout.x[layout.x[:, 0] < H, 0])            # the boundary points along one face
    ax.plot(nodes, np.full_like(nodes, H), "|", color=INK, ms=5, mew=0.6, label="boundary points")
    ax.plot(np.full_like(nodes, H), nodes, "_", color=INK, ms=5, mew=0.6)
    ax.text(0.055, 0.085, "bisector plane", fontsize=8, color=MUTED, rotation=45, rotation_mode="anchor")
    ax.set(xlim=(0, 0.56), ylim=(0, 0.56), aspect="equal", xlabel="$x$", ylabel="$y$")
    ax.set_title("(a) section through an edge", loc="left", fontsize=8.5)
    ax.legend(frameon=False, fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2, handletextpad=0.3,
              columnspacing=0.7)
    # (b) The same on a logarithmic scale of the depth
    ax = axs[1]
    ax.plot(depths, np.zeros_like(depths), "o", color=BLUE, ms=4.5, mec="white", mew=0.5)
    ax.plot(tail, np.zeros_like(tail), "o", mfc="white", mec=ORANGE, ms=3.5, mew=0.8)
    for d, lab in ((depths[0], "$d_0$"), (depths[-1], r"$d_{\min}$")):
        ax.annotate(lab, (d, 0), (d, 0.16), ha="center", fontsize=8)
    ax.annotate("tail, fixed weights", (tail[5], 0), (tail[5], -0.24), ha="center", fontsize=8, color=ORANGE)
    ax.annotate("free", (depths[4], 0), (depths[4], -0.24), ha="center", fontsize=8, color=BLUE)
    ax.set(xscale="log", xlim=(3e-6, 0.5), ylim=(-0.5, 0.5), yticks=[], xlabel="depth below the edge")
    ax.set_box_aspect(1)
    ax.spines["left"].set_visible(False)
    ax.set_title("(b) depths on the bisector", loc="left", fontsize=8.5)
    # (c) The stretched variable along a line
    ax = axs[2]
    u = np.linspace(-1, 1, 400)
    for d, ls, col, lab in ((0.2, "-", GREEN, "$d = 0.2$"), (0.2 / 16, "--", BLUE, "$d = 0.2/16$"), (0.2 / 256, "-.", ORANGE, "$d = 0.2/256$")):
        lam = np.log((2 * H - d) / d)
        ax.plot(u, H * np.tanh(0.5 * lam * u), ls, color=col, lw=1.4, label=lab)
    ax.set(xlabel="stretched variable $u$", ylabel="edge position $t$", xlim=(-1, 1), ylim=(-0.55, 0.55))
    ax.set_box_aspect(0.95)
    ax.axhline(H, color=MUTED, lw=0.5, ls=":")
    ax.axhline(-H, color=MUTED, lw=0.5, ls=":")
    ax.legend(frameon=False, fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2, handlelength=1.6,
              handletextpad=0.35, columnspacing=0.7)
    ax.set_title("(c) parameter of a line", loc="left", fontsize=8.5)
    finish(fig, axs, fname, grid=False, tight="tight")


def almond_layout(n_mid=6, n_F=3, fname="fig_almond_layout"):
    """Figure 8: the source layout of the almond, n_mid:n_F = 6:3"""
    y, (tpts, *_) = almond.source_rows(n_mid, n_F)
    L = almond.L
    # Outline of the mid-plane section: the front half from semi_axes(), the rear half an ellipse
    xf = np.linspace(0, L * almond.T1, 400)
    af = almond.semi_axes(xf)[0]
    psi = np.linspace(0, np.pi / 2, 200)
    xr, ar = -almond.A_ * np.sin(psi), almond.B_ * np.cos(psi)
    X = np.concatenate((xr[::-1], xf)) * 1e3
    Y = np.concatenate((ar[::-1], af)) * 1e3
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.3), gridspec_kw=dict(width_ratios=[2.2, 1]))
    ax = axs[0]
    ax.fill_between(X, 0, Y, color=BODY, lw=0)
    ax.plot(X, Y, color=INK, lw=1.2)
    ax.plot(X, almond.ETA_F * Y, color=ORANGE, lw=0.9, ls=(0, (4, 2)), label="focal curve")
    ax.plot(y[:, 0] * 1e3, y[:, 1] * 1e3, ".", color=BLUE, ms=2.2, label="sources")
    ax.set(xlabel="$x$ [mm]", ylabel="$y$ [mm]", aspect="equal", ylim=(0, 62))
    ax.legend(frameon=False, fontsize=8, loc="upper right", markerscale=2.5)
    ax.set_title("(a) mid-plane $z = 0$, upper half", loc="left", fontsize=8.5)
    x0 = L * almond.T1 * 1e3
    ax.plot([x0 - 6, x0 - 6, x0 + 1, x0 + 1, x0 - 6], [0, 5, 5, 0, 0], color=MUTED, lw=0.6)
    ax = axs[1]
    sel = X > x0 - 6.5
    ax.fill_between(X[sel] - x0, 0, Y[sel], color=BODY, lw=0)
    ax.plot(X[sel] - x0, Y[sel], color=INK, lw=1.2)
    ax.plot(X[sel] - x0, almond.ETA_F * Y[sel], color=ORANGE, lw=0.9, ls=(0, (4, 2)))
    near = y[:, 0] * 1e3 > x0 - 6.5
    ax.plot(y[near, 0] * 1e3 - x0, y[near, 1] * 1e3, ".", color=BLUE, ms=3, label="free")
    ax.plot(tpts[:, 0] * 1e3 - x0, tpts[:, 1] * 1e3, "o", mfc="white", mec=ORANGE, ms=2.2, mew=0.5, label="tail")
    ax.set(xlabel="offset from the tip [mm]", ylabel="$y$ [mm]", xlim=(-6, 0.3), ylim=(0, 3.4))
    ax.legend(frameon=False, fontsize=8, loc="upper right", markerscale=1.6)
    ax.set_title("(b) the tip", loc="left", fontsize=8.5)
    finish(fig, axs, fname, grid=False, tight="tight")


cube_layout()
almond_layout()
