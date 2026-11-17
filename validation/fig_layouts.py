"""Print figures of the two source layouts, drawn from the coordinates that
the solvers use.

    python3 fig_layouts.py      # writes fig_cube_layout.pdf/.png and fig_almond_layout.pdf/.png

Cube: cross-section through an edge (free lines, tail lines, base points,
collocation nodes) and the stretched variable along a line. Almond: source
rows and tip tail in the mid-plane, with the focal curve.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pec_almond_skeleton as a
import pec_cube_lightning3d as m

BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, BODY = "#1f1f1e", "#6b6a64", "#eeedea"


def style():
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})


def cube(J=6, fname="fig_cube_layout"):
    H = m.H
    depths = 0.2 * 0.5 ** np.arange(J + 1)
    tail = depths[-1] * 0.5 ** np.arange(1, 10)
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[1, 1, 1.05]))
    # (a) cross-section of the cube through the middle of an edge, corner region
    ax = axs[0]
    ax.fill([0, H, H, 0], [0, 0, H, H], color=BODY, lw=0)
    ax.plot([0, H, H], [H, H, 0], color=INK, lw=1.2)
    ax.plot([0, H], [0, H], color=MUTED, lw=0.7, ls=(0, (4, 3)))
    c = np.sqrt(0.5)
    ax.plot(H - depths * c, H - depths * c, "o", color=BLUE, ms=4, mec="white", mew=0.5, label="free lines")
    ax.plot(H - tail * c, H - tail * c, "o", mfc="white", mec=ORANGE, ms=3, mew=0.8, label="tail lines")
    s = (np.arange(6) + 0.5) / 6 * 2 * (H - 0.2) - (H - 0.2)
    s = s[s > 0]
    ax.plot(s, np.full_like(s, H - 0.2), "s", color=GREEN, ms=3.5, mec="white", mew=0.4, label="base points")
    ax.plot(np.full_like(s, H - 0.2), s, "s", color=GREEN, ms=3.5, mec="white", mew=0.4)
    nodes, _ = m.graded_nodes_capped(8, 3, 0.5, 0.15)
    nodes = nodes[nodes > 0]
    ax.plot(nodes, np.full_like(nodes, H), "|", color=INK, ms=5, mew=0.6, label="collocation")
    ax.plot(np.full_like(nodes, H), nodes, "_", color=INK, ms=5, mew=0.6)
    ax.text(0.085, 0.115, "bisector plane", fontsize=7.5, color=MUTED, rotation=45, rotation_mode="anchor")
    ax.set(xlim=(0, 0.56), ylim=(0, 0.56), aspect="equal", xlabel="$x$", ylabel="$y$")
    ax.set_title("(a) section through an edge", loc="left", fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, handletextpad=0.3,
              columnspacing=0.8)
    # (b) the same on a logarithmic scale of the depth
    ax = axs[1]
    ax.plot(depths, np.zeros_like(depths), "o", color=BLUE, ms=4.5, mec="white", mew=0.5)
    ax.plot(tail, np.zeros_like(tail), "o", mfc="white", mec=ORANGE, ms=3.5, mew=0.8)
    for j, d in enumerate(depths):
        if j in (0, J):
            ax.annotate(f"$d_{{{j if j == 0 else 'J'}}}$", (d, 0), (d, 0.16), ha="center", fontsize=8)
    ax.annotate("weights $2^{-ie}$", (tail[4], 0), (tail[4], -0.24), ha="center", fontsize=8, color=ORANGE)
    ax.annotate("free", (depths[3], 0), (depths[3], -0.24), ha="center", fontsize=8, color=BLUE)
    ax.set(xscale="log", xlim=(3e-6, 0.5), ylim=(-0.5, 0.5), yticks=[], xlabel="depth below the edge")
    ax.spines["left"].set_visible(False)
    ax.set_title("(b) depths on the bisector", loc="left", fontsize=8.5)
    # (c) stretched variable along a line
    ax = axs[2]
    u = np.linspace(-1, 1, 400)
    for d, ls, col in ((0.2, "-", GREEN), (0.2 / 16, "--", BLUE), (0.2 / 256, "-.", ORANGE)):
        lam = np.log((2 * H - d) / d)
        ax.plot(u, H * np.tanh(0.5 * lam * u), ls, color=col, lw=1.4, label=f"$d = {d:.2g}$")
    ax.set(xlabel="stretched variable $u$", ylabel="position $t$ along the edge", xlim=(-1, 1), ylim=(-0.52, 0.52))
    ax.axhline(H, color=MUTED, lw=0.5, ls=":")
    ax.axhline(-H, color=MUTED, lw=0.5, ls=":")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title("(c) parameter of a line", loc="left", fontsize=8.5)
    for ax in axs:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout(w_pad=1.0)
    fig.savefig(fname + ".pdf", bbox_inches="tight")
    fig.savefig(fname + ".png", dpi=200, bbox_inches="tight")
    print("wrote", fname + ".pdf/.png")


def almond(n_mid=6, J=3, fname="fig_almond_layout"):
    y = a.source_rows(n_mid, J, 0.5)
    tpts, row, rr = a.tip_tail(n_mid, J, 0.5)
    L = a.L
    xo = np.linspace(-a.A_ if hasattr(a, "A_") else y[:, 0].min(), L * 0.58333, 600)
    # outline of the mid-plane section: front part from ab(), rear part an ellipse
    xf = np.linspace(0, L * 0.58333, 400)
    af = a.ab(xf)[0]
    psi = np.linspace(0, np.pi / 2, 200)
    xr, ar = -a.A_ * np.sin(psi), a.B_ * np.cos(psi)
    X = np.concatenate((xr[::-1], xf)) * 1e3
    Y = np.concatenate((ar[::-1], af)) * 1e3
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.5), gridspec_kw=dict(width_ratios=[2.3, 1]))
    ax = axs[0]
    ax.fill_between(X, 0, Y, color=BODY, lw=0)
    ax.plot(X, Y, color=INK, lw=1.2)
    ax.plot(X, a.ETA_F * Y, color=ORANGE, lw=0.9, ls=(0, (4, 2)), label="focal curve")
    ax.plot(y[:, 0] * 1e3, y[:, 1] * 1e3, ".", color=BLUE, ms=2.2, label="sources")
    ax.set(xlabel="$x$ [mm]", ylabel="$y$ [mm]", aspect="equal", ylim=(0, 62))
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", markerscale=2.5)
    ax.set_title("(a) mid-plane $z = 0$, upper half", loc="left", fontsize=8.5)
    x0 = L * 0.58333 * 1e3
    ax.plot([x0 - 6, x0 - 6, x0 + 1, x0 + 1, x0 - 6], [0, 5, 5, 0, 0], color=MUTED, lw=0.6)
    ax = axs[1]
    sel = X > x0 - 6.5
    ax.fill_between(X[sel] - x0, 0, Y[sel], color=BODY, lw=0)
    ax.plot(X[sel] - x0, Y[sel], color=INK, lw=1.2)
    ax.plot(X[sel] - x0, a.ETA_F * Y[sel], color=ORANGE, lw=0.9, ls=(0, (4, 2)))
    near = y[:, 0] * 1e3 > x0 - 6.5
    ax.plot(y[near, 0] * 1e3 - x0, y[near, 1] * 1e3, ".", color=BLUE, ms=3, label="free")
    ax.plot(tpts[:, 0] * 1e3 - x0, tpts[:, 1] * 1e3, "o", mfc="white", mec=ORANGE, ms=2.2, mew=0.5, label="tail")
    ax.set(xlabel="distance from the tip [mm]", xlim=(-6, 0.3), ylim=(0, 3.4))
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", markerscale=1.6)
    ax.set_title("(b) the tip", loc="left", fontsize=8.5)
    for ax in axs:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout(w_pad=1.0)
    fig.savefig(fname + ".pdf", bbox_inches="tight")
    fig.savefig(fname + ".png", dpi=200, bbox_inches="tight")
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    style()
    cube()
    almond()
