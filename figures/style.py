"""The look of the figures of the paper, and where their data live"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from mfs.core import rel  # noqa: E402, F401

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
BLUE, NAVY, ORANGE, RUST, GREEN, PURPLE = "#2a78d6", "#184f96", "#eb6834", "#9c3b12", "#1baf7a", "#8a5cc7"
INK, MUTED, GRID, BODY = "#1f1f1e", "#6b6a64", "#e4e3df", "#eeedea"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})


def data(name):
    return os.path.join(DATA, name)


def envelope(points):
    """The (time, error) runs that are at least 3% more accurate than every faster run"""
    out, best = [], np.inf
    for t, e in sorted(points):
        if e < 0.97 * best:
            out.append((t, e))
            best = e
    return np.array(out)


def curve(ax, x, y, col, mk, ls, label, ms=5, lw=1.5):
    ax.loglog(x, y, linestyle=ls, color=col, marker=mk, ms=ms, lw=lw, mec="white", mew=0.7, label=label)


def finish(fig, axs, name, rect=(0, 0, 1, 1), grid=True, tight=None):
    """Grid, no top and right spines; writes name.pdf and name.png next to this file"""
    for ax in np.ravel(axs):
        if grid:
            ax.grid(True, which="major", color=GRID, lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout(rect=rect)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), dpi=200, bbox_inches=tight)
    print(f"wrote figures/{name}.pdf")
