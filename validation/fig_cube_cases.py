"""Cube beyond the standard case: error against wall time for a general
incident wave and for a side of two wavelengths, from the data of
cube_cases.py (MFS and the best Bempp runs with the mirror symmetry, timed
alternately).

    python3 fig_cube_cases.py      # writes fig_cube_cases.pdf/.png
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

from fig_cost_all import BLUE, INK, MUTED, curve

PANELS = (("k1o", "(a) side one wavelength, general wave"), ("k2", "(b) side two wavelengths, standard wave"))

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0), sharey=True)
    for ax, (case, title) in zip(axs, PANELS):
        if os.path.exists(f"cube_cases_{case}.npz"):
            d = np.load(f"cube_cases_{case}.npz")
            curve(ax, d["mfs"][:, 5], d["mfs"][:, 6], BLUE, "o", "-", "MFS")
            env = d["env"][d["env"][:, 6] < 2e-2]
            curve(ax, env[:, 5], env[:, 6], "#9c3b12", "X", (0, (3, 1)), "Bempp, graded, with symmetry")
        ax.set_title(title, loc="left", fontsize=8.5)
        ax.grid(True, which="major", color="#e4e3df", lw=0.6)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("wall time [s]  (4 cores)")
    axs[0].set_ylabel("relative error of the scattered field")
    axs[0].legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=2.6)
    fig.tight_layout()
    fig.savefig("fig_cube_cases.pdf")
    fig.savefig("fig_cube_cases.png", dpi=200)
    print("wrote fig_cube_cases.pdf/.png")
