"""Print figure of the sphere study from results_pec_sphere.json (written by
pec_sphere_convergence.py): error against sqrt(N) for four source radii, and
error, boundary residual and coefficient norm against the source radius, for
k = 2 pi and 4 pi.

    python3 pec_sphere_convergence.py && python3 fig_sphere.py
    python3 fig_sphere.py short        # k = 2 pi only, both panels in one row (fig_sphere_short.pdf, used by the paper)
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLS = ("#2a78d6", "#eb6834", "#1baf7a", "#8a5cc7")
MARK = ("o", "D", "s", "^")
DASH = ("-", "-.", "--", ":")
INK, MUTED = "#1f1f1e", "#6b6a64"

if __name__ == "__main__":
    res = json.load(open("results_pec_sphere.json"))
    ks = sorted({r["k"] for r in res})
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    short = "short" in sys.argv
    if short:
        ks = ks[:1]
        fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.9), squeeze=False)
        axs = axs.T
    else:
        fig, axs = plt.subplots(2, len(ks), figsize=(7.4, 5.6), squeeze=False)
    for col, k in enumerate(ks):
        rn = [r for r in res if r["k"] == k and r["sweep"] == "N"]
        rr = sorted((r for r in res if r["k"] == k and r["sweep"] == "R"), key=lambda r: r["R"])
        ax = axs[0, col]
        for i, R in enumerate(sorted({r["R"] for r in rn})):
            s = sorted((r for r in rn if r["R"] == R), key=lambda r: r["N"])
            ax.semilogy([np.sqrt(r["N"]) for r in s], [r["mid"] for r in s], linestyle=DASH[i % 4], color=COLS[i % 4],
                        marker=MARK[i % 4], ms=4.5, lw=1.4, mec="white", mew=0.6, label=f"$R = {R}$")
        ax.set(xlabel=r"$\sqrt{N}$", ylabel="relative error against Mie" if col == 0 else "")
        ax.set_title(f"$k = {k / np.pi:.0f}\\pi$: error against the number of sources", loc="left", fontsize=8.5)
        ax.legend(frameon=False, fontsize=7.5, loc="lower left")
        ax = axs[1, col]
        R = [r["R"] for r in rr]
        for i, (key, lab) in enumerate((("mid", "error of the field"), ("bc", "boundary residual"))):
            ax.semilogy(R, [r[key] for r in rr], linestyle=DASH[i], color=COLS[i], marker=MARK[i], ms=4.5, lw=1.4,
                        mec="white", mew=0.6, label=lab)
        ax.semilogy(R, [r["alpha_norm"] for r in rr], linestyle=":", color=MUTED, marker="s", ms=4, lw=1.2,
                    mec="white", mew=0.5, label="coefficient norm")
        ax.set(xlabel="source radius $R$")
        ax.set_title(f"$k = {k / np.pi:.0f}\\pi$, $N = {rr[0]['N']}$: dependence on the source radius", loc="left", fontsize=8.5)
        ax.legend(frameon=False, fontsize=7.5, loc="upper center")
    for ax in axs.ravel():
        ax.grid(True, which="major", color="#e4e3df", lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    name = "fig_sphere_short" if short else "fig_sphere"
    fig.savefig(name + ".pdf")
    fig.savefig(name + ".png", dpi=200)
    print(f"wrote {name}.pdf/.png")
