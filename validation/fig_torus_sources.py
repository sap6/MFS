"""Torus: how small may the source torus be?

Boundary residual and coefficient norm |alpha| against the minor radius a_s of
a concentric source torus, for three tube radii a/R = 0.25, 0.5, 0.8 (k = 2 pi,
wave along the axis), with the same number of source rings and dipoles per
ring throughout. The vertical lines mark a_s = R - sqrt(R^2 - a^2), the
distance from the centre circle of the tube to the focal ring, where the
continued field is singular (DERIVATIONS.md, E): below it the coefficients
grow without bound.

    python3 fig_torus_sources.py run    # compute; writes fig_torus_sources.npz
    python3 fig_torus_sources.py        # table and figure from the saved data
"""
import sys

import numpy as np

from pec_torus import Torus, scatter, wave

K = 2 * np.pi
TUBES = (0.25, 0.5, 0.8)
FRACS = (0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05)
NS, NPHI = 40, 160
STYLE = {0.25: ("#2a78d6", "o", "-"), 0.5: ("#eb6834", "D", "-."), 0.8: ("#1baf7a", "s", "--")}
INK, MUTED = "#1f1f1e", "#6b6a64"


def run():
    d, p = wave()
    out = {}
    for a in TUBES:
        tor = Torus(1.0, a)
        rows = []
        for f in FRACS:
            s = scatter(tor, K, d, p, NS, NPHI, a_s=f * a, ntheta_col=3 * NS)
            i = s.info
            rows.append((f, i["bc"], i["alpha_norm"], i["opt"]))
            print(f"a/R={a:4.2f} a_s/a={f:4.2f} residual={i['bc']:.1e} |alpha|={i['alpha_norm']:.1e} optical={i['opt']:.1e}", flush=True)
        out[f"a{a:g}"] = np.array(rows)
    np.savez("fig_torus_sources.npz", **out)


def report():
    d = np.load("fig_torus_sources.npz")
    for a in TUBES:
        thr = (1 - np.sqrt(1 - a * a)) / a
        print(f"a/R = {a}: focal-ring threshold a_s/a = {thr:.3f}")
        for f, bc, al, op in d[f"a{a:g}"]:
            print(f"   a_s/a={f:4.2f}{'  ' if f > thr else ' *'} residual={bc:.1e} |alpha|={al:.1e} optical={op:.1e}")
    return d


def plot(d, fname="fig_torus_sources"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    for a in TUBES:
        col, mk, ls = STYLE[a]
        r = d[f"a{a:g}"]
        thr = (1 - np.sqrt(1 - a * a)) / a
        for ax, j in zip(axes, (1, 2)):
            ax.semilogy(r[:, 0], r[:, j], color=col, marker=mk, ls=ls, lw=1.5, ms=5, label=f"a/R = {a:g}")
            ax.axvline(thr, color=col, lw=0.9, ls=":")
    axes[0].set_ylabel("boundary residual", color=INK)
    axes[1].set_ylabel(r"coefficient norm $\|\alpha\|$", color=INK)
    for ax in axes:
        ax.set_xlabel(r"source torus radius $a_s/a$", color=INK)
        ax.grid(True, which="major", color="#d9d8d2", lw=0.6)
        ax.tick_params(colors=MUTED)
        for s in ax.spines.values():
            s.set_color(MUTED)
    axes[0].legend(frameon=False, fontsize=8.5)
    axes[1].set_title(r"dotted: $a_s = R-\sqrt{R^2-a^2}$ (focal ring)", fontsize=9, color=MUTED)
    axes[0].set_title(f"{3 * NS * NPHI} unknowns, k = 2$\\pi$", fontsize=9, color=MUTED)
    fig.tight_layout()
    fig.savefig(fname + ".pdf")
    fig.savefig(fname + ".png", dpi=200)
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    if "run" in sys.argv:
        run()
    plot(report())
