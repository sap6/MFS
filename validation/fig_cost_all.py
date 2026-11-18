"""One figure with the four cost comparisons (sphere, torus, almond, cube):
error against wall time for the MFS and for Bempp, from the saved data of
the fig_cost_*.py scripts. An alternative to the four separate figures and
the three small tables, if the paper has to be shorter.

    python3 fig_cost_all.py      # writes fig_cost_all.pdf/.png
"""
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

BLUE, NAVY, GREEN, ORANGE, PINK, PURPLE = "#2a78d6", "#184f96", "#1baf7a", "#eb6834", "#e87ba4", "#8a5cc7"
INK, MUTED = "#1f1f1e", "#6b6a64"


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def lower_envelope(points):
    """The runs that no faster run beats, as an array (time, error)."""
    env, emin = [], np.inf
    for t, e in sorted(points):
        if e < 0.97 * emin:
            env.append((t, e))
            emin = e
    return np.array(env)


def curve(ax, t, e, col, mk, ls, label):
    ax.loglog(t, e, linestyle=ls, color=col, marker=mk, ms=5, lw=1.5, mec="white", mew=0.7, label=label)


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.4))
    ax = axs[0, 0]
    d = np.load("fig_cost_sphere.npz")
    curve(ax, d["mfs"][:, 1], d["mfs"][:, 2], BLUE, "o", "-", "MFS")
    curve(ax, d["bem"][:, 1], d["bem"][:, 2], ORANGE, "D", "-.", "Bempp")
    ax.set_title("(a) sphere, two wavelengths across (against Mie)", loc="left", fontsize=8.5)
    ax = axs[0, 1]
    d = np.load("fig_cost_torus.npz")
    curve(ax, d["mfs"][:, 1], d["mfs"][:, 2], BLUE, "o", "-", "MFS")
    curve(ax, d["bem"][:, 1], d["bem"][:, 2], ORANGE, "D", "-.", "Bempp, with symmetry")      # both use the rotational symmetry
    ax.set_title("(b) torus", loc="left", fontsize=8.5)
    ax = axs[1, 0]
    d = np.load("fig_cost_almond.npz")
    ref = d[f"mfs_E{len(d['mfs']) - 1}"]
    curve(ax, d["mfs"][:-1, 2], [rel(d[f"mfs_E{i}"], ref) for i in range(len(d["mfs"]) - 1)], BLUE, "o", "-", "MFS")
    curve(ax, d["bem"][:, 1], [rel(d[f"bem_E{i}"], ref) for i in range(len(d["bem"]))], ORANGE, "D", "-.", "Bempp, with symmetry")
    ax.set_title("(c) almond, 1.19 GHz, nose-on", loc="left", fontsize=8.5)
    ax.set_xticks([2, 5, 10, 20])
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, pos: f"{v:g}"))
    ax = axs[1, 1]
    d = np.load("fig_cost_cube.npz")
    ref = d["ref_E"]                                    # the finest level, J = 12
    jref = int(d["ref_info"][0])
    for key, col, mk, ls, lab in (("dense", BLUE, "o", "-", "MFS"), ("tail", NAVY, "h", (0, (6, 1.5)), "MFS, one depth per octave")):
        rows = [r for r in d[key] if int(r[0]) >= 4 and not (key == "dense" and int(r[0]) == jref)]
        curve(ax, [r[1] for r in rows], [rel(d[f"{key}_E{int(r[0])}"], ref) for r in rows], col, mk, ls, lab)
    bem = sorted((np.load(f) for f in glob.glob("bem_cube_h*_q4.npz")), key=lambda b: int(b["ndof"]))
    bem = [b for b in bem if int(b["ndof"]) > 1000]
    curve(ax, [float(b["t_solve"]) for b in bem], [rel(b["E"], ref) for b in bem], ORANGE, "D", "-.", "Bempp, quasi-uniform")
    if os.path.exists("bem_grading_scan_k1_final.npz"):  # with the symmetry: best grading and quadrature at each time, re-timed
        env = np.load("bem_grading_scan_k1_final.npz")["rows"][:, [5, 6]]
    elif os.path.exists("bem_grading_scan_k1.npz"):
        env = lower_envelope([(r[5], r[6]) for r in np.load("bem_grading_scan_k1.npz")["rows"]])
    else:
        best = []
        for tag in ("", "_q6_8", "_q8_10"):
            best += [(float(b["t_solve"]), rel(b["E"], ref)) for b in (np.load(g) for g in glob.glob(f"bem_cube_sym_n*_g2{tag}.npz"))]
        env = lower_envelope(best)
    curve(ax, env[:, 0], env[:, 1], "#9c3b12", "X", (0, (3, 1)), "Bempp, graded, with symmetry")
    ax.set_title("(d) cube, side one wavelength", loc="left", fontsize=8.5)
    ax.set_ylim(8e-7, 3e-2)
    for ax in axs.ravel():
        ax.grid(True, which="major", color="#e4e3df", lw=0.6)
        ax.legend(frameon=False, fontsize=6.8, loc="upper right" if ax is axs[1, 0] else "lower left", handlelength=2.6)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axs[1]:
        ax.set_xlabel("wall time [s]  (4 cores)")
    for ax in axs[:, 0]:
        ax.set_ylabel("relative error of the scattered field")
    fig.tight_layout()
    fig.savefig("fig_cost_all.pdf")
    fig.savefig("fig_cost_all.png", dpi=200)
    print("wrote fig_cost_all.pdf/.png")
