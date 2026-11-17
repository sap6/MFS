"""One figure with the four cost comparisons (sphere, torus, almond, cube):
error against wall time for the MFS and for Bempp, from the saved data of
the fig_cost_*.py scripts. An alternative to the four separate figures and
the three small tables, if the paper has to be shorter.

    python3 fig_cost_all.py      # writes fig_cost_all.pdf/.png
"""
import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, NAVY, GREEN, ORANGE, PINK, PURPLE = "#2a78d6", "#184f96", "#1baf7a", "#eb6834", "#e87ba4", "#8a5cc7"
INK, MUTED = "#1f1f1e", "#6b6a64"


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def curve(ax, t, e, col, mk, ls, label):
    ax.loglog(t, e, linestyle=ls, color=col, marker=mk, ms=5, lw=1.5, mec="white", mew=0.7, label=label)


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.2))
    ax = axs[0, 0]
    d = np.load("fig_cost_sphere.npz")
    curve(ax, d["mfs"][:, 1], d["mfs"][:, 2], BLUE, "o", "-", "MFS")
    curve(ax, d["bem"][:, 1], d["bem"][:, 2], ORANGE, "D", "-.", "Bempp")
    ax.set_title("(a) sphere, two wavelengths across (against Mie)", loc="left", fontsize=8.5)
    ax = axs[0, 1]
    d = np.load("fig_cost_torus.npz")
    curve(ax, d["mfs"][:, 1], d["mfs"][:, 2], BLUE, "o", "-", "MFS, mode by mode")
    curve(ax, d["dense"][:, 1], d["dense"][:, 2], GREEN, "s", "--", "MFS, full matrix")
    curve(ax, d["bem"][:, 1], d["bem"][:, 2], ORANGE, "D", "-.", "Bempp")
    ax.set_title("(b) torus", loc="left", fontsize=8.5)
    ax = axs[1, 0]
    d = np.load("fig_cost_almond.npz")
    ref = d[f"mfs_E{len(d['mfs']) - 1}"]
    curve(ax, d["mfs"][:-1, 2], [rel(d[f"mfs_E{i}"], ref) for i in range(len(d["mfs"]) - 1)], BLUE, "o", "-", "MFS")
    curve(ax, d["bem"][:, 1], [rel(d[f"bem_E{i}"], ref) for i in range(len(d["bem"]))], ORANGE, "D", "-.", "Bempp")
    ax.set_title("(c) almond, 1.19 GHz, nose-on", loc="left", fontsize=8.5)
    ax = axs[1, 1]
    d = np.load("fig_cost_cube.npz")
    ref = d["dense_E10"]
    for key, col, mk, ls, lab in (("tail", BLUE, "o", "-", "MFS, one depth per octave"), ("dense", NAVY, "h", (0, (6, 1.5)), "MFS, two depths per octave")):
        rows = [r for r in d[key] if not (key == "dense" and int(r[0]) == 10)]
        curve(ax, [r[1] for r in rows], [rel(d[f"{key}_E{int(r[0])}"], ref) for r in rows], col, mk, ls, lab)
    for tag, col, mk, ls, lab in (("_q4", ORANGE, "D", "-.", "Bempp, quasi-uniform"),):
        bem = sorted((np.load(f) for f in glob.glob("bem_cube_h*_q4.npz")), key=lambda b: int(b["ndof"]))
        bem = [b for b in bem if int(b["ndof"]) > 1000]
        curve(ax, [float(b["t_solve"]) for b in bem], [rel(b["E"], ref) for b in bem], col, mk, ls, lab)
    for tag, col, mk, ls, lab in (("", PINK, "v", (0, (5, 1, 1, 1)), "Bempp, graded, quadrature 4/6"), ("_q6_8", PURPLE, "P", (0, (1, 1)), "Bempp, graded, quadrature 6/8")):
        gr = sorted((np.load(f) for f in glob.glob(f"bem_cube_graded_n*_g2{tag}.npz")), key=lambda b: int(b["ndof"]))
        curve(ax, [float(b["t_solve"]) for b in gr], [rel(b["E"], ref) for b in gr], col, mk, ls, lab)
    best = []
    for tag in ("", "_q6_8", "_q8_10"):
        best += [(float(b["t_solve"]), rel(b["E"], ref)) for b in (np.load(g) for g in glob.glob(f"bem_cube_sym_n*_g2{tag}.npz"))]
    best.sort()
    env, emin = [], np.inf
    for t, e in best:                                   # lower envelope over the quadrature settings
        if e < 0.97 * emin:
            env.append((t, e))
            emin = e
    if env:
        env = np.array(env)
        curve(ax, env[:, 0], env[:, 1], "#9c3b12", "X", (0, (3, 1)), "Bempp, graded, with mirror symmetry")
    ax.set_title("(d) cube, side one wavelength", loc="left", fontsize=8.5)
    for ax in axs.ravel():
        ax.grid(True, which="major", color="#e4e3df", lw=0.6)
        ax.legend(frameon=False, fontsize=7, loc="lower left")
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
