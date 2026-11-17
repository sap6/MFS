"""Convergence law on the cube: error against the square root of the number
of unknowns, for the MFS layouts and for Bempp (NOTES.md item 59).

    python3 fig_cube_convergence.py

N is the number of unknowns of the whole cube (for the MFS all eight parity
classes together, for Bempp the degrees of freedom of the mesh). Bempp: for
each graded mesh the smallest error over all gradings and quadrature orders
of bem_grading_scan_k1.npz (if present), else the stored g = 2 runs.
"""
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fig_cost_cube import INK, MUTED, STYLE


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def data():
    d = np.load("fig_cost_cube.npz")
    jref = int(d["dense"][-1][0])
    ref = d[f"dense_E{jref}"]
    out = {}
    for name, key in (("tail", "MFS+tail"), ("dense", "MFS dense"), ("notail", "MFS no tail"), ("plain", "plain MFS")):
        rows = [(4 * r[2], rel(d[f"{name}_E{int(r[0])}"], ref)) for r in d[name] if not (name == "dense" and int(r[0]) == jref)]
        out[key] = np.array(rows)
    rows = sorted((int(b["ndof"]), rel(b["E"], ref)) for b in (np.load(f) for f in glob.glob("bem_cube_h*_q4.npz")) if int(b["ndof"]) > 1000)
    out["Bempp"] = np.array(rows)
    best = {}
    if os.path.exists("bem_grading_scan_k1.npz"):
        for g, reg, sing, n, N, t, err in np.load("bem_grading_scan_k1.npz")["rows"]:
            best[int(N)] = min(best.get(int(N), 1.0), err)
    else:
        for f in glob.glob("bem_cube_sym_n*_g2*.npz"):
            b = np.load(f)
            if "_obl" not in f:
                best[int(b["ndof"])] = min(best.get(int(b["ndof"]), 1.0), rel(b["E"], ref))
    out["Bempp sym"] = np.array(sorted(best.items()))
    return out


def fit(c, lo=1):
    x, y = np.sqrt(c[lo:, 0]), np.log(c[lo:, 1])
    p = np.polyfit(x, y, 1)
    return p, np.exp(np.abs(y - np.polyval(p, x)).max())


if __name__ == "__main__":
    c = data()
    p, dev = fit(c["MFS dense"])
    print(f"two depths per octave, J >= 4: error = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N)), largest deviation a factor {dev:.2f}")
    for key in ("MFS dense", "MFS+tail", "Bempp sym"):
        v = c[key]
        print(f"{key}: local order -dlog(err)/dlog(N): " + " ".join(f"{-np.log(v[i + 1, 1] / v[i, 1]) / np.log(v[i + 1, 0] / v[i, 0]):.1f}" for i in range(len(v) - 1)))
        print("     N:", " ".join(str(int(n)) for n in v[:, 0]), " err:", " ".join(f"{e:.2e}" for e in v[:, 1]))
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    labels = {"MFS+tail": "MFS with tails, one depth per octave", "MFS dense": "MFS with tails, two depths per octave",
              "MFS no tail": "MFS, edge sources, no tail", "plain MFS": "plain MFS (inner cube)",
              "Bempp": "Bempp, quasi-uniform meshes", "Bempp sym": "Bempp, graded meshes (best setting per mesh)"}
    xs = np.array([np.sqrt(c["MFS dense"][1, 0]) - 4, np.sqrt(c["MFS dense"][-1, 0]) + 8])
    ax.semilogy(xs, np.exp(np.polyval(p, xs)), color=MUTED, lw=0.8, ls="-", zorder=1)
    ax.annotate(rf"${np.exp(p[1]):.2f}\,\exp({p[0]:.3f}\sqrt{{N}})$", (xs[1] - 42, np.exp(np.polyval(p, xs[1] - 42)) * 3.2), fontsize=8, color=MUTED, ha="left")
    for key in ("plain MFS", "MFS no tail", "MFS+tail", "MFS dense", "Bempp", "Bempp sym"):
        col, mk, ls, _ = STYLE[key]
        v = c[key]
        ax.semilogy(np.sqrt(v[:, 0]), v[:, 1], linestyle=ls, color=col, marker=mk, ms=5, lw=1.4, mec="white", mew=0.7, label=labels[key], zorder=3)
    ax.set_xlabel(r"$\sqrt{N}$,  $N$ = number of unknowns")
    ax.set_ylabel("relative error of the scattered field")
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, columnspacing=1.2)
    fig.tight_layout()
    fig.savefig("fig_cube_convergence.pdf", bbox_inches="tight")
    fig.savefig("fig_cube_convergence.png", dpi=200, bbox_inches="tight")
    print("wrote fig_cube_convergence.pdf/.png")
