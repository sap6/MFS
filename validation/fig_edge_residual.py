"""Boundary residual near an edge of the cube, with and without the singular
tail: |tangential (E_sc + E_inc)| at points of the face y = 1/2 as a function
of their distance to the edge x = y = 1/2, away from the vertices (rms over
five positions along the edge and the two tangential components).

    python3 fig_edge_residual.py [J]     # writes fig_edge_residual.pdf/.png and .npz
"""
import sys

import numpy as np

import cube_mfs as cm
import pec_cube_lightning3d as m

BLUE, GREEN, INK, MUTED = "#2a78d6", "#1baf7a", "#1f1f1e", "#6b6a64"


def residual(J, tail, rho, zs):
    cube = cm.Cube(2 * J, classes=cm.classes_of(m.D, m.P), tail=tail, **cm.dense(J))     # the layout of the paper (two depths per octave)
    coef, res = cube.solve(m.D, m.P)
    x = np.array([(m.H - r, m.H, z) for z in zs for r in rho])
    T = np.tile(np.array([[1.0, 0, 0], [0, 0, 1.0]]), (len(x), 1, 1))
    A = cube.matrix(x, T, list(coef))
    E = sum(A[i][:, :cube.unknowns[c]] @ coef[c] for i, c in enumerate(coef)).reshape(-1, 2)
    inc = np.einsum("mab,mb->ma", T, m.incident(x))
    r = np.sqrt(np.mean((np.abs(E + inc) ** 2).reshape(len(zs), len(rho), 2), axis=(0, 2)))
    s = np.sqrt(np.mean((np.abs(E) ** 2).reshape(len(zs), len(rho), 2), axis=(0, 2)))
    return r, s, cube.depths[-1]


if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    rho = np.logspace(-6, np.log10(0.3), 120)
    zs = (0.03, 0.08, 0.13, 0.18, 0.23)
    r1, s1, dJ = residual(J, 9, rho, zs)
    r0, s0, _ = residual(J, 0, rho, zs)
    np.savez("fig_edge_residual.npz", rho=rho, tail=r1, notail=r0, field=s1, dJ=dJ, J=J)
    for lo, hi in ((1e-6, 1e-5), (1e-5, 1e-4), (1e-4, 1e-3), (1e-3, 1e-2), (1e-2, 1e-1), (1e-1, 0.3)):
        sel = (rho >= lo) & (rho < hi)
        print(f"rho in [{lo:.0e}, {hi:.0e}): residual with tail {np.sqrt(np.mean(r1[sel] ** 2)):.2e}, without {np.sqrt(np.mean(r0[sel] ** 2)):.2e}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.loglog(rho, r0, "--", color=GREEN, lw=1.6, label="edge sources, no tail")
    ax.loglog(rho, r1, "-", color=BLUE, lw=1.6, label="edge sources and tail")
    ax.axvline(dJ, color=MUTED, lw=0.7, ls=":")
    ax.text(dJ * 1.15, ax.get_ylim()[0] * 2.5, r"$d_{\min}$", color=MUTED, fontsize=8)
    ax.set(xlabel="distance to the edge", ylabel="tangential residual  $|\\mathbf{t}\\cdot(\\mathbf{E}_N+\\mathbf{E}^{\\rm inc})|$")
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_title(f"Unit cube, $J = {J}$: residual on a face near an edge", loc="left", fontsize=9)
    fig.tight_layout()
    fig.savefig("fig_edge_residual.pdf")
    fig.savefig("fig_edge_residual.png", dpi=200)
    print("wrote fig_edge_residual.pdf/.png")
