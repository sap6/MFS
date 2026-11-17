"""Unit sphere, k = 2 pi: error of the scattered field against the Mie series
as a function of the distance from the surface, for the MFS (N = 400 dipoles
on the inner sphere R_s = 0.5) and for Bempp (EFIE, RWG, h = 0.07, default
quadrature). The MFS field is a finite sum of dipole fields and is evaluated
in the same way at every distance; the boundary element potential becomes
nearly singular within about one element of the surface. NOTES.md item 69.

    python3 sphere_nearfield.py run    # solves both and writes sphere_nearfield.npz
    python3 sphere_nearfield.py        # table and fig_sphere_nearfield.pdf/.png
"""
import sys

import numpy as np

import mie
from fig_cost_sphere import K, mfs
from pec_sphere_convergence import D, P, fibonacci_sphere, rel

DELTA = np.logspace(-4, 0, 17)
NPTS = 200


def run():
    import bem_env  # noqa: F401
    import bempp_cl.api as bempp
    import bem_reference as br
    shell = fibonacci_sphere(NPTS)
    pts = np.concatenate([(1 + d) * shell for d in DELTA])
    exact = mie.scattered_field(pts, K, D, P)
    out = {"delta": DELTA}
    for N in (400, 1200):
        E, t = mfs(N, pts)
        out[f"mfs{N}"] = np.array([rel(E[i * NPTS:(i + 1) * NPTS], exact[i * NPTS:(i + 1) * NPTS]) for i in range(len(DELTA))])
        print(f"MFS N = {N}: " + " ".join(f"{e:.1e}" for e in out[f"mfs{N}"]), flush=True)
    for h in (0.1, 0.07):
        E, ndof = br.scattered_field(bempp.shapes.sphere(h=h), K, pts)
        out[f"bem{ndof}"] = np.array([rel(E[i * NPTS:(i + 1) * NPTS], exact[i * NPTS:(i + 1) * NPTS]) for i in range(len(DELTA))])
        print(f"Bempp {ndof} dofs (h = {h}): " + " ".join(f"{e:.1e}" for e in out[f"bem{ndof}"]), flush=True)
    np.savez("sphere_nearfield.npz", **out)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    d = np.load("sphere_nearfield.npz")
    print("distance   " + "   ".join(k for k in d if k != "delta"))
    for i, dl in enumerate(d["delta"]):
        print(f"{dl:8.1e}   " + "   ".join(f"{d[k][i]:.1e}" for k in d if k != "delta"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from fig_cost_all import BLUE, NAVY, ORANGE, INK, MUTED
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(figsize=(5.0, 3.3))
    keys = [k for k in d if k != "delta"]
    style = {"mfs": ((BLUE, "o", "-"), (NAVY, "h", (0, (6, 1.5)))), "bem": ((ORANGE, "D", "-."), ("#9c3b12", "X", (0, (3, 1))))}
    count = {"mfs": 0, "bem": 0}
    for k in keys:
        kind = k[:3]
        col, mk, ls = style[kind][count[kind]]
        count[kind] += 1
        label = f"MFS, $N = {int(k[3:])}$" if kind == "mfs" else f"Bempp, {int(k[3:])} unknowns"
        ax.loglog(d["delta"], d[k], linestyle=ls, color=col, marker=mk, ms=4.5, lw=1.4, mec="white", mew=0.6, label=label)
    ax.set_xlabel("distance from the surface")
    ax.set_ylabel("relative error against Mie")
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.36))
    fig.tight_layout()
    fig.savefig("fig_sphere_nearfield.pdf")
    fig.savefig("fig_sphere_nearfield.png", dpi=200)
    print("wrote fig_sphere_nearfield.pdf/.png")
