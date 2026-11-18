"""Error against wall time for the unit sphere (k = 2 pi), against the Mie
series: MFS with dipoles on the inner sphere R = 0.5, and Bempp (EFIE, RWG).

    python3 fig_cost_sphere.py run    # solve and time both on this machine; writes fig_cost_sphere.npz
    python3 fig_cost_sphere.py        # table and figure from the saved data

Times include assembly, solve and evaluation at 400 points on |x| = 2. No
symmetry is used by either method. Run alone on the machine.
"""
import sys
import time

import numpy as np
from lsq import qr_solve

import mie
from pec_sphere_convergence import D, P, dyadic, evaluate, fibonacci_sphere, rel, tangent_basis

K = 2 * np.pi
NS = (50, 100, 150, 200, 300, 400, 600, 800, 1200)
HS = (0.2, 0.14, 0.1, 0.07, 0.056)     # 0.056: about 15500 unknowns, the largest full matrix under the 4 GB of the OpenCL driver
STYLE = {"MFS": ("#2a78d6", "o", "-", "MFS, dipoles on an inner sphere"), "Bempp": ("#eb6834", "D", "-.", "Bempp (EFIE, RWG)")}
INK, MUTED = "#1f1f1e", "#6b6a64"


def mfs(N, pts, R=0.5):
    t0 = time.time()
    x, y = fibonacci_sphere(2 * N), R * fibonacci_sphere(N)
    T = tangent_basis(x)
    A = np.einsum("mab,mnbc->manc", T, dyadic(x, y, K)).reshape(4 * N, 3 * N)
    b = -np.einsum("mab,mb->ma", T, mie.plane_wave(x, K, D, P)).reshape(-1)
    alpha = qr_solve(np.asfortranarray(A), b, cond=1e-14)[0].reshape(N, 3)   # same accuracy as the SVD solver, 4 times faster
    E = evaluate(pts, y, alpha, K)
    return E, time.time() - t0


def run():
    import bem_env  # noqa: F401  (memory for the OpenCL driver; must come before bempp)
    import bempp_cl.api as bempp
    import bem_reference as br
    pts = 2.0 * fibonacci_sphere(400)
    exact = mie.scattered_field(pts, K, D, P)
    mfs(50, pts)
    rows = []
    for N in NS:
        E, t = mfs(N, pts)
        rows.append((3 * N, t, rel(E, exact)))
        print(f"MFS   unknowns={3 * N:5d} time={t:7.2f}s error={rows[-1][2]:.2e}", flush=True)
    out = {"mfs": np.array(rows)}
    br.scattered_field(bempp.shapes.sphere(h=0.4), K, pts)            # warm-up (compilation)
    rows = []
    print("Bempp assembles with", bempp.DEFAULT_DEVICE_INTERFACE, bempp.DEFAULT_PRECISION, flush=True)
    for h in HS:
        t0, split = time.time(), []
        E, ndof = br.scattered_field(bempp.shapes.sphere(h=h), K, pts, split)
        rows.append((ndof, time.time() - t0, rel(E, exact)))
        print(f"Bempp unknowns={ndof:5d} time={rows[-1][1]:7.2f}s error={rows[-1][2]:.2e}  (assembly, LU, evaluation: "
              + ", ".join(f"{v:.1f}" for v in split[0]) + " s)", flush=True)
    out["bem"] = np.array(rows)
    np.savez("fig_cost_sphere.npz", **out)


def report():
    d = np.load("fig_cost_sphere.npz")
    print("method   unknowns   time[s]   error vs Mie")
    for name, key in (("MFS", "mfs"), ("Bempp", "bem")):
        for unk, t, e in d[key]:
            print(f"{name:6s} {int(unk):9d} {t:9.2f} {e:12.2e}")
    return {"MFS": d["mfs"][:, 1:], "Bempp": d["bem"][:, 1:]}


def plot(curves, fname="fig_cost_sphere"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for key, (col, mk, ls, label) in STYLE.items():
        c = curves[key]
        ax.loglog(c[:, 0], c[:, 1], ls, color=col, marker=mk, ms=6, lw=1.6, mec="white", mew=0.8, label=label)
        ax.annotate(key, (c[-1, 0] * 1.2, c[-1, 1]), ha="left", va="center", fontsize=8, color=INK)
    ax.set_xlabel("wall time [s]  (4 cores)")
    ax.set_ylabel("relative error of the field on |x| = 2 (Mie)")
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_title("Unit sphere, diameter = two wavelengths", loc="left", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(fname + ".pdf")
    fig.savefig(fname + ".png", dpi=200)
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    plot(report())
