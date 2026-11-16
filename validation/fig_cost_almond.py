"""Error against wall time for the NASA almond at 1.19 GHz, nose-on incidence:
MFS with skeleton sources against Bempp (EFIE, RWG, dense LU).

    python3 fig_cost_almond.py run    # solve and time both on this machine; writes fig_cost_almond.npz
    python3 fig_cost_almond.py        # table and figure from the saved data

Times include assembly, solve and evaluation at 400 points on |x| = 0.3 m. The
MFS uses the two mirror symmetries of this incidence; Bempp does not. The
error is measured against the finest MFS layout. Run alone on the machine.
"""
import sys
import time

import numpy as np

import pec_almond_skeleton as a
from pec_sphere_convergence import fibonacci_sphere

FREQ = 1.19
LAYOUTS = ((3, 0), (4, 1), (5, 2), (6, 3), (8, 5))        # (n_mid, J); the last is the reference
BEM = ((28, 26), (40, 36), (60, 54))
STYLE = {"MFS": ("#2a78d6", "o", "-", "MFS, skeleton sources"), "Bempp": ("#eb6834", "D", "-.", "Bempp (EFIE, RWG)")}
INK, MUTED = "#1f1f1e", "#6b6a64"


def run():
    k = 2 * np.pi * FREQ * 1e9 / 299792458.0
    pts = 0.3 * fibonacci_sphere(400)
    a.solve(k, a.source_rows(2, 0, 0.5))                   # warm-up
    out, rows = {}, []
    for n_mid, J in LAYOUTS:
        t0 = time.time()
        y = a.source_rows(n_mid, J, 0.5)
        c, info = a.solve(k, y)
        E = (a.matrix(pts, y, k) @ c).reshape(-1, 3)
        t = time.time() - t0
        ext, sca = a.cross_sections(k, y, c)
        rows.append((len(y), info["unknowns"], t, abs(ext - sca) / sca))
        out[f"mfs_E{len(rows) - 1}"] = E
        print(f"MFS n_mid={n_mid} J={J}: sources={len(y)} unknowns={info['unknowns']} rows={info['rows']} "
              f"time={t:.1f}s optical={rows[-1][3]:.1e}", flush=True)
    out["mfs"] = np.array(rows)
    a.bem_reference(k, pts, ((10, 10),))                   # warm-up: Bempp compiles its kernels on first use
    brows = []
    for i, (ndof, Eb, t) in enumerate(a.bem_reference(k, pts, BEM)):
        brows.append((ndof, t))
        out[f"bem_E{i}"] = Eb
        print(f"Bempp dofs={ndof} time={t:.1f}s", flush=True)
    out["bem"] = np.array(brows)
    np.savez("fig_cost_almond.npz", **out)


def rel(x, y):
    return np.linalg.norm(x - y) / np.linalg.norm(y)


def report():
    d = np.load("fig_cost_almond.npz")
    ref = d[f"mfs_E{len(d['mfs']) - 1}"]
    curves = {"MFS": [], "Bempp": []}
    print("method   unknowns   time[s]   error vs finest MFS   optical")
    for i, (n, unk, t, opt) in enumerate(d["mfs"][:-1]):
        e = rel(d[f"mfs_E{i}"], ref)
        curves["MFS"].append((t, e))
        print(f"MFS    {int(unk):9d} {t:9.1f} {e:14.2e} {opt:16.1e}")
    for i, (ndof, t) in enumerate(d["bem"]):
        e = rel(d[f"bem_E{i}"], ref)
        curves["Bempp"].append((t, e))
        print(f"Bempp  {int(ndof):9d} {t:9.1f} {e:14.2e}")
    return {k_: np.array(v) for k_, v in curves.items()}


def plot(curves, fname="fig_cost_almond"):
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
    ax.set_ylabel("relative error of the field on |x| = 0.3 m")
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_title("NASA almond, 1.19 GHz, nose-on", loc="left", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(fname + ".pdf")
    fig.savefig(fname + ".png", dpi=200)
    print("wrote", fname + ".pdf/.png")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    plot(report())
