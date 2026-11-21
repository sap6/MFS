"""
Figure 8: the regular tetrahedron of unit edge length at k = 2π.
  (a) The six triangular walls of the skeleton, with edge lines on one wall and corner rays.
  (b) Error of the scattered field against √N_u (N_u the unknowns of the whole tetrahedron),
      measured against the finest layout J = 10 and checked against graded symmetric Bempp.
  (c) Error against wall time on two cores, comparing MFS (two and one depths per octave) with
      the lower envelope of symmetric graded Bempp across quadrature orders 6/8 and 8/10.

    python figures/tetrahedron.py          # draw from data/tetrahedron_convergence.npz
    python figures/tetrahedron.py --run    # compute MFS and Bempp runs (about 10 minutes)
"""
import gc
import subprocess
import sys
import time

if "--run" in sys.argv:
    subprocess.check_call([sys.executable, __file__, "--run-mfs"])
    for n_edge, q_reg in ((16, 6), (24, 6), (32, 6), (16, 8), (24, 8), (32, 8), (44, 8)):
        subprocess.check_call([sys.executable, __file__, "--bem-one", str(n_edge), str(q_reg)])
    sys.argv.remove("--run")

if "--bem-one" in sys.argv:
    # Graded symmetric Bempp-cl (quadrature q_reg / q_reg+2, grading g = 1.8), one fresh process per run
    import os
    import numpy as np
    from mfs import bem, core
    idx_arg = sys.argv.index("--bem-one")
    n_edge = int(sys.argv[idx_arg + 1])
    q_reg = int(sys.argv[idx_arg + 2]) if idx_arg + 2 < len(sys.argv) else 8
    bem.bempp.GLOBAL_PARAMETERS.quadrature.regular, bem.bempp.GLOBAL_PARAMETERS.quadrature.singular = q_reg, q_reg + 2
    npz_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "tetrahedron_convergence.npz")
    out = dict(np.load(npz_path))
    pts = out["pts"]
    d_w, p_w = core.wave()
    bem.tetrahedron(6, 1.8, 2 * np.pi, d_w, p_w, pts[:4])              # warm-up OpenCL kernels
    t0 = time.time()
    Eb, Nb = bem.tetrahedron(n_edge, 1.8, 2 * np.pi, d_w, p_w, pts)
    dt = time.time() - t0
    err_b = core.rel(Eb, out["ref_E"])
    out[f"bem_E{n_edge}_q{q_reg}"] = Eb
    if q_reg == 8:
        out[f"bem_E{n_edge}"] = Eb
    old_bem = out.get("bem", np.zeros((0, 5)))
    if old_bem.ndim == 2 and old_bem.shape[1] == 4:
        old_bem = np.column_stack((old_bem[:, 0], np.full(len(old_bem), 8.0), old_bem[:, 1:]))
    bem_rows = [r for r in old_bem.tolist() if not (int(r[0]) == n_edge and int(r[1]) == q_reg)]
    bem_rows.append([n_edge, q_reg, Nb, dt, err_b])
    out["bem"] = np.array(sorted(bem_rows))
    np.savez(npz_path, **out)
    print(f"bempp  n={n_edge:2d} q={q_reg}/{q_reg+2}: N_u={Nb:5d}, t={dt:5.1f}s, vs_J10={err_b:.2e}", flush=True)
    sys.exit(0)

import numpy as np
from matplotlib.patches import Polygon
from style import BLUE, GREEN, INK, MUTED, NAVY, ORANGE, PURPLE, RUST, curve, data, envelope, finish, plt, rel

from mfs import core, tetrahedron

if "--run-mfs" in sys.argv:
    pts = 2.0 * core.fibonacci(400)
    tetrahedron.solve(2, pts=pts)                                      # warm-up: compilation is not timed
    out = {"pts": pts}
    for name, Js, kw in (("dense", range(3, 11), {}), ("one", range(3, 9), dict(q=0.5)),
                         ("notail", (4, 6, 8), dict(tails=False)), ("plain", (4, 6, 8, 10), None)):
        rows = []
        for J in Js:
            gc.collect()
            E, t, info = (tetrahedron.solve(None, pts=pts, n_base=J) if kw is None
                          else tetrahedron.solve(J, pts=pts, **kw))
            rows.append((J, t, info["unknowns"][0], info["rows"], info["optical"], info["ext"]))
            out[f"{name}_E{J}"] = E
            print(f"{name:6s} J={J:2d}: unk/cls={info['unknowns'][0]:4d}, "
                  f"t={t:5.1f}s, opt={info['optical']:.2e}, ext={info['ext']:.6f}", flush=True)
        out[name] = np.array(rows)
    out["ref_E"] = out["dense_E10"]
    out["bem"] = np.zeros((0, 5))
    np.savez(data("tetrahedron_convergence.npz"), **out)
    sys.exit(0)

# Load stored results
res = np.load(data("tetrahedron_convergence.npz"))
ref = res["ref_E"]
specs = {"dense": lambda J: dict(J=J), "one": lambda J: dict(J=J, q=0.5),
         "notail": lambda J: dict(J=J, tails=False), "plain": lambda J: dict(J=None, n_base=J)}
curves = {name: np.array([(sum(tetrahedron.Tetrahedron(classes=[], **specs[name](int(r[0]))).unknowns.values()),
                           rel(res[f"{name}_E{int(r[0])}"], ref), float(r[1]))
                          for r in res[name] if int(r[0]) < 10 or name == "plain"])
          for name in specs}
bem_arr = res["bem"]
if bem_arr.shape[1] == 4:
    bem_arr = np.column_stack((bem_arr[:, 0], np.full(len(bem_arr), 8.0), bem_arr[:, 1:]))
best_by_N = {}
for _, _, Nb, _, eb in bem_arr:
    best_by_N[int(Nb)] = min(best_by_N.get(int(Nb), 1.0), float(eb))
curves["graded"] = np.array(sorted(best_by_N.items()))
bem_time_env = envelope([(float(r[3]), float(r[4])) for r in bem_arr])

# Root-exponential fit to levels J = 3..8 (J = 9 approaches the J = 10 reference)
c = curves["dense"][:-1]
p = np.polyfit(np.sqrt(c[:, 0]), np.log(c[:, 1]), 1)
print(f"tetrahedron fit: error = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N_u)), max factor "
      f"{np.exp(np.abs(np.log(c[:, 1]) - np.polyval(p, np.sqrt(c[:, 0]))).max()):.2f}")

# Three-panel figure: (a) 3D skeleton + sources, (b) error vs sqrt(N_u), (c) error vs wall time [s]
fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.15), gridspec_kw=dict(width_ratios=[0.92, 1.16, 1.12]))

# Panel (a): orthographic view of the regular tetrahedron and its six triangular skeleton walls
ax = axs[0]
az, el = np.radians(-34.0), np.radians(22.0)
R_az = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
R_el = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
P2 = (R_el @ R_az)[:2]
V = tetrahedron.VERTS @ P2.T
O = np.zeros(2)
edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
for i, j in edges:
    col = "#cfe0f5" if (i, j) == (0, 1) else "#eae9e4"
    al = 0.65 if (i, j) == (0, 1) else 0.30
    ax.add_patch(Polygon([O, V[i], V[j]], closed=True, facecolor=col, edgecolor="none", alpha=al, zorder=1))
for i in range(4):
    ax.plot([0, V[i, 0]], [0, V[i, 1]], ls=(0, (3, 2)), lw=0.75, color=MUTED, zorder=2)

t_hat = np.array([0.0, 1.0, 1.0]) / np.sqrt(2.0)
b_hat = np.array([-1.0, 0.0, 0.0])
mid0 = np.array([tetrahedron.A_TET, 0.0, 0.0])
for d_ln in 0.10 * (2 ** -0.5) ** np.arange(0, 9, 2):
    ts = np.linspace(-(0.5 - 1.5 * d_ln), 0.5 - 1.5 * d_ln, 40)
    pts3 = mid0[None, :] + d_ln * b_hat[None, :] + ts[:, None] * t_hat[None, :]
    xy = pts3 @ P2.T
    ax.plot(xy[:, 0], xy[:, 1], color=BLUE, lw=0.95, zorder=3)
v0 = tetrahedron.VERTS[0]
for rd in [np.ones(3) / np.sqrt(3.0)] + [(t_hat - (fr / np.sqrt(2.0)) * b_hat) / np.linalg.norm(t_hat - (fr / np.sqrt(2.0)) * b_hat)
                                         for fr in (0.8, 0.4, 0.2)]:
    seg = np.vstack((v0 - 0.08 * rd, v0 - 0.008 * rd)) @ P2.T
    ax.plot(seg[:, 0], seg[:, 1], color=ORANGE, lw=0.95, zorder=4)

for i, j in edges:
    ls_e = (0, (3, 2)) if (i, j) == (2, 3) else "-"
    lw_e = 0.8 if (i, j) == (2, 3) else 1.15
    ax.plot([V[i, 0], V[j, 0]], [V[i, 1], V[j, 1]], ls=ls_e, lw=lw_e, color=INK, zorder=5)
ax.plot(0, 0, "o", color=MUTED, ms=3.0, zorder=5)
ax.text(V[0, 0] - 0.01, V[0, 1] + 0.05, "corner\nrays", fontsize=7.5, color=ORANGE, ha="center", va="bottom")
ax.text(0.5 * (V[0, 0] + V[1, 0]) + 0.04, 0.5 * (V[0, 1] + V[1, 1]) - 0.03, "edge\nlines", fontsize=7.5, color=BLUE, va="top")
ax.text(-0.02, -0.05, r"$\mathcal{M}$", fontsize=8.2, color=MUTED, ha="right", va="top")
ax.set(xlim=(-0.56, 0.60), ylim=(-0.48, 0.63), aspect="equal")
ax.axis("off")
ax.set_title("(a) skeleton and wall sources", loc="left", fontsize=8.5)

# Panel (b): convergence against sqrt(N_u)
ax = axs[1]
xs = np.array([np.sqrt(curves["dense"][0, 0]) - 4, np.sqrt(curves["dense"][-1, 0]) + 6])
ax.semilogy(xs, np.exp(np.polyval(p, xs)), color=MUTED, lw=0.85, zorder=1)
ax.annotate(rf"${np.exp(p[1]):.2f}\,\exp({p[0]:.4f}\sqrt{{N_u}})$",
            xy=(121.0, np.exp(np.polyval(p, 121.0))), xytext=(85.0, 6.8e-3),
            fontsize=7.4, color=MUTED, ha="left", va="bottom",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.55, shrinkA=1, shrinkB=2))
for key, (col, mk, ls, label) in {"plain": (PURPLE, "^", ":", "plain MFS (inner copy)"),
                                  "notail": (GREEN, "s", "--", "MFS without tails"),
                                  "one": (NAVY, "h", (0, (6, 1.5)), "MFS (1 depth/octave)"),
                                  "dense": (BLUE, "o", "-", "MFS (2 depths/octave)"),
                                  "graded": (RUST, "X", (0, (3, 1)), "Bempp (graded, symmetry)")}.items():
    v = curves[key]
    curve(ax, np.sqrt(v[:, 0]), v[:, 1], col, mk, ls, label, ms=4.4, lw=1.35, logx=False, zorder=3)
ax.set(xlabel=r"$\sqrt{N_u}$ (whole-body unknowns)", ylabel="relative error of the scattered field",
       ylim=(4e-5, 1.0))
ax.set_title(r"(b) error against $\sqrt{N_u}$", loc="left", fontsize=8.5)

# Panel (c): error against wall time [s]
ax = axs[2]
curve(ax, curves["one"][:-1, 2], curves["one"][:-1, 1], NAVY, "h", (0, (6, 1.5)), ms=4.4, lw=1.35, zorder=3)
curve(ax, curves["dense"][:, 2], curves["dense"][:, 1], BLUE, "o", "-", ms=4.4, lw=1.35, zorder=4)
curve(ax, bem_time_env[:, 0], bem_time_env[:, 1], RUST, "X", (0, (3, 1)), ms=4.4, lw=1.35, zorder=3)
ax.set(xlabel="wall time [s]  (2 cores)", xlim=(2.5, 135), ylim=(4e-5, 1.0))
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_title("(c) error against wall time", loc="left", fontsize=8.5)

handles, labels = axs[1].get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, fontsize=7.8, loc="lower center", bbox_to_anchor=(0.56, -0.02),
           ncol=3, columnspacing=1.2, handletextpad=0.45)
for a in (axs[1], axs[2]):
    a.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.28, wspace=0.28)
finish(fig, None, "fig_tetrahedron", tight="tight")

