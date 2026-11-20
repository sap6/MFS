"""
Figure 8: the regular tetrahedron of unit edge length at k = 2π.
  (a) The six triangular walls of the skeleton, with edge lines on one wall and corner rays.
  (b) Error of the scattered field against √N_u (N_u the unknowns of the whole tetrahedron),
      measured against the finest layout J = 10 and checked against graded symmetric Bempp.

    python figures/tetrahedron.py          # draw from data/tetrahedron_convergence.npz
    python figures/tetrahedron.py --run    # compute MFS and Bempp runs (about 6 minutes)
"""
import gc
import subprocess
import sys
import time

if "--run" in sys.argv:
    subprocess.check_call([sys.executable, __file__, "--run-mfs"])
    for n_edge in (16, 24, 32, 44):
        subprocess.check_call([sys.executable, __file__, "--bem-one", str(n_edge)])
    sys.argv.remove("--run")

if "--bem-one" in sys.argv:
    # Graded symmetric Bempp-cl (quadrature orders 6 and 8, grading g = 1.8), one fresh process per mesh
    import os
    import numpy as np
    from mfs import bem, core
    bem.bempp.GLOBAL_PARAMETERS.quadrature.regular, bem.bempp.GLOBAL_PARAMETERS.quadrature.singular = 6, 8
    n_edge = int(sys.argv[sys.argv.index("--bem-one") + 1])
    npz_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "tetrahedron_convergence.npz")
    out = dict(np.load(npz_path))
    pts = out["pts"]
    d_w, p_w = core.wave()
    bem.tetrahedron(6, 1.8, 2 * np.pi, d_w, p_w, pts[:4])              # warm-up OpenCL kernels
    t0 = time.time()
    Eb, Nb = bem.tetrahedron(n_edge, 1.8, 2 * np.pi, d_w, p_w, pts)
    dt = time.time() - t0
    err_b = core.rel(Eb, out["ref_E"])
    out[f"bem_E{n_edge}"] = Eb
    bem_rows = [r for r in out.get("bem", np.zeros((0, 4))).tolist() if int(r[0]) != n_edge]
    bem_rows.append([n_edge, Nb, dt, err_b])
    out["bem"] = np.array(sorted(bem_rows))
    np.savez(npz_path, **out)
    print(f"bempp  n={n_edge:2d}: N_u={Nb:5d}, t={dt:5.1f}s, vs_J10={err_b:.2e}", flush=True)
    sys.exit(0)

import numpy as np
from matplotlib.patches import Polygon
from style import BLUE, GREEN, INK, MUTED, NAVY, ORANGE, PURPLE, RUST, curve, data, finish, plt, rel

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
            print(f"{name:6s} J={J:2d}: unk/cls={info['unknowns'][0]:4d} (N_u={4*info['unknowns'][0]:5d}), "
                  f"t={t:5.1f}s, opt={info['optical']:.2e}, ext={info['ext']:.6f}", flush=True)
        out[name] = np.array(rows)
    out["ref_E"] = out["dense_E10"]
    out["bem"] = np.zeros((0, 4))
    np.savez(data("tetrahedron_convergence.npz"), **out)
    sys.exit(0)

# Load stored results
res = np.load(data("tetrahedron_convergence.npz"))
ref = res["ref_E"]
specs = {"dense": lambda J: dict(J=J), "one": lambda J: dict(J=J, q=0.5),
         "notail": lambda J: dict(J=J, tails=False), "plain": lambda J: dict(J=None, n_base=J)}
curves = {name: np.array([(sum(tetrahedron.Tetrahedron(classes=[], **specs[name](int(r[0]))).unknowns.values()),
                           rel(res[f"{name}_E{int(r[0])}"], ref)) for r in res[name] if int(r[0]) < 10 or name == "plain"])
          for name in specs}
curves["graded"] = np.column_stack((res["bem"][:, 1], res["bem"][:, 3]))

# Root-exponential fit to levels J = 3..8 (J = 9 approaches the J = 10 reference)
c = curves["dense"][:-1]
p = np.polyfit(np.sqrt(c[:, 0]), np.log(c[:, 1]), 1)
print(f"tetrahedron fit: error = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N_u)), max factor "
      f"{np.exp(np.abs(np.log(c[:, 1]) - np.polyval(p, np.sqrt(c[:, 0]))).max()):.2f}")

# Two-panel figure: (a) 3D skeleton + sources on one wall, (b) convergence against sqrt(N_u)
fig, axs = plt.subplots(1, 2, figsize=(6.7, 3.1), gridspec_kw=dict(width_ratios=[0.95, 1.38]))

# Panel (a): orthographic view of the regular tetrahedron and its six triangular skeleton walls
ax = axs[0]
az, el = np.radians(-34.0), np.radians(22.0)
R_az = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
R_el = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
P2 = (R_el @ R_az)[:2]                                                # 2 x 3 orthographic projection
V = tetrahedron.VERTS @ P2.T
O = np.zeros(2)
edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
for i, j in edges:
    col = "#cfe0f5" if (i, j) == (0, 1) else "#eae9e4"
    al = 0.65 if (i, j) == (0, 1) else 0.30
    ax.add_patch(Polygon([O, V[i], V[j]], closed=True, facecolor=col, edgecolor="none", alpha=al, zorder=1))
for i in range(4):
    ax.plot([0, V[i, 0]], [0, V[i, 1]], ls=(0, (3, 2)), lw=0.75, color=MUTED, zorder=2)

# Edge lines and corner rays on the highlighted wall (0, v_0, v_1)
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

# Tetrahedron wireframe (back edge (2, 3) dashed)
for i, j in edges:
    ls_e = (0, (3, 2)) if (i, j) == (2, 3) else "-"
    lw_e = 0.8 if (i, j) == (2, 3) else 1.15
    ax.plot([V[i, 0], V[j, 0]], [V[i, 1], V[j, 1]], ls=ls_e, lw=lw_e, color=INK, zorder=5)
ax.plot(0, 0, "o", color=MUTED, ms=3.0, zorder=5)
ax.text(V[0, 0] + 0.02, V[0, 1] + 0.02, "corner\nrays", fontsize=7.5, color=ORANGE, va="bottom")
ax.text(0.5 * (V[0, 0] + V[1, 0]) + 0.03, 0.5 * (V[0, 1] + V[1, 1]) - 0.02, "edge\nlines", fontsize=7.5, color=BLUE, va="top")
ax.text(-0.02, -0.05, r"$\mathcal{M}$", fontsize=8.2, color=MUTED, ha="right", va="top")
ax.set(xlim=(-0.54, 0.64), ylim=(-0.52, 0.56), aspect="equal")
ax.axis("off")
ax.set_title("(a) skeleton and wall sources", loc="left", fontsize=8.5)

# Panel (b): convergence against sqrt(N_u)
ax = axs[1]
xs = np.array([np.sqrt(curves["dense"][0, 0]) - 4, np.sqrt(curves["dense"][-1, 0]) + 6])
ax.semilogy(xs, np.exp(np.polyval(p, xs)), color=MUTED, lw=0.8, zorder=1)
ax.annotate(rf"${np.exp(p[1]):.2f}\,\exp({p[0]:.4f}\sqrt{{N_u}})$",
            (99.0, 3.4e-3), fontsize=8.0, color=MUTED)
for key, (col, mk, ls, label) in {"plain": (PURPLE, "^", ":", "plain MFS (inner tetrahedron)"),
                                  "notail": (GREEN, "s", "--", "MFS without the tails"),
                                  "one": (NAVY, "h", (0, (6, 1.5)), "MFS (one depth per octave)"),
                                  "dense": (BLUE, "o", "-", "MFS (two depths per octave)"),
                                  "graded": (RUST, "X", (0, (3, 1)), "Bempp, graded meshes")}.items():
    v = curves[key]
    curve(ax, np.sqrt(v[:, 0]), v[:, 1], col, mk, ls, label, lw=1.4, logx=False, zorder=3)
ax.set(xlabel=r"$\sqrt{N_u}$,  $N_u$ = unknowns of the whole tetrahedron", ylabel="relative error of the scattered field")
ax.set_title(r"(b) regular tetrahedron, $k = 2\pi$: error against $\sqrt{N_u}$", loc="left", fontsize=8.5)
ax.legend(frameon=False, fontsize=7.6, loc="upper center", bbox_to_anchor=(0.46, -0.28), ncol=2, columnspacing=1.0)
ax.grid(True, which="major", color="#e4e3df", lw=0.6)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.24, wspace=0.26)
finish(fig, None, "fig_tetrahedron", tight="tight")
