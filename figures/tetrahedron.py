"""
Figure 8: the regular tetrahedron of unit edge length at k = 2π.
  (a) 3D view showing the three-kite graded boundary collocation on F_0 and the interior skeleton M.
  (b) Unfolded skeleton wall triangle O v_0 v_1 showing the inner-tetrahedron cut, lightning edge lines
      with vertex-clustered t(u) nodes, 9-line edge tail, vertex dipoles, and corner-tail rays.
  (c) Error of the scattered field against √N_u (N_u the unknowns of the whole tetrahedron),
      measured against the finest layout J = 10 and checked against graded symmetric Bempp.
  (d) Error against wall time on two cores, comparing MFS (two and one depths per octave) with
      the lower envelope of symmetric graded Bempp across quadrature orders 6/8 and 8/10.

    python figures/tetrahedron.py          # draw from data/tetrahedron_convergence.npz
    python figures/tetrahedron.py --run    # compute MFS and Bempp runs
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
    # Graded symmetric Bempp-cl (quadrature q_reg / q_reg+2, grading g = 1.8), one fresh process per run;
    # keeps the minimum wall time across repeated alternating runs
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
    prev = [r for r in old_bem.tolist() if int(r[0]) == n_edge and int(r[1]) == q_reg]
    best_dt = min(dt, float(prev[0][3])) if prev else dt
    bem_rows = [r for r in old_bem.tolist() if not (int(r[0]) == n_edge and int(r[1]) == q_reg)]
    bem_rows.append([n_edge, q_reg, Nb, best_dt, err_b])
    out["bem"] = np.array(sorted(bem_rows))
    np.savez(npz_path, **out)
    print(f"bempp  n={n_edge:2d} q={q_reg}/{q_reg+2}: N_u={Nb:5d}, t={dt:5.2f}s (best={best_dt:5.2f}s), vs_J10={err_b:.3e}", flush=True)
    sys.exit(0)

import os
import numpy as np
from matplotlib.patches import ConnectionPatch, Polygon
from style import BLUE, GREEN, INK, MUTED, NAVY, ORANGE, PURPLE, RUST, curve, data, envelope, finish, plt, rel

from mfs import core, lines, tetrahedron

if "--run-mfs" in sys.argv or "--retime-mfs" in sys.argv:
    retime = "--retime-mfs" in sys.argv
    npz_path = data("tetrahedron_convergence.npz")
    old = dict(np.load(npz_path)) if os.path.exists(npz_path) else {}
    pts = old["pts"] if "pts" in old else 2.0 * core.fibonacci(400)
    tetrahedron.solve(2, pts=pts)                                      # warm-up: compilation is not timed
    out = dict(old) if retime else {"pts": pts}
    # Preserve previously computed Bempp E fields if recomputing MFS reference, and update their errors vs new ref_E
    for k_old, v_old in old.items():
        if k_old.startswith(("bem_E", "bem_uniform_E")):
            out[k_old] = v_old
    for name, Js, kw in (("dense", range(3, 11), {}), ("one", range(3, 9), dict(q=0.5)),
                         ("notail", (4, 6, 8), dict(tails=False)), ("plain", (4, 6, 8, 10), None)):
        old_map = {int(r[0]): float(r[1]) for r in old.get(name, np.zeros((0, 6)))} if retime else {}
        rows = []
        for J in Js:
            lines._TREES.clear()
            gc.collect()
            E, t, info = (tetrahedron.solve(None, pts=pts, n_base=J) if kw is None
                          else tetrahedron.solve(J, pts=pts, **kw))
            best_t = min(t, old_map[J]) if J in old_map else t
            rows.append((J, best_t, info["unknowns"][0], info["rows"], info["optical"], info["ext"]))
            out[f"{name}_E{J}"] = E
            print(f"{name:6s} J={J:2d}: unk/cls={info['unknowns'][0]:4d}, "
                  f"t={t:5.2f}s (best={best_t:5.2f}s), opt={info['optical']:.2e}, ext={info['ext']:.6f}", flush=True)
        out[name] = np.array(rows)
    out["ref_E"] = out["dense_E10"]
    # Refresh stored Bempp error columns against new ref_E if Bempp fields are already present
    if "bem_uniform" in old and len(old["bem_uniform"]):
        out["bem_uniform"] = np.array([
            [r[0], r[1], r[2], core.rel(out[f"bem_uniform_E{int(r[0])}"], out["ref_E"])
             if f"bem_uniform_E{int(r[0])}" in out else r[3]]
            for r in old["bem_uniform"].tolist()
        ])
    if "--reset-bem" not in sys.argv and "bem" in old and len(old["bem"]):
        bem_old = old["bem"]
        if bem_old.shape[1] == 4:
            bem_old = np.column_stack((bem_old[:, 0], np.full(len(bem_old), 8.0), bem_old[:, 1:]))
        bem_upd = []
        for r in bem_old.tolist():
            ne_i, qr_i = int(r[0]), int(r[1])
            key_b = f"bem_E{ne_i}_q{qr_i}" if f"bem_E{ne_i}_q{qr_i}" in out else f"bem_E{ne_i}"
            eb = core.rel(out[key_b], out["ref_E"]) if key_b in out else float(r[4])
            bem_upd.append([ne_i, qr_i, int(r[2]), float(r[3]), eb])
        out["bem"] = np.array(sorted(bem_upd))
    else:
        out["bem"] = np.zeros((0, 5))
    np.savez(npz_path, **out)
    sys.exit(0)

# Load stored results
res = np.load(data("tetrahedron_convergence.npz"))
ref = res["ref_E"]
specs = {"dense": lambda J: dict(J=J), "one": lambda J: dict(J=J, q=0.5),
         "notail": lambda J: dict(J=J, tails=False), "plain": lambda J: dict(J=None, n_base=J)}
curves = {name: np.array([(sum(tetrahedron.Tetrahedron(classes=[], **specs[name](int(r[0]))).unknowns.values()),
                           rel(res[f"{name}_E{int(r[0])}"], ref), float(r[1]), int(r[0]))
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
bem_unif = res["bem_uniform"] if "bem_uniform" in res else np.zeros((0, 4))
if len(bem_unif):
    curves["uniform"] = np.column_stack((bem_unif[:, 1], bem_unif[:, 3]))

# Root-exponential fit to levels J = 3..8 (J = 9 approaches the J = 10 reference)
c = curves["dense"][:-1]
p = np.polyfit(np.sqrt(c[:, 0]), np.log(c[:, 1]), 1)
print(f"tetrahedron fit: error = {np.exp(p[1]):.3g} exp({p[0]:.4f} sqrt(N_u)), max factor "
      f"{np.exp(np.abs(np.log(c[:, 1]) - np.polyval(p, np.sqrt(c[:, 0]))).max()):.2f}")

# Four-panel layout: stacked left column ((a) 3D face kites & skeleton, (b) unfolded wall O v_0 v_1),
# plus full-height (c) error vs sqrt(N_u) and (d) error vs wall time [s]
fig = plt.figure(figsize=(7.6, 3.55))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.14, 1.10], wspace=0.28)
gs_left = gs[0].subgridspec(2, 1, height_ratios=[1.12, 0.96], hspace=0.22)
ax_3d = fig.add_subplot(gs_left[0])
ax_wall = fig.add_subplot(gs_left[1])
ax_conv = fig.add_subplot(gs[1])
ax_time = fig.add_subplot(gs[2])

# Panel (a): 3D view showing graded three-kite boundary collocation on front face + interior skeleton walls
az, el = np.radians(-38.0), np.radians(20.0)
R_az = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
R_el = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
P2 = (R_el @ R_az)[:2]
V = tetrahedron.VERTS @ P2.T
O = np.zeros(2)
edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
for i, j in edges:
    col = "#cfe0f5" if (i, j) == (0, 1) else "#eae9e4"
    al = 0.62 if (i, j) == (0, 1) else 0.25
    ax_3d.add_patch(Polygon([O, V[i], V[j]], closed=True, facecolor=col, edgecolor="none", alpha=al, zorder=1))
for i in range(4):
    ax_3d.plot([0, V[i, 0]], [0, V[i, 1]], ls=(0, (2.5, 1.8)), lw=0.65, color=MUTED, zorder=2)

# Show graded 3-kite panels & collocation points on front visible face (0, 1, 2)
F_vis = np.array([tetrahedron.VERTS[0], tetrahedron.VERTS[1], tetrahedron.VERTS[2]])
c_vis = F_vis.mean(axis=0)
brk_vis = np.array([0.0, 0.06, 0.14, 0.32, 1.0])
for k_i in range(3):
    Ak = F_vis[k_i]
    M1 = 0.5 * (F_vis[k_i] + F_vis[(k_i + 1) % 3])
    M2 = 0.5 * (F_vis[k_i] + F_vis[(k_i + 2) % 3])
    # Kite boundary to centroid
    seg_c = np.vstack((M1, c_vis, M2)) @ P2.T
    ax_3d.plot(seg_c[:, 0], seg_c[:, 1], color="#8e8d88", lw=0.6, zorder=3)
    for b_v in brk_vis[1:-1]:
        s1 = np.array([(1 - b_v) * (1 - v) * Ak + (1 - b_v) * v * M1 + b_v * (1 - v) * M2 + b_v * v * c_vis
                       for v in np.linspace(0, 1, 10)]) @ P2.T
        s2 = np.array([(1 - u) * (1 - b_v) * Ak + (1 - u) * b_v * M1 + u * (1 - b_v) * M2 + u * b_v * c_vis
                       for u in np.linspace(0, 1, 10)]) @ P2.T
        ax_3d.plot(s1[:, 0], s1[:, 1], color="#b8b7b2", lw=0.38, zorder=3)
        ax_3d.plot(s2[:, 0], s2[:, 1], color="#b8b7b2", lw=0.38, zorder=3)
    # Graded collocation nodes (representative level)
    xi_n = np.array([0.025, 0.10, 0.23, 0.50, 0.82])
    Un, Vn = [g.ravel() for g in np.meshgrid(xi_n, xi_n)]
    p_kite = (((1 - Un) * (1 - Vn))[:, None] * Ak + ((1 - Un) * Vn)[:, None] * M1
              + (Un * (1 - Vn))[:, None] * M2 + (Un * Vn)[:, None] * c_vis) @ P2.T
    ax_3d.plot(p_kite[:, 0], p_kite[:, 1], ".", color=INK, ms=1.65, zorder=4)

for i, j in edges:
    ls_e = (0, (3, 2)) if 3 in (i, j) else "-"
    lw_e = 0.72 if 3 in (i, j) else 1.05
    ax_3d.plot([V[i, 0], V[j, 0]], [V[i, 1], V[j, 1]], ls=ls_e, lw=lw_e, color=INK, zorder=5)
ax_3d.plot(0, 0, "o", color=MUTED, ms=2.5, zorder=5)
ax_3d.text(-0.04, 0.03, r"$\mathcal{M}$", fontsize=7.5, color=BLUE, ha="right", va="bottom")
ax_3d.text(0.28, 0.28, "3-kite\ngraded $\\Gamma$", fontsize=6.8, color=INK, ha="left", va="center")
ax_3d.set(xlim=(-0.54, 0.62), ylim=(-0.47, 0.57), aspect="equal")
ax_3d.axis("off")
ax_3d.set_title(r"(a) 3D kites on $\Gamma$ & skeleton $\mathcal{M}$", loc="left", fontsize=8.1, pad=2)

# Panel (b): Unfolded 2D view of one triangular skeleton wall O v_0 v_1 (|t| <= 1/2, depth b in [0, A_TET])
h_wall = tetrahedron.A_TET
ax_wall.add_patch(Polygon([[-0.5, 0.0], [0.5, 0.0], [0.0, h_wall]], closed=True, facecolor="#f4f3ef", edgecolor="none", zorder=0))
# Inner 0.56-scaled tetrahedron intersection at depth (1 - 0.56) * h_wall = 0.156
b_in = (1.0 - 0.56) * h_wall
ax_wall.add_patch(Polygon([[-0.28, b_in], [0.28, b_in], [0.0, h_wall]], closed=True, facecolor="#e3ece8", edgecolor=GREEN,
                          lw=0.65, ls="--", zorder=1))
ax_wall.text(-0.028, 0.245, "inner base", fontsize=6.5, color=GREEN, ha="center", va="center")

# Free lightning edge lines d_j = 0.14 * (2^{-1/2})^j with stretched t(u) nodes clustering towards corners
d_vis = 0.14 * (2.0 ** -0.5) ** np.arange(0, 8)
u_samp = np.linspace(-0.88, 0.88, 9)
for idx_d, dv in enumerate(d_vis):
    tmax = 0.5 - 1.5 * dv
    ax_wall.plot([-tmax, tmax], [dv, dv], color=BLUE, lw=0.80, zorder=3)
    lam_v = np.log((1.0 - 1.5 * dv) / (1.5 * dv))
    t_nodes = 0.5 * np.tanh(0.5 * lam_v * u_samp)
    t_nodes = t_nodes[np.abs(t_nodes) <= tmax]
    ax_wall.plot(t_nodes, np.full_like(t_nodes, dv), "o", color=BLUE, ms=1.45, zorder=4)
    # Vertex dipoles along side rays O v_0 and O v_1
    frac = dv / h_wall
    for sgn_t in (-1.0, 1.0):
        ax_wall.plot(sgn_t * 0.5 * (1.0 - frac), dv, "s", color=NAVY, ms=1.7, zorder=4)

# 9-line singular edge tail below d_min + radial corner-tail rays at both vertices
ax_wall.plot([-0.485, 0.485], [0.006, 0.006], color=ORANGE, lw=1.1, ls=(0, (2, 1)), zorder=4)
for sgn_t in (-1.0, 1.0):
    for ang_f in (0.8, 0.4, 0.2):
        # Direction in (t, b) wall coordinates from corner (sgn_t * 0.5, 0)
        dt_r = -sgn_t * 1.0
        db_r = ang_f
        nrm = np.hypot(dt_r, db_r)
        ax_wall.plot([sgn_t * 0.5 + 0.012 * dt_r / nrm, sgn_t * 0.5 + 0.095 * dt_r / nrm],
                     [0.012 * db_r / nrm, 0.095 * db_r / nrm], color=ORANGE, lw=0.85, zorder=5)

# Dyadic boundary kite panel ticks along bottom edge v_1 -- v_0
for b_t in (0.0, 0.03, 0.07, 0.16, 0.50):
    for sgn_t in (-1.0, 1.0):
        tx = sgn_t * (0.5 - 0.5 * b_t)
        ax_wall.plot([tx, tx], [-0.014, 0.0], color=MUTED, lw=0.5, zorder=3)

# Zoom inset in upper-right of ax_wall magnifying corner v_0 = (0.5, 0) (3 corner rays + tail + line cut-off)
zw_t, zh_t = 0.095, 0.068
ax_wall.plot([0.5 - zw_t, 0.5, 0.5, 0.5 - zw_t, 0.5 - zw_t], [0, 0, zh_t, zh_t, 0], color=RUST, lw=0.55, zorder=6)
ins_w = ax_wall.inset_axes([0.58, 0.48, 0.42, 0.45])
ins_w.set_facecolor("white")
ins_w.add_patch(Polygon([[0.5 - zw_t, 0.0], [0.5, 0.0], [0.5 - 0.5 * zh_t / h_wall, zh_t], [0.5 - zw_t, zh_t]],
                        closed=True, facecolor="#f5f4f0", edgecolor="none"))
for dv in d_vis[d_vis <= zh_t]:
    tm = 0.5 - 1.5 * dv
    if tm > 0.5 - zw_t:
        ins_w.plot([0.5 - zw_t, tm], [dv, dv], color=BLUE, lw=0.8)
        ins_w.plot(tm, dv, "o", color=BLUE, ms=1.9)
    ins_w.plot(0.5 * (1.0 - dv / h_wall), dv, "s", color=NAVY, ms=2.0)
for dt_tail in 0.14 * (2.0 ** -0.5) ** 7 * 0.5 ** np.arange(1, 5):
    ins_w.plot([0.5 - zw_t, 0.5 - 1.5 * d_vis[-1]], [dt_tail, dt_tail], color=ORANGE, lw=0.65, ls=":")
for ang_f in (0.8, 0.4, 0.2):
    nrm = np.hypot(1.0, ang_f)
    r_pts = 0.072 * (0.68 ** np.arange(5))
    ins_w.plot(0.5 - r_pts / nrm, r_pts * ang_f / nrm, "o-", color=ORANGE, lw=0.65, ms=1.8, mfc="white", mew=0.55)
ins_w.plot([0.5 - zw_t, 0.5], [0.0, 0.0], color=INK, lw=1.0)
ins_w.plot([0.5, 0.5 - 0.5 * zh_t / h_wall], [0.0, zh_t], color=MUTED, lw=0.65, ls=(0, (2, 1.5)))
ins_w.text(0.501, zh_t - 0.005, "corner\nzoom", fontsize=5.9, color=RUST, ha="right", va="top")
ins_w.set(xlim=(0.5 - zw_t, 0.506), ylim=(-0.006, zh_t), xticks=[], yticks=[], aspect="equal")
for sp in ins_w.spines.values():
    sp.set_edgecolor(RUST)
    sp.set_linewidth(0.55)
ax_wall.add_artist(ConnectionPatch(xyA=(0.5 - zw_t, zh_t), coordsA=ax_wall.transData,
                                   xyB=(0.5 - zw_t, -0.006), coordsB=ins_w.transData,
                                   color=RUST, lw=0.5, ls=":", zorder=6))

ax_wall.plot([-0.5, 0.5], [0.0, 0.0], color=INK, lw=1.15, zorder=5)
ax_wall.plot([-0.5, 0.0, 0.5], [0.0, h_wall, 0.0], color=MUTED, lw=0.75, ls=(0, (3, 2)), zorder=2)
ax_wall.plot(0.0, h_wall, "o", color=MUTED, ms=2.4, zorder=5)
ax_wall.text(0.0, h_wall + 0.018, "$O$", fontsize=6.8, color=MUTED, ha="center", va="bottom")
ax_wall.text(-0.50, -0.025, "$v_1$", fontsize=6.8, color=INK, ha="center", va="top")
ax_wall.text(0.50, -0.025, "$v_0$", fontsize=6.8, color=INK, ha="center", va="top")
ax_wall.text(0.0, -0.028, r"graded $\Gamma$ + 9-line tail", fontsize=6.7, color=ORANGE, ha="center", va="top")
ax_wall.text(-0.42, 0.085, "lines\n$t(u)$", fontsize=6.5, color=BLUE, ha="center", va="bottom")
ax_wall.text(-0.30, 0.255, "vertex\ndipoles", fontsize=6.3, color=NAVY, ha="center", va="bottom")
ax_wall.set(xlim=(-0.56, 0.59), ylim=(-0.095, 0.45), aspect="equal")
ax_wall.axis("off")
ax_wall.set_title(r"(b) unfolded wall $O v_0 v_1$", loc="left", fontsize=8.1, pad=1)

# Panel (c): convergence against sqrt(N_u)
xs = np.array([np.sqrt(curves["dense"][0, 0]) - 4, np.sqrt(curves["dense"][-1, 0]) + 6])
ax_conv.semilogy(xs, np.exp(np.polyval(p, xs)), color=MUTED, lw=0.85, zorder=1)
ax_conv.annotate(rf"${np.exp(p[1]):.2f}\,\exp({p[0]:.4f}\sqrt{{N_u}})$",
                 xy=(108.0, np.exp(np.polyval(p, 108.0))), xytext=(88.0, 5.5e-3),
                 fontsize=7.0, color=MUTED, ha="left", va="bottom",
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.55, shrinkA=1, shrinkB=2))
for key, (col, mk, ls, label) in {
    "plain": (PURPLE, "^", ":", "plain MFS (inner copy)"),
    "notail": (GREEN, "s", "--", "MFS without tails"),
    "one": (NAVY, "h", (0, (6, 1.5)), "MFS (1 depth/octave)"),
    "dense": (BLUE, "o", "-", "MFS (2 depths/octave)"),
    **({"uniform": (ORANGE, "D", "-.", "Bempp (quasi-uniform)")} if "uniform" in curves else {}),
    "graded": (RUST, "X", (0, (3, 1)), "Bempp (graded, symmetry)"),
}.items():
    v = curves[key]
    curve(ax_conv, np.sqrt(v[:, 0]), v[:, 1], col, mk, ls, label, ms=4.2, lw=1.35, logx=False, zorder=3)
d_c = curves["dense"]
b_c = curves["graded"]
ax_conv.text(np.sqrt(d_c[0, 0]) + 2.5, d_c[0, 1] * 1.12, f"$J={int(d_c[0, 3])}$", fontsize=6.8, color=BLUE)
ax_conv.text(np.sqrt(d_c[-1, 0]) - 2.0, d_c[-1, 1] * 0.58, f"$J={int(d_c[-1, 3])}$", fontsize=6.8, color=BLUE, ha="right")
ax_conv.text(np.sqrt(b_c[-1, 0]) + 1.5, b_c[-1, 1] * 1.18, "8 GB\nlimit", fontsize=6.2, color=RUST, va="bottom", linespacing=0.95)
ax_conv.set(xlabel=r"$\sqrt{N_u}$ (whole-body unknowns)", ylabel="relative error of the scattered field",
            xlim=(14, 152), ylim=(4e-5, 1.0))
ax_conv.set_title(r"(c) error against $\sqrt{N_u}$", loc="left", fontsize=8.4)

# Panel (d): error against wall time [s]
m_tet_env = envelope([(float(r[2]), float(r[1])) for r in np.vstack((curves["one"], curves["dense"]))])
ax_time.loglog(m_tet_env[:, 0], m_tet_env[:, 1], "-", color=BLUE, lw=3.6, alpha=0.18, zorder=2)
if len(bem_unif):
    curve(ax_time, bem_unif[:, 2], bem_unif[:, 3], ORANGE, "D", "-.", ms=4.1, lw=1.3, zorder=3)
curve(ax_time, curves["one"][:, 2], curves["one"][:, 1], NAVY, "h", (0, (6, 1.5)), ms=4.2, lw=1.35, zorder=3)
curve(ax_time, curves["dense"][:, 2], curves["dense"][:, 1], BLUE, "o", "-", ms=4.2, lw=1.35, zorder=4)
curve(ax_time, bem_time_env[:, 0], bem_time_env[:, 1], RUST, "X", (0, (3, 1)), ms=4.2, lw=1.35, zorder=3)
ax_time.text(d_c[0, 2] * 1.15, d_c[0, 1] * 1.10, f"$J={int(d_c[0, 3])}$", fontsize=6.8, color=BLUE)
ax_time.text(d_c[-1, 2] * 0.88, d_c[-1, 1] * 0.58, f"$J={int(d_c[-1, 3])}$", fontsize=6.8, color=BLUE, ha="right")
ax_time.text(bem_time_env[-1, 0] * 1.02, bem_time_env[-1, 1] * 1.25, "8 GB\nlimit", fontsize=6.2, color=RUST, ha="right", va="bottom", linespacing=0.95)
ax_time.set(xlabel="wall time [s]  (2 cores)", xlim=(1.6, 390), ylim=(4e-5, 1.0))
ax_time.set_xticks([2, 5, 10, 20, 50, 100, 250])
ax_time.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:g}"))
ax_time.xaxis.set_minor_formatter(plt.NullFormatter())
ax_time.set_title("(d) error against wall time", loc="left", fontsize=8.4)

handles, labels = ax_conv.get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.54, -0.02),
           ncol=3, columnspacing=1.05, handletextpad=0.4)
for a in (ax_conv, ax_time):
    a.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.25)
finish(fig, None, "fig_tetrahedron", tight="tight")
