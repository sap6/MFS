"""
Figures 3 and 9: the source and boundary collocation layouts of the cube and
of the NASA almond, with zoom callouts and log-scale hierarchies making the
lightning / geometric clustering visually explicit.

    python figures/layouts.py
"""
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle
from style import BLUE, BODY, GREEN, INK, MUTED, ORANGE, RUST, finish, plt

from mfs import almond, cube


def _cube_panel_breaks(levels, hmax=cube.H / 4):
    H = cube.H
    dist = np.concatenate(([0.0], H * 0.5 ** np.arange(levels, -1, -1)))
    brk = [0.0]
    for a, b in zip(dist[:-1], dist[1:]):
        n = int(np.ceil((b - a) / hmax))
        brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
    return H - np.array(brk)[::-1]


def cube_layout(J=6, fname="fig_cube_layout"):
    """Figure 3: the source and boundary collocation layout of the cube at level J"""
    H, layout = cube.H, cube.Cube(J, classes=[])
    depths = layout.depths
    tail = depths[-1] * 0.5 ** np.arange(1, 10)
    ends = _cube_panel_breaks(max(2, int(np.ceil(np.log2(H / depths[-1])))), H / 4)
    nodes = np.unique(layout.x[layout.x[:, 0] < H, 0])
    c = np.sqrt(0.5)

    fig, axs = plt.subplots(1, 3, figsize=(7.5, 3.2), gridspec_kw=dict(width_ratios=[1.08, 1.15, 1.02]))

    # (a) Cross-section through the middle of an edge + corner zoom inset resolving lightning & tail
    ax = axs[0]
    ax.fill([0, H, H, 0], [0, 0, H, H], color=BODY, lw=0, zorder=0)
    ax.plot([0, H, H], [H, H, 0], color=INK, lw=1.25, zorder=2)
    ax.plot([0, H], [0, H], color=MUTED, lw=0.75, ls=(0, (4, 3)), zorder=1)

    # Panel dividers + Gauss boundary collocation points on both faces
    for e in ends:
        ax.plot([e, e], [H, H + 0.017], color=MUTED, lw=0.55, zorder=2)
        ax.plot([H, H + 0.017], [e, e], color=MUTED, lw=0.55, zorder=2)
    ax.plot(nodes, np.full_like(nodes, H + 0.006), ".", color=INK, ms=2.2, zorder=3, label="boundary nodes (4/panel)")
    ax.plot(np.full_like(nodes, H + 0.006), nodes, ".", color=INK, ms=2.2, zorder=3)

    # Free lines, tail lines, and inner-cube base dipoles
    ax.plot(H - depths * c, H - depths * c, "o", color=BLUE, ms=3.8, mec="white", mew=0.45, zorder=4, label=r"free lines ($\times 2^{-1/2}$)")
    ax.plot(H - tail * c, H - tail * c, "o", mfc="white", mec=ORANGE, ms=3.0, mew=0.75, zorder=5, label=r"9-line tail ($\times 1/2$)")
    s = (np.arange(8) + 0.5) / 8 * 2 * (H - 0.2) - (H - 0.2)
    s = s[s > 0]
    ax.plot(s, np.full_like(s, H - 0.2), "s", color=GREEN, ms=3.5, mec="white", mew=0.4, zorder=3, label="base dipoles")
    ax.plot(np.full_like(s, H - 0.2), s, "s", color=GREEN, ms=3.5, mec="white", mew=0.4, zorder=3)

    # Zoom box around (0.5, 0.5) + x16 inset placed in upper-left conductor region away from diagonal & base grid
    zw = 0.024
    ax.add_patch(Rectangle((H - zw, H - zw), zw + 0.004, zw + 0.004, fill=False, edgecolor=RUST, lw=0.75, zorder=6))
    ins = ax.inset_axes([0.03, 0.46, 0.38, 0.38])
    ins.set_facecolor("white")
    ins.fill([H - zw, H, H, H - zw], [H - zw, H - zw, H, H], color="#f5f4f0", lw=0)
    ins.plot([H - zw, H, H], [H, H, H - zw], color=INK, lw=1.1, zorder=2)
    ins.plot([H - zw, H], [H - zw, H], color=MUTED, lw=0.65, ls=(0, (3, 2)), zorder=1)
    for e in ends[ends >= H - zw]:
        ins.plot([e, e], [H, H + 0.0022], color=MUTED, lw=0.55)
        ins.plot([H, H + 0.0022], [e, e], color=MUTED, lw=0.55)
    nd_z = nodes[nodes >= H - zw]
    ins.plot(nd_z, np.full_like(nd_z, H + 0.0008), ".", color=INK, ms=2.4, zorder=3)
    ins.plot(np.full_like(nd_z, H + 0.0008), nd_z, ".", color=INK, ms=2.4, zorder=3)
    d_z = depths[depths * c <= zw]
    t_z = tail[tail * c <= zw]
    ins.plot(H - d_z * c, H - d_z * c, "o", color=BLUE, ms=3.6, mec="white", mew=0.45, zorder=4)
    ins.plot(H - t_z * c, H - t_z * c, "o", mfc="white", mec=ORANGE, ms=3.0, mew=0.75, zorder=5)
    ins.set(xlim=(H - zw, H + 0.0032), ylim=(H - zw, H + 0.0032), xticks=[], yticks=[], aspect="equal")
    for sp in ins.spines.values():
        sp.set_edgecolor(RUST)
        sp.set_linewidth(0.7)
    ins.text(H - 0.06 * zw, H - 0.90 * zw, r"$\times 16$ zoom", fontsize=6.4, color=RUST, ha="right", va="bottom")
    ax.add_patch(ConnectionPatch(xyA=(H - zw, H - 0.25 * zw), coordsA=ax.transData,
                                 xyB=(1.0, 0.72), coordsB=ins.transAxes, color=RUST, lw=0.55, ls=":"))

    ax.text(0.055, 0.088, r"skeleton $\mathcal{M}$", fontsize=7.2, color=MUTED, rotation=45, rotation_mode="anchor")
    ax.text(0.325, 0.07, "conductor", fontsize=7.3, color=MUTED)
    ax.text(0.503, 0.520, "edge", fontsize=7.5, color=INK)
    ax.set(xlim=(0, 0.565), ylim=(0, 0.565), aspect="equal", xlabel="$x$", ylabel="$y$", xticks=[0, 0.25, 0.5], yticks=[0, 0.25, 0.5])
    ax.set_title("(a) section & corner zoom", loc="left", fontsize=8.4)
    ax.legend(frameon=False, fontsize=7.1, loc="upper center", bbox_to_anchor=(0.48, -0.30), ncol=2,
              handletextpad=0.25, columnspacing=0.55)

    # (b) Matched logarithmic scale: dyadic boundary panels (top) vs lightning + tail lines (bottom)
    ax = axs[1]
    pan_d = H - ends[::-1]
    nd_d = H - nodes[::-1]
    y_bnd, y_src = 0.24, -0.14
    for k_p, (a_p, b_p) in enumerate(zip(pan_d[:-1], pan_d[1:])):
        if a_p > 0:
            col_p = "#e3e2dd" if k_p % 2 == 0 else "#f2f1ed"
            ax.fill_between([a_p, b_p], y_bnd - 0.055, y_bnd + 0.055, color=col_p, edgecolor=MUTED, lw=0.45, zorder=1)
    for p_v in pan_d[pan_d > 0]:
        ax.plot([p_v, p_v], [y_bnd - 0.065, y_bnd + 0.065], color=MUTED, lw=0.55, zorder=2)
    ax.plot(nd_d, np.full_like(nd_d, y_bnd), ".", color=INK, ms=2.6, zorder=3)
    ax.text(6.5e-3, y_bnd + 0.10, r"boundary $\Gamma$: dyadic panels (4 nodes)", ha="center", fontsize=6.9, color=INK)

    ax.plot(depths, np.full_like(depths, y_src), "o", color=BLUE, ms=3.8, mec="white", mew=0.45, zorder=4)
    ax.plot(tail, np.full_like(tail, y_src), "o", mfc="white", mec=ORANGE, ms=3.4, mew=0.85, zorder=4)
    ax.axvline(depths[-1], color=MUTED, lw=0.7, ls=":", zorder=0)
    ax.plot([depths[0], depths[0]], [y_src + 0.05, y_src + 0.13], color=MUTED, lw=0.55)
    ax.text(depths[0], y_src + 0.15, "$d_0$", ha="center", fontsize=7.6, color=INK)
    ax.text(depths[-1] * 1.18, 0.05, r"$d_{\min}$", ha="left", va="center", fontsize=7.6, color=MUTED)
    ax.text(5.2e-2, y_src - 0.13, "free lines\n(2 / octave)", ha="center", va="top", fontsize=7.1, color=BLUE)
    ax.text(6.5e-5, y_src - 0.13, "singular tail\n(1 / oct, $r^e$)", ha="center", va="top", fontsize=7.1, color=ORANGE)

    ax.set(xscale="log", xlim=(2.5e-6, 0.65), ylim=(-0.55, 0.46), yticks=[], xlabel=r"distance / depth from edge (log scale)")
    ax.set_box_aspect(0.92)
    ax.spines["left"].set_visible(False)
    ax.set_title("(b) matched log-scale grading", loc="left", fontsize=8.4)

    # (c) Stretched variable t(u) along an edge line, with discrete sample projection showing vertex clustering
    ax = axs[2]
    u = np.linspace(-1, 1, 400)
    for d, ls, col, lab in ((0.2, "-", GREEN, "$d = 0.2$"),
                            (0.2 / 16, "--", BLUE, "$d = 0.2/16$"),
                            (0.2 / 256, "-.", ORANGE, "$d = 0.2/256$")):
        lam = np.log((2 * H - d) / d)
        ax.plot(u, H * np.tanh(0.5 * lam * u), ls, color=col, lw=1.35, label=lab, zorder=2)

    d_demo = 0.2 / 16
    lam_d = np.log((2 * H - d_demo) / d_demo)
    u_samp = np.linspace(-0.92, 0.92, 15)
    t_samp = H * np.tanh(0.5 * lam_d * u_samp)
    ax.plot(u_samp, t_samp, "o", color=BLUE, ms=2.8, mec="white", mew=0.35, zorder=4)
    ax.plot(np.full_like(t_samp, 1.07), t_samp, "_", color=BLUE, ms=4.2, mew=0.75, clip_on=False, zorder=4)
    ax.text(1.07, 0.54, "$t$-nodes", fontsize=6.6, color=BLUE, ha="center", clip_on=False)

    ax.axhline(H, color=MUTED, lw=0.55, ls=":")
    ax.axhline(-H, color=MUTED, lw=0.55, ls=":")
    ax.text(-0.96, 0.40, "vertex $+1/2$", fontsize=6.7, color=MUTED)
    ax.text(0.18, -0.45, "vertex $-1/2$", fontsize=6.7, color=MUTED)
    ax.set(xlabel="stretched variable $u$", ylabel="edge coordinate $t$", xlim=(-1, 1.15), ylim=(-0.55, 0.55),
           xticks=[-1, 0, 1])
    ax.set_box_aspect(0.92)
    ax.legend(frameon=False, fontsize=7.1, loc="upper center", bbox_to_anchor=(0.46, -0.30), ncol=2,
              handlelength=1.5, handletextpad=0.3, columnspacing=0.6)
    ax.set_title("(c) vertex clustering $t(u)$", loc="left", fontsize=8.4)
    for ax in axs:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.29, wspace=0.34)
    finish(fig, None, fname, grid=False, tight="tight")


def almond_layout(n_mid=6, n_F=3, fname="fig_almond_layout"):
    """Figure 9: mid-plane source rows, 3:1 cross-section boundary & rim grading, and conical tip tail"""
    L = almond.L
    etas = almond.row_etas(n_mid, n_F)
    y_all, (tpts, trow, trr) = almond.source_rows(n_mid, n_F)
    xb, wb, _ = almond.collocation(s=0.35, hmax=0.012)

    pts_mid, pts_foc = [], []
    x0_cut, xt_tip = L * (almond.T1 - 0.004), L * almond.T1
    for i, eta in enumerate(etas):
        target = pts_mid if i < n_mid else pts_foc
        x = 0.5 * 0.7 * almond.surf_dist(0.0, eta * almond.semi_axes(0.0)[0])
        while x < x0_cut:
            yv = eta * almond.semi_axes(x)[0]
            target.append((x, yv))
            x += 0.7 * almond.surf_dist(x, yv)
        Amu = np.sqrt(almond.A_ ** 2 - almond.B_ ** 2 * (1 - eta ** 2))
        psi = 0.5 * 0.7 * almond.surf_dist(0.0, eta * almond.B_) / Amu
        while psi < np.pi / 2:
            xr_i, yr_i = -Amu * np.sin(psi), eta * almond.B_ * np.cos(psi)
            ds = 0.7 * almond.surf_dist(xr_i, yr_i)
            if yr_i > 0.3 * ds:
                target.append((xr_i, yr_i))
            psi += ds / np.sqrt((Amu * np.cos(psi)) ** 2 + (eta * almond.B_ * np.sin(psi)) ** 2)
    pts_mid, pts_foc = np.array(pts_mid) * 1e3, np.array(pts_foc) * 1e3

    xf = np.linspace(0, L * almond.T1, 400)
    af = almond.semi_axes(xf)[0]
    psi = np.linspace(0, np.pi / 2, 200)
    xr, ar = -almond.A_ * np.sin(psi), almond.B_ * np.cos(psi)
    X = np.concatenate((xr[::-1], xf)) * 1e3
    Y = np.concatenate((ar[::-1], af)) * 1e3
    x_tip_mm = xt_tip * 1e3

    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.85), gridspec_kw=dict(width_ratios=[1.52, 1.12, 1.14]))

    # (a) Mid-plane z = 0: inner Chebyshev rows vs outer geometric focal rows
    ax = axs[0]
    ax.fill_between(X, 0, Y, color=BODY, lw=0, zorder=0)
    ax.plot(X, Y, color=INK, lw=1.2, zorder=3)
    ax.plot(X, almond.ETA_F * Y, color=ORANGE, lw=1.0, ls=(0, (4, 2)), zorder=4, label=r"focal curve ($\eta_F = 0.943$)")
    ax.plot(pts_mid[:, 0], pts_mid[:, 1], ".", color=GREEN, ms=2.3, zorder=2, label=f"{n_mid} Chebyshev inner rows")
    ax.plot(pts_foc[:, 0], pts_foc[:, 1], ".", color=BLUE, ms=2.4, zorder=3, label=f"{n_F + 1} geometric rim rows")
    ax.axvline(0.0, color=MUTED, lw=0.65, ls=":", zorder=1)
    ax.text(3.5, 53.0, "cut (b)", fontsize=6.9, color=MUTED)
    ax.plot([x_tip_mm - 5.5, x_tip_mm - 5.5, x_tip_mm + 0.4, x_tip_mm + 0.4, x_tip_mm - 5.5],
            [0, 3.9, 3.9, 0, 0], color=RUST, lw=0.65, zorder=5)
    ax.text(x_tip_mm - 2.0, 21.0, "tip (c)", fontsize=6.9, color=RUST, ha="right")
    ax.set(xlabel="$x$ [mm]", ylabel="$y$ [mm]", aspect="equal", xlim=(-112, 156), ylim=(0, 62))
    ax.set_anchor("N")
    ax.legend(frameon=False, fontsize=7.0, loc="upper center", bbox_to_anchor=(0.48, -0.82), ncol=2,
              markerscale=2.2, handletextpad=0.3, columnspacing=0.55)
    ax.set_title("(a) mid-plane $z = 0$ (layout $6{:}3$)", loc="left", fontsize=8.4)

    # (b) Cross-section at x = 0 (3:1 ellipse) + magnified rim bar resolving geometric rows & boundary grading
    ax = axs[1]
    a0, b0 = almond.B_ * 1e3, almond.C_ * 1e3
    ph_grid = np.linspace(0, np.pi / 2, 200)
    ax.fill_between(a0 * np.cos(ph_grid), 0, b0 * np.sin(ph_grid), color=BODY, lw=0, zorder=0)
    ax.plot(a0 * np.cos(ph_grid), b0 * np.sin(ph_grid), color=INK, lw=1.2, zorder=2)

    phis_col = np.unique(np.round(np.arctan2(xb[:, 2] / almond.C_, xb[:, 1] / almond.B_), 6))
    ax.plot(a0 * np.cos(phis_col), b0 * np.sin(phis_col), "o", color=INK, ms=2.6, zorder=4, label=r"boundary rings on $\Gamma$")
    etas_arr = np.array(etas)
    ax.plot(etas_arr[:n_mid] * a0, np.zeros(n_mid), "o", color=GREEN, ms=3.1, zorder=4)
    ax.plot(etas_arr[n_mid:] * a0, np.zeros(n_F + 1), "o", color=BLUE, ms=3.1, zorder=5)
    ax.plot(almond.ETA_F * a0, 0.0, "*", color=ORANGE, ms=5.5, zorder=6)

    ins_b = ax.inset_axes([0.10, 0.56, 0.86, 0.40])
    ins_b.set_facecolor("white")
    ph_z = np.linspace(0, np.arccos(0.81), 120)
    ins_b.fill_between(np.cos(ph_z), 0, np.sin(ph_z) / 3.0, color="#f5f4f0", lw=0)
    ins_b.plot(np.cos(ph_z), np.sin(ph_z) / 3.0, color=INK, lw=1.05)
    ph_cz = phis_col[np.cos(phis_col) >= 0.81]
    ins_b.plot(np.cos(ph_cz), np.sin(ph_cz) / 3.0, "o", color=INK, ms=2.6, zorder=4)
    ins_b.plot(etas_arr[n_mid:], np.zeros(n_F + 1), "o", color=BLUE, ms=3.6, mec="white", mew=0.4, zorder=5)
    ins_b.axvline(almond.ETA_F, color=ORANGE, lw=0.85, ls=(0, (3, 1.5)), zorder=3)
    ins_b.text(almond.ETA_F + 0.006, 0.148, r"$\eta_F$", fontsize=6.8, color=ORANGE, ha="left")
    ins_b.text(0.816, 0.048, r"rim zoom ($\times 1/2$)", fontsize=6.3, color=BLUE)
    ins_b.set(xlim=(0.81, 1.012), ylim=(-0.035, 0.21), xticks=[0.85, 0.94, 1.0], yticks=[])
    ins_b.tick_params(labelsize=6.3, pad=1.5)
    for sp in ("top", "right", "left"):
        ins_b.spines[sp].set_visible(False)

    ax.set(xlabel="$y$ [mm]", ylabel="$z$ [mm]", xlim=(0, 52), ylim=(-2, 42))
    ax.legend(frameon=False, fontsize=7.0, loc="upper center", bbox_to_anchor=(0.48, -0.27), ncol=1, handletextpad=0.3)
    ax.set_title("(b) section $x = 0$ & rim zoom", loc="left", fontsize=8.4)

    # (c) Enlargement of the conical tip + log-distance tail strip (1 mm down to 10^{-4} mm)
    ax = axs[2]
    sel = X > x_tip_mm - 5.6
    ax.fill_between(X[sel] - x_tip_mm, 0, Y[sel], color=BODY, lw=0, zorder=0)
    ax.plot(X[sel] - x_tip_mm, Y[sel], color=INK, lw=1.15, zorder=3)
    ax.plot(X[sel] - x_tip_mm, almond.ETA_F * Y[sel], color=ORANGE, lw=0.9, ls=(0, (3, 1.5)), zorder=4)

    xb_tip = xb[xb[:, 0] * 1e3 > x_tip_mm - 5.5]
    ux_tip = np.unique(np.round(xb_tip[:, 0], 9))
    rim_b = np.array([xb_tip[np.isclose(xb_tip[:, 0], ux)][np.argmin(xb_tip[np.isclose(xb_tip[:, 0], ux), 2])] for ux in ux_tip])
    ax.plot(rim_b[:, 0] * 1e3 - x_tip_mm, almond.semi_axes(rim_b[:, 0])[0] * 1e3 + 0.08,
            ".", color=INK, ms=2.8, zorder=5, label=r"$\Gamma$ rings ($\geq 0.5$ mm)")

    near_m = pts_mid[:, 0] > x_tip_mm - 5.6
    near_f = pts_foc[:, 0] > x_tip_mm - 5.6
    ax.plot(pts_mid[near_m, 0] - x_tip_mm, pts_mid[near_m, 1], ".", color=GREEN, ms=2.8, zorder=2)
    ax.plot(pts_foc[near_f, 0] - x_tip_mm, pts_foc[near_f, 1], ".", color=BLUE, ms=2.8, zorder=3, label=r"free ($\geq 1$ mm)")
    ax.plot(tpts[:, 0] * 1e3 - x_tip_mm, tpts[:, 1] * 1e3, "o", mfc="white", mec=ORANGE, ms=2.2, mew=0.55,
            zorder=4, label=r"tail ($r^{\nu+2}$)")
    xcut = -0.004 * L * 1e3
    ax.plot([xcut, xcut], [0.0, 2.35], color=MUTED, lw=0.65, ls=":", zorder=1)
    ax.text(xcut - 0.12, 1.85, "1 mm\ncut-off", ha="right", va="center", fontsize=6.7, color=MUTED)

    ins_c = ax.inset_axes([0.44, 0.67, 0.54, 0.25])
    ins_c.set_facecolor("white")
    r_bnd = x_tip_mm - rim_b[:, 0] * 1e3
    r_row_tail = (xt_tip - tpts[trow == n_mid + n_F, 0]) * 1e3
    r_row_free = x_tip_mm - pts_foc[(pts_foc[:, 0] > x_tip_mm - 5.5) & (pts_foc[:, 1] > 0.85 * almond.ETA_F * almond.semi_axes(pts_foc[:, 0] * 1e-3)[0] * 1e3), 0]
    ins_c.plot(r_bnd, np.full_like(r_bnd, 0.26), ".", color=INK, ms=2.6, zorder=4)
    ins_c.plot(r_row_free[::2], np.full_like(r_row_free[::2], -0.14), "o", color=BLUE, ms=2.5, mec="white", mew=0.3)
    ins_c.plot(r_row_tail[::11], np.full_like(r_row_tail[::11], -0.14), "o", mfc="white", mec=ORANGE, ms=2.5, mew=0.65)
    ins_c.axvline(abs(xcut), color=MUTED, lw=0.6, ls=":")
    ins_c.set(xscale="log", xlim=(8e-5, 6.2), ylim=(-0.48, 0.62), xticks=[1e-4, 1e-2, 1.0], yticks=[])
    ins_c.tick_params(labelsize=6.4, pad=1.2)
    ins_c.set_title(r"$\Gamma$ & outer row vs $r$ [mm]", fontsize=6.5, color=ORANGE, loc="left", pad=1.8)
    for sp in ("top", "right", "left"):
        ins_c.spines[sp].set_visible(False)

    ax.set(xlabel="offset from tip [mm]", ylabel="$y$ [mm]", xlim=(-5.5, 0.35), ylim=(0, 4.85))
    ax.legend(frameon=False, fontsize=6.9, loc="upper center", bbox_to_anchor=(0.48, -0.27), ncol=2,
              markerscale=1.5, handletextpad=0.25, columnspacing=0.5)
    ax.set_title("(c) conical tip & 4-decade tail", loc="left", fontsize=8.4)
    for ax in axs:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.89, bottom=0.30, wspace=0.30)
    finish(fig, None, fname, grid=False, tight="tight")


cube_layout()
almond_layout()
