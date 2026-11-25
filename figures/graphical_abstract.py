"""
Graphical abstract for the Journal of Computational Physics:
  (a) Source & graded boundary placement on the medial skeleton of a polyhedral edge, with corner zoom;
  (b) Root-exponential convergence on convex polyhedra (cube and regular tetrahedron) vs inner-surface stall;
  (c) Monostatic radar cross-section of the NASA almond at 10.25 GHz (8.6 wavelengths) vs reference.

    python figures/graphical_abstract.py
"""
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle
from style import BLUE, BODY, GREEN, INK, MUTED, NAVY, ORANGE, PURPLE, RUST, curve, data, finish, plt, rel

from mfs import cube, tetrahedron

fig, axs = plt.subplots(1, 3, figsize=(8.2, 2.95), gridspec_kw=dict(width_ratios=[1.02, 1.16, 1.16]))

# (a) Skeleton wall, graded boundary collocation, and x16 corner zoom resolving lightning + tail
ax = axs[0]
H = cube.H
layout = cube.Cube(5, classes=[])
depths = layout.depths
tail = depths[-1] * 0.5 ** np.arange(1, 7)
nodes = cube.graded_nodes(max(2, int(np.ceil(np.log2(H / depths[-1])))), 4, H / 4)[0]
ax.fill([0, H, H, 0], [0, 0, H, H], color=BODY, lw=0, zorder=0)
ax.plot([0, H, H], [H, H, 0], color=INK, lw=1.25, zorder=2)
ax.plot([0, H], [0, H], color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.plot(nodes, np.full_like(nodes, H + 0.006), ".", color=INK, ms=2.1, zorder=3, label=r"graded $\Gamma$ nodes")
ax.plot(np.full_like(nodes, H + 0.006), nodes, ".", color=INK, ms=2.1, zorder=3)

c = np.sqrt(0.5)
ax.plot(H - depths * c, H - depths * c, "o", color=BLUE, ms=3.9, mec="white", mew=0.45, zorder=4, label="free lines")
ax.plot(H - tail * c, H - tail * c, "o", mfc="white", mec=ORANGE, ms=3.2, mew=0.8, zorder=5, label=r"tail ($2^{-ie}$)")
s = np.array([0.075, 0.225])
ax.plot(s, np.full_like(s, H - 0.2), "s", color=GREEN, ms=3.5, mec="white", mew=0.4, zorder=3, label="base dipoles")
ax.plot(np.full_like(s, H - 0.2), s, "s", color=GREEN, ms=3.5, mec="white", mew=0.4, zorder=3)

zw = 0.032
ax.add_patch(Rectangle((H - zw, H - zw), zw + 0.005, zw + 0.005, fill=False, edgecolor=RUST, lw=0.7, zorder=6))
ins = ax.inset_axes([0.04, 0.45, 0.40, 0.40])
ins.set_facecolor("white")
ins.fill([H - zw, H, H, H - zw], [H - zw, H - zw, H, H], color="#f5f4f0", lw=0)
ins.plot([H - zw, H, H], [H, H, H - zw], color=INK, lw=1.05, zorder=2)
ins.plot([H - zw, H], [H - zw, H], color=MUTED, lw=0.65, ls=(0, (3, 2)), zorder=1)
nd_z = nodes[nodes >= H - zw]
ins.plot(nd_z, np.full_like(nd_z, H + 0.001), ".", color=INK, ms=2.3, zorder=3)
ins.plot(np.full_like(nd_z, H + 0.001), nd_z, ".", color=INK, ms=2.3, zorder=3)
d_z = depths[depths * c <= zw]
t_z = tail[tail * c <= zw]
ins.plot(H - d_z * c, H - d_z * c, "o", color=BLUE, ms=3.6, mec="white", mew=0.4, zorder=4)
ins.plot(H - t_z * c, H - t_z * c, "o", mfc="white", mec=ORANGE, ms=3.0, mew=0.75, zorder=5)
ins.set(xlim=(H - zw, H + 0.004), ylim=(H - zw, H + 0.004), xticks=[], yticks=[], aspect="equal")
for sp in ins.spines.values():
    sp.set_edgecolor(RUST)
    sp.set_linewidth(0.65)
ins.text(H - 0.06 * zw, H - 0.90 * zw, "zoom", fontsize=6.4, color=RUST, ha="right", va="bottom")
ax.add_patch(ConnectionPatch(xyA=(H - zw, H - 0.25 * zw), coordsA=ax.transData,
                             xyB=(1.0, 0.72), coordsB=ins.transAxes, color=RUST, lw=0.5, ls=":"))

ax.text(0.052, 0.088, r"skeleton $\mathcal{M}$", fontsize=7.4, color=MUTED, rotation=45, rotation_mode="anchor")
ax.text(0.31, 0.07, "conductor", fontsize=7.4, color=MUTED)
ax.text(0.503, 0.518, "edge", fontsize=7.6, color=INK)
ax.set(xlim=(0, 0.565), ylim=(0, 0.565), aspect="equal", xticks=[0, 0.25, 0.5], yticks=[0, 0.25, 0.5],
       xlabel="$x$", ylabel="$y$")
ax.set_title("1. Skeleton + singular tail", loc="left", fontsize=8.6, fontweight="bold")
ax.legend(frameon=False, fontsize=7.1, loc="upper center", bbox_to_anchor=(0.48, -0.28), ncol=2,
          handletextpad=0.25, columnspacing=0.6)
ax.grid(False)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (b) Polyhedral convergence: cube and regular tetrahedron at k = 2 pi
ax = axs[1]
rc = np.load(data("cube_convergence.npz"))
ref_c = rc["ref_E"]
c_plain = np.array([(4 * r[2], rel(rc[f"plain_E{int(r[0])}"], ref_c)) for r in rc["plain"] if 4 <= r[0] < 12])
c_notail = np.array([(4 * r[2], rel(rc[f"notail_E{int(r[0])}"], ref_c)) for r in rc["notail"] if 4 <= r[0] < 12])
c_dense = np.array([(4 * r[2], rel(rc[f"dense_E{int(r[0])}"], ref_c)) for r in rc["dense"] if 3 <= r[0] <= 11])
cu = np.load(data("cube_bempp_uniform.npz"))
c_unif = np.array([(n, rel(E, ref_c)) for n, E in zip(cu["ndof"], cu["E"]) if n > 1000])
best_b = {}
for *_, N, _, err in np.load(data("cube_cases_k1.npz"))["bem"]:
    best_b[int(N)] = min(best_b.get(int(N), 1.0), err)
c_bem = np.array(sorted(best_b.items()))

rt = np.load(data("tetrahedron_convergence.npz"))
ref_t = rt["ref_E"]
t_dense = np.array([(sum(tetrahedron.Tetrahedron(int(r[0]), classes=[]).unknowns.values()),
                     rel(rt[f"dense_E{int(r[0])}"], ref_t)) for r in rt["dense"] if 3 <= r[0] <= 9])

curve(ax, np.sqrt(c_plain[:, 0]), c_plain[:, 1], PURPLE, "^", ":", "inner surface", ms=4.0, lw=1.25, logx=False)
curve(ax, np.sqrt(c_notail[:, 0]), c_notail[:, 1], GREEN, "s", "--", "no tail", ms=3.9, lw=1.25, logx=False)
curve(ax, np.sqrt(c_unif[:, 0]), c_unif[:, 1], ORANGE, "D", "-.", "Bempp (uniform)", ms=3.9, lw=1.25, logx=False)
curve(ax, np.sqrt(c_bem[:, 0]), c_bem[:, 1], RUST, "X", (0, (3, 1)), "Bempp (graded)", ms=4.1, lw=1.3, logx=False)
curve(ax, np.sqrt(t_dense[:, 0]), t_dense[:, 1], NAVY, "h", (0, (5, 1.5)), "MFS, tetrahedron", ms=4.1, lw=1.35, logx=False)
curve(ax, np.sqrt(c_dense[:, 0]), c_dense[:, 1], BLUE, "o", "-", "MFS, cube", ms=4.2, lw=1.45, logx=False)
ax.text(np.sqrt(c_bem[-1, 0]), c_bem[-1, 1] * 1.48, "8 GB limit", fontsize=6.3, color=RUST, ha="center", va="bottom")
ax.set(xlabel=r"$\sqrt{N_u}$ (whole-body unknowns)", ylabel="relative field error", ylim=(5e-6, 0.8))
ax.set_title(r"2. Polyhedra ($k = 2\pi$): $\exp(-c\sqrt{N_u})$", loc="left", fontsize=8.6, fontweight="bold")
ax.legend(frameon=False, fontsize=6.8, loc="upper center", bbox_to_anchor=(0.48, -0.28), ncol=2,
          handlelength=1.55, handletextpad=0.32, columnspacing=0.55)

# (c) NASA almond monostatic RCS at 10.25 GHz (8.6 wavelengths)
ax = axs[2]
ra = np.load(data("almond_rcs_10.25GHz.npz"))
d_hh = np.loadtxt(data("austin_almond/ref_rcs.III.A.s1.fx4.H.txt"))
keep = np.isclose(np.round(d_hh[:, 2], 2) * 2, np.round(np.round(d_hh[:, 2], 2) * 2))
phi, ref_hh = np.round(d_hh[keep, 2], 2), d_hh[keep, 3]
d_vv = np.loadtxt(data("austin_almond/ref_rcs.III.A.s1.fx4.V.txt"))
keep_v = np.isclose(np.round(d_vv[:, 2], 2) * 2, np.round(np.round(d_vv[:, 2], 2) * 2))
ref_vv = d_vv[keep_v, 3]

ax.plot(phi, 10 * np.log10(ra["sigma_VV"]), "-", color=NAVY, lw=1.2, label="MFS VV (0.02 dB)", zorder=2)
ax.plot(phi, ref_vv, ls=(0, (2.2, 1.6)), color=RUST, lw=0.88, label="Austin ref. VV", zorder=3)
ax.plot(phi, 10 * np.log10(ra["sigma_HH"]), "-", color=BLUE, lw=1.25, label="MFS HH (0.12 dB)", zorder=4)
ax.plot(phi, ref_hh, ls=(0, (2.2, 1.6)), color=ORANGE, lw=0.90, label="Austin ref. HH", zorder=5)
ax.set(xlim=(0, 180), ylim=(-72, -16), xticks=[0, 60, 120, 180], xlabel="azimuth [deg]", ylabel="monostatic RCS [dBsm]")
ax.set_title("3. NASA almond at 10.25 GHz", loc="left", fontsize=8.6, fontweight="bold")
ax.legend(frameon=False, fontsize=7.1, loc="upper center", bbox_to_anchor=(0.48, -0.28), ncol=2,
          handlelength=1.6, handletextpad=0.35, columnspacing=0.65)

for ax in (axs[1], axs[2]):
    ax.grid(True, which="major", color="#e4e3df", lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.06, right=0.99, top=0.89, bottom=0.30, wspace=0.33)
finish(fig, None, "fig_graphical_abstract", tight="tight")
