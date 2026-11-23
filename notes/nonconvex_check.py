"""
3D non-convex L-shaped conducting polyhedron prototype (Deliverable 2):
    Omega_L = ([-0.5, 0.5]^2 \\ (0, 0.5]^2) x [-0.25, 0.25],   k = 2 * pi,
illuminated head-on into the 270-deg re-entrant concave notch by
    d = (-1/sqrt(2), -1/sqrt(2), 0),   p = (0, 0, 1),
exploiting exact z -> -z reflection symmetry (class c = (1, 1, -1)) identically in
both MFS (mfs.lines + mfs.core) and Bempp-cl (mfs.bem.Mirror).

Stages (each < 600s under mfs-safe-run):
    mfs-safe-run python notes/nonconvex_check.py --part cone
    mfs-safe-run python notes/nonconvex_check.py --part sanity
    mfs-safe-run python notes/nonconvex_check.py --part ref
    mfs-safe-run python notes/nonconvex_check.py --part mfs
    mfs-safe-run python notes/nonconvex_check.py --part bem
    python notes/nonconvex_check.py
"""
import gc
import os
import resource
import sys
import time

import numpy as np
import scipy.linalg as sla
from scipy.linalg import blas
from scipy.special import lpmv

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

from mfs import cone, core, lines  # noqa: E402
from mfs.cube import NU as NU_CUBE, TAU as TAU_CUBE, TAIL_E, RULE  # noqa: E402

NPZ = os.path.join(REPO, "data", "nonconvex_check.npz")
K = 2.0 * np.pi
D_W = np.array([-1.0 / np.sqrt(2.0), -1.0 / np.sqrt(2.0), 0.0])
P_W = np.array([0.0, 0.0, 1.0])
REFL = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, -1.0]])
CLS = (1, 1, -1)
SGN = np.array([[1.0, -1.0]])
PTS = 2.0 * core.fibonacci(400)
V_SAD = np.array([0.0, 0.0, 0.25])

# 2D L-polygon vertices in CCW order: 5 salient 90-deg corners + (0, 0) re-entrant 270-deg corner
L_POLY = np.array([
    [-0.5, -0.5], [0.5, -0.5], [0.5, 0.0], [0.0, 0.0], [0.0, 0.5], [-0.5, 0.5],
])

# All 18 physical boundary edges of Omega_L (6 vertical, 6 top z=+0.25, 6 bottom z=-0.25)
PHYS_EDGES = []
for _i in range(len(L_POLY)):
    _p0, _p1 = L_POLY[_i], L_POLY[(_i + 1) % len(L_POLY)]
    PHYS_EDGES.append((np.array([_p0[0], _p0[1], -0.25]), np.array([_p0[0], _p0[1], 0.25])))
    for _z in (-0.25, 0.25):
        PHYS_EDGES.append((np.array([_p0[0], _p0[1], _z]), np.array([_p1[0], _p1[1], _z])))


def clearance_L(pts):
    """Exact signed interior clearance inside Omega_L (> 0 strictly inside, <= 0 on/outside)."""
    pts = np.atleast_2d(pts)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    d_box = np.minimum.reduce([x + 0.5, 0.5 - x, y + 0.5, 0.5 - y, z + 0.25, 0.25 - z])
    d_notch = np.hypot(np.minimum(x, 0.0), np.minimum(y, 0.0))
    return np.minimum(d_box, d_notch)


def dist_phys_edges(pts):
    """Exact Euclidean distance from pts (M, 3) to the nearest of the 18 physical edges of Omega_L."""
    dmin = np.full(len(pts), np.inf)
    for a, b in PHYS_EDGES:
        ab = b - a
        t = np.clip(((pts - a) @ ab) / float(ab @ ab), 0.0, 1.0)
        dmin = np.minimum(dmin, np.linalg.norm(pts - (a[None, :] + t[:, None] * ab[None, :]), axis=1))
    return dmin


# ---------------------------------------------------------------------------
# 1. Exact Non-Convex Saddle Corner Exponents at (0, 0, +0.25) (~30 lines)
# ---------------------------------------------------------------------------

def saddle_cone_exponents(n_pole=14):
    """Leading (nu_sad, tau_sad) of the 5-octant exterior cone at (0, 0, +0.25) via lpmv(0, nu, -x.s)."""
    s2 = 1.0 / np.sqrt(2.0)
    rays = [([1., 0., 0.], [0., -s2, -s2]), ([0., 1., 0.], [-s2, 0., -s2]), ([0., 0., -1.], [-s2, -s2, 0.])]
    poles, t0s, t1s = [], [], []
    for v, t in [(np.array(v), np.array(t)) for v, t in rays]:
        for dl in 0.65 * np.exp(-2.8 * (np.sqrt(n_pole) - np.sqrt(np.arange(1, n_pole + 1)))):
            p = np.cos(dl) * v + np.sin(dl) * t
            tp = cone.unit(t - (t @ p) * p)
            poles.append(p); t0s.append(tp); t1s.append(np.cross(p, tp))
    for ph in np.linspace(-0.2 * np.pi, -1.3 * np.pi, max(6, n_pole // 2)):
        p = cone.unit(np.array([np.cos(ph), np.sin(ph), -0.9]))
        tp = cone.unit(np.array([0., 0., 1.]) - p[2] * p)
        poles.append(p); t0s.append(tp); t1s.append(np.cross(p, tp))
    poles, t0s, t1s = np.array(poles), np.array(t0s), np.array(t1s)
    assert np.all((poles[:, 2] < -1e-6) & ((poles[:, 0] < -1e-6) | (poles[:, 1] < -1e-6)))

    g, w = np.polynomial.legendre.leggauss(4)
    ed = np.concatenate(([0.0], 0.5 * 0.6 ** np.arange(20, -1, -1)))
    u01 = (0.5 * (ed[:-1, None] + ed[1:, None]) + 0.5 * (ed[1:, None] - ed[:-1, None]) * g).ravel()
    w01 = (0.5 * (ed[1:, None] - ed[:-1, None]) * w).ravel()
    th, wb1 = 0.5 * np.pi * np.concatenate((u01, 1.0 - u01[::-1])), 0.5 * np.pi * np.concatenate((w01, w01[::-1]))
    # 3*3 - 2*2 = 5 boundary quarter-great-circle arcs of the 3 conducting octants (z<0 and (x<0 or y<0))
    arcs = [([1., 0., 0.], [0., -1., 0.], [0., 0., 1.]), ([0., -1., 0.], [-1., 0., 0.], [0., 0., 1.]),
            ([-1., 0., 0.], [0., 1., 0.], [0., 0., 1.]), ([0., 1., 0.], [0., 0., -1.], [1., 0., 0.]),
            ([0., 0., -1.], [1., 0., 0.], [0., 1., 0.])]
    xb = np.vstack([np.outer(np.cos(th), u0) + np.outer(np.sin(th), u1) for u0, u1, _ in arcs])
    nb, wb = np.vstack([np.tile(nf, (len(th), 1)) for _, _, nf in arcs]), np.tile(wb1, 5)
    rho = np.min(np.arccos(np.clip(xb @ np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]]).T, -1, 1)), axis=1)
    pts = cone.unit(np.random.default_rng(0).normal(size=(6000, 3)))
    xi = pts[(pts[:, 2] > 0.08) | ((pts[:, 0] > 0.08) & (pts[:, 1] > 0.08))][:600]

    cb, ci = np.clip(-(xb @ poles.T), -1 + 1e-15, 1 - 1e-15), np.clip(-(xi @ poles.T), -1 + 1e-15, 1 - 1e-15)
    P_fn = lambda nu, c: lpmv(0, nu, c)
    dP_fn = lambda nu, c: nu * (P_fn(nu - 1.0, c) - c * P_fn(nu, c)) / (1.0 - c * c)
    out = []
    for neum in (False, True):
        wt = np.sqrt(wb * (rho if neum else 1.0))
        def sig(nu):
            pi, dpi = P_fn(nu, ci), dP_fn(nu, ci)
            C = np.stack((pi, -dpi * (xi @ t0s.T), -dpi * (xi @ t1s.T)), axis=2).reshape(len(xi), -1) / np.sqrt(len(xi))
            p, dp = P_fn(nu, cb), dP_fn(nu, cb)
            if not neum:
                B = np.stack((p, -dp * (xb @ t0s.T), -dp * (xb @ t1s.T)), axis=2).reshape(len(xb), -1) * wt[:, None]
            else:
                ns, d2p = nb @ poles.T, (nu * (nu + 1.0) * p - 2.0 * cb * dp) / (cb * cb - 1.0)
                B = np.stack((-dp * ns, d2p * ns * (xb @ t0s.T) - dp * (nb @ t0s.T),
                              d2p * ns * (xb @ t1s.T) - dp * (nb @ t1s.T)), axis=2).reshape(len(xb), -1) * wt[:, None]
            return cone.smallest_sigma(B, C * (20.0 if neum else 1.0))
        exps = cone.minima(sig, 0.2, 1.35, 0.04, tol=0.15)
        out.append((float(exps[0]), float(sig(exps[0]))))
    return out[0], out[1], len(poles)


# ---------------------------------------------------------------------------
# 2. Boundary Collocation on z >= 0 (3 top squares + 8 vertical half-faces)
# ---------------------------------------------------------------------------

def _intervals_1d(lo, hi, grade_lo, grade_hi, levels, hmax=0.125):
    """1D panel endpoints on [lo, hi], dyadically graded only towards requested physical-edge ends."""
    L = hi - lo
    half = 0.5 * L if (grade_lo and grade_hi) else L

    def _one_end(span):
        dist = np.concatenate(([0.0], span * 0.5 ** np.arange(levels, -1, -1)))
        dist = dist[dist <= span + 1e-15]
        brk = [0.0]
        for a, b in zip(dist[:-1], dist[1:]):
            n = max(1, int(np.ceil((b - a) / hmax)))
            brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
        return np.array(brk)

    if grade_lo and grade_hi:
        d = _one_end(half)
        pts = np.concatenate((lo + d, hi - d[:-1][::-1]))
        d_end = np.minimum(0.5 * (pts[:-1] + pts[1:]) - lo, hi - 0.5 * (pts[:-1] + pts[1:]))
    elif grade_lo:
        d = _one_end(L)
        pts = lo + d
        d_end = 0.5 * (pts[:-1] + pts[1:]) - lo
    elif grade_hi:
        d = _one_end(L)
        pts = hi - d[::-1]
        d_end = hi - 0.5 * (pts[:-1] + pts[1:])
    else:
        n = max(1, int(np.ceil(L / hmax)))
        pts = np.linspace(lo, hi, n + 1)
        d_end = np.full(n, L)
    return pts[:-1], pts[1:], d_end


def _patch_quad(u_spec, v_spec, levels, ng=4, hmax=0.25):
    """Anisotropic Gauss-Legendre rule on [u0, u1] x [v0, v1] graded towards flagged ends."""
    ua, ub, du = _intervals_1d(*u_spec, levels, hmax)
    va, vb, dv = _intervals_1d(*v_spec, levels, hmax)
    ng_t = max(3, ng - 1)
    rules = {n: np.polynomial.legendre.leggauss(n) for n in {ng, ng_t}}
    us, vs, ws = [], [], []
    for p in range(len(ua)):
        for q in range(len(va)):
            ratio = du[p] / max(dv[q], 1e-30)
            if 0.65 <= ratio <= 1.55:
                nu, nv = ng, ng
            elif du[p] > dv[q]:
                nu, nv = ng_t, ng
            else:
                nu, nv = ng, ng_t
            gu, wu = rules[nu]
            gv, wv = rules[nv]
            up = 0.5 * (ua[p] + ub[p]) + 0.5 * (ub[p] - ua[p]) * gu
            vq = 0.5 * (va[q] + vb[q]) + 0.5 * (vb[q] - va[q]) * gv
            wp = 0.5 * (ub[p] - ua[p]) * wu
            wq = 0.5 * (vb[q] - va[q]) * wv
            U, V = np.meshgrid(up, vq, indexing="ij")
            us.append(U.ravel())
            vs.append(V.ravel())
            ws.append(np.outer(wp, wq).ravel())
    return np.concatenate(us), np.concatenate(vs), np.concatenate(ws)


def build_collocation(d_min, nu_sad, ng=4, hmax=0.25):
    """
    Collocation points x, tangent frames T, and row weights rw on z >= 0:
    3 top 0.5 x 0.5 squares (z = 0.25) + 8 vertical 0.5 x 0.25 half-faces,
    graded dyadically towards all true boundary edges and saddle corner (0,0,0.25),
    NOT across smooth coplanar seams (x=0, y<0 or y=0, x<0, or z=0 midplane).
    """
    levels = max(2, int(np.ceil(np.log2(0.20 / d_min))))
    xs, ns, ws = [], [], []

    def add_face(pts, normal, w):
        xs.append(pts)
        ns.append(np.tile(normal, (len(pts), 1)))
        ws.append(w)

    # Top z = +0.25: split smooth interior seams (x=0,y<0 and y=0,x<0) at +-0.125 so collocation grades
    # dyadically ONLY towards the 6 physical top edges and the saddle corner (0, 0, 0.25).
    top_specs = [
        ((-0.5, -0.125, True, False), (-0.5, -0.125, True, False)),  # B00 outer corner
        ((-0.5, -0.125, True, False), (-0.125, 0.0, False, False)),  # B00 x=-0.5 strip (y=0 smooth)
        ((-0.125, 0.0, False, False), (-0.5, -0.125, True, False)),  # B00 y=-0.5 strip (x=0 smooth)
        ((-0.125, 0.0, False, True), (-0.125, 0.0, False, True)),    # B00 saddle corner (0, 0, 0.25)
        ((-0.5, -0.125, True, False), (0.0, 0.5, False, True)),      # B01 outer part (y=0 smooth)
        ((-0.125, 0.0, False, True), (0.0, 0.5, True, True)),        # B01 notch edge x=0 & saddle (0,0,0.25)
        ((0.0, 0.5, False, True), (-0.5, -0.125, True, False)),      # B10 outer part (x=0 smooth)
        ((0.0, 0.5, True, True), (-0.125, 0.0, False, True)),        # B10 notch edge y=0 & saddle (0,0,0.25)
    ]
    for x_sp, y_sp in top_specs:
        u, v, w = _patch_quad(x_sp, y_sp, levels, ng, hmax)
        add_face(np.column_stack((u, v, np.full_like(u, 0.25))), [0.0, 0.0, 1.0], w)

    # 8 vertical 0.5 x 0.25 half-faces on z in [0, 0.25] (z=0 is smooth mirror seam -> grade z only at +0.25)
    z_sp = (0.0, 0.25, False, True)
    vert_specs = [
        # (fixed_axis, fixed_val, u_spec, outward_normal)
        (0, -0.5, (-0.5, 0.0, True, False), [-1.0, 0.0, 0.0]),   # x = -0.5, y in [-0.5, 0] (y=0 smooth)
        (0, -0.5, (0.0, 0.5, False, True), [-1.0, 0.0, 0.0]),    # x = -0.5, y in [0, 0.5]  (y=0 smooth)
        (1, -0.5, (-0.5, 0.0, True, False), [0.0, -1.0, 0.0]),   # y = -0.5, x in [-0.5, 0] (x=0 smooth)
        (1, -0.5, (0.0, 0.5, False, True), [0.0, -1.0, 0.0]),    # y = -0.5, x in [0, 0.5]  (x=0 smooth)
        (0, 0.5, (-0.5, 0.0, True, True), [1.0, 0.0, 0.0]),      # x = +0.5, y in [-0.5, 0]
        (1, 0.5, (-0.5, 0.0, True, True), [0.0, 1.0, 0.0]),      # y = +0.5, x in [-0.5, 0]
        (0, 0.0, (0.0, 0.5, True, True), [1.0, 0.0, 0.0]),       # notch face x = 0, y in [0, 0.5] (+x into notch!)
        (1, 0.0, (0.0, 0.5, True, True), [0.0, 1.0, 0.0]),       # notch face y = 0, x in [0, 0.5] (+y into notch!)
    ]
    for ax_f, val_f, u_sp, n_out in vert_specs:
        u, z, w = _patch_quad(u_sp, z_sp, levels, ng, hmax)
        pts = np.zeros((len(u), 3))
        pts[:, ax_f] = val_f
        pts[:, 1 - ax_f] = u
        pts[:, 2] = z
        add_face(pts, n_out, w)

    x = np.ascontiguousarray(np.vstack(xs))
    n_vec = np.vstack(ns)
    area = np.concatenate(ws)
    assert abs(float(area.sum()) - 1.75) < 1e-13, f"expected half-area 1.75, got {area.sum()}"
    # Verify outward normals: x + eps * n_vec is strictly outside Omega_L, x - eps * n_vec is strictly inside
    assert np.all(clearance_L(x + 1e-6 * n_vec) <= 0.0) and np.all(clearance_L(x - 1e-6 * n_vec) > 5e-7)
    T = np.ascontiguousarray(core.tangent_basis(n_vec))
    rho_edge = dist_phys_edges(x)
    r_sad = np.linalg.norm(x - V_SAD[None, :], axis=1)
    p_sad = max(0.0, 1.0 - 2.0 * float(nu_sad))
    rw = np.sqrt(area * rho_edge * (r_sad ** p_sad))
    return x, T, rw


# ---------------------------------------------------------------------------
# 3. Three MFS Source Families on Omega_L
# ---------------------------------------------------------------------------

def inner_eroded_points(n_div, margin=0.07):
    """
    Family 1 ('inner'): electric dipoles on the z > 0 faces of an inset L-block at
    exact margin 0.07 inside Omega_L: ([-0.43, 0.43]^2 \\ (-0.07, 0.43]^2) x [-0.18, 0.18].
    """
    lo, hi, c_in, zh = -0.5 + margin, 0.5 - margin, -margin, 0.25 - margin
    h = 0.5 / n_div
    pts = []

    def grid_1d(a, b):
        m = max(1, int(np.round((b - a) / h)))
        return a + (np.arange(m) + 0.5) / m * (b - a)

    # Top L-face at z = +zh (three rectangles covering [-0.43, 0.43]^2 \ (-0.07, 0.43]^2)
    for xa, xb, ya, yb in ((lo, c_in, lo, c_in), (lo, c_in, c_in, hi), (c_in, hi, lo, c_in)):
        gx, gy = grid_1d(xa, xb), grid_1d(ya, yb)
        U, V = [g.ravel() for g in np.meshgrid(gx, gy, indexing="ij")]
        pts.append(np.column_stack((U, V, np.full_like(U, zh))))

    # 6 vertical walls of the inset L-block on z in (0, zh)
    gz = grid_1d(0.0, zh)
    v_walls = [
        (0, lo, lo, hi),     # x = -0.43, y in [-0.43, 0.43]
        (1, lo, lo, hi),     # y = -0.43, x in [-0.43, 0.43]
        (0, hi, lo, c_in),   # x = +0.43, y in [-0.43, -0.07]
        (1, hi, lo, c_in),   # y = +0.43, x in [-0.43, -0.07]
        (0, c_in, c_in, hi), # inset notch wall x = -0.07, y in [-0.07, 0.43]
        (1, c_in, c_in, hi), # inset notch wall y = -0.07, x in [-0.07, 0.43]
    ]
    for ax_f, val_f, ua, ub in v_walls:
        gu = grid_1d(ua, ub)
        U, Z = [g.ravel() for g in np.meshgrid(gu, gz, indexing="ij")]
        p = np.zeros((len(U), 3))
        p[:, ax_f], p[:, 1 - ax_f], p[:, 2] = val_f, U, Z
        pts.append(p)

    pe = np.ascontiguousarray(np.vstack(pts))
    cl = clearance_L(pe)
    assert np.all(pe[:, 2] > 1e-6) and np.all(cl >= margin - 1e-12), f"min clearance {cl.min()}"
    return pe


class LShapeMFS:
    """
    MFS solver on Omega_L in the z -> -z reflection sector CLS = (1, 1, -1):
      family == 'inner'          : eroded inner L-surface at margin 0.07 (parameter J = n_div)
      family == 'notail'         : medial sheets + 11 salient lines + 1 re-entrant diagonal sheet + vertex rays, no tails
      family == 'full'           : 'notail' + 9-line Meixner edge tails + 5 salient corner tails + 10-ray saddle tail
      family == 'no_saddle_tail' : 'full' without the 10-ray saddle-corner tail (ablation at J = 5)
    """
    def __init__(self, family, J, nu_sad, tau_sad, k=K, d0=0.08, c_end=0.72):
        self.family, self.J, self.k = family, J, float(k)
        self.nu_sad, self.tau_sad = float(nu_sad), float(tau_sad)
        s2 = 1.0 / np.sqrt(2.0)

        if family == "inner":
            self.depths = np.zeros(0)
            self.pe = inner_eroded_points(int(J), margin=0.07)
            self.pm = np.zeros((0, 3))
            self.lines = []
            self.v_tails = []
            d_coll = 0.08 * 2.0 ** (-4)   # well-resolved J=4 graded boundary norm for all inner N_u
        else:
            assert family in ("notail", "full", "no_saddle_tail")
            tails = family in ("full", "no_saddle_tail")
            sad_tail = (family == "full")
            self.depths = d0 * 2.0 ** (-0.5 * np.arange(2 * J + 1))
            d_coll = float(self.depths[-1])
            n0 = max(4, 2 * J - 4)

            # Base dipoles on inset inner-L core (margin 0.08, meeting d_0=0.08) + medial spine x=-0.25, y=-0.25
            b_surf = inner_eroded_points(6, margin=d0)
            med_pts = []
            for u_m in np.linspace(-0.33, 0.33, 6):
                for z_m in (0.045, 0.11):
                    med_pts.append([-0.25, u_m, z_m])
                    if abs(u_m + 0.25) > 0.04:
                        med_pts.append([u_m, -0.25, z_m])
            base = np.vstack((b_surf, np.array(med_pts)))
            assert np.all(clearance_L(base) >= d0 - 1e-12)

            # Vertex rays (both electric & magnetic dipoles) at 5 salient corners + 3 conducting octants at V_SAD
            # Offset by (d_j / sqrt(2)) along each face normal => clearance == d_j / sqrt(2) = 0.7071 * d_j > 0.7 * d_j!
            salient_corners = [
                ([-0.5, -0.5, 0.25], [1.0, 1.0, -1.0]),
                ([-0.5, 0.5, 0.25], [1.0, -1.0, -1.0]),
                ([0.5, -0.5, 0.25], [-1.0, 1.0, -1.0]),
                ([0.0, 0.5, 0.25], [-1.0, -1.0, -1.0]),
                ([0.5, 0.0, 0.25], [-1.0, -1.0, -1.0]),
            ]
            sad_oct_dirs = [[-1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [1.0, -1.0, -1.0]]
            v_rays = [np.array(vc)[None, :] + (self.depths * s2)[:, None] * np.array(vd)[None, :]
                      for vc, vd in salient_corners] + [
                      V_SAD[None, :] + (self.depths * s2)[:, None] * np.array(vd)[None, :]
                      for vd in sad_oct_dirs]
            v_pts = np.vstack(v_rays)
            for idx_d, dj in enumerate(self.depths):
                pts_dj = v_pts[idx_d::len(self.depths)]
                assert np.all(clearance_L(pts_dj) >= 0.7 * dj), f"vertex ray clearance < 0.7*d_j at d={dj}"

            self.pe = np.ascontiguousarray(np.vstack((base, v_pts)))
            self.pm = np.ascontiguousarray(v_pts)

            # Edge definitions on z >= 0:
            #   (mid_edge, t_hat, b_hat, H_e, is_vertical, has_meixner_tail)
            edge_specs = [
                # 5 salient 90-deg vertical edges (z in [-0.25, 0.25], H_e = 0.25, z_mid = 0)
                ([-0.5, -0.5, 0.0], [0., 0., 1.], [s2, s2, 0.], 0.25, True, True),
                ([-0.5, 0.5, 0.0], [0., 0., 1.], [s2, -s2, 0.], 0.25, True, True),
                ([0.5, -0.5, 0.0], [0., 0., 1.], [-s2, s2, 0.], 0.25, True, True),
                ([0.0, 0.5, 0.0], [0., 0., 1.], [-s2, -s2, 0.], 0.25, True, True),
                ([0.5, 0.0, 0.0], [0., 0., 1.], [-s2, -s2, 0.], 0.25, True, True),
                # 1 re-entrant 270-deg vertical edge (0, 0, z) -> interior diagonal sheet (-1,-1,0)/sqrt(2), gamma=2 (smooth)
                ([0.0, 0.0, 0.0], [0., 0., 1.], [-s2, -s2, 0.], 0.25, True, False),
                # 6 top horizontal salient 90-deg edges on z = +0.25 (2 of length 1.0, 4 of length 0.5)
                ([-0.5, 0.0, 0.25], [0., 1., 0.], [s2, 0., -s2], 0.50, False, True),
                ([0.0, -0.5, 0.25], [1., 0., 0.], [0., s2, -s2], 0.50, False, True),
                ([-0.25, 0.5, 0.25], [1., 0., 0.], [0., -s2, -s2], 0.25, False, True),
                ([0.5, -0.25, 0.25], [0., 1., 0.], [-s2, 0., -s2], 0.25, False, True),
                ([0.0, 0.25, 0.25], [0., 1., 0.], [-s2, 0., -s2], 0.25, False, True),
                ([0.25, 0.0, 0.25], [1., 0., 0.], [0., -s2, -s2], 0.25, False, True),
            ]

            self.lines = []
            r_vert = np.array([0], np.int64)
            r_horiz = np.array([0, 1], np.int64)
            last_by_edge = {}
            for idx_d, dj in enumerate(self.depths):
                for idx_e, (mid_e, th, bh, H_e, is_v, has_tail) in enumerate(edge_specs):
                    if (not has_tail) and (idx_d % 2 == 1):
                        continue
                    mid_e, th, bh = np.array(mid_e), np.array(th), np.array(bh)
                    L_e = 2.0 * H_e
                    lam = float(np.log((L_e - c_end * dj) / (c_end * dj)))
                    n = max(3, int(round(L_e * n0))) + int(np.ceil(2.0 * np.log(L_e / dj)))
                    F = np.ascontiguousarray(np.vstack((th, bh, np.cross(th, bh))))
                    mids = np.ascontiguousarray((mid_e + dj * bh)[None, :])
                    t_max = H_e * np.tanh(0.5 * lam)
                    ts_chk = np.linspace(-t_max, t_max, 9)[:, None] * th[None, :] + mids
                    min_cl = float(clearance_L(ts_chk).min())
                    assert min_cl >= 0.7 * dj, f"line clearance {min_cl:.6f} < 0.7*{dj:.6f} on edge {mid_e}"
                    entry = (H_e, th, mids, lam, float(dj), n, F, np.ones((2, 1, 1, 3)), False,
                             r_vert if is_v else r_horiz, is_v)
                    self.lines.append(entry)
                    last_by_edge[idx_e] = entry

            if tails:
                d_min = float(self.depths[-1])
                r_st = 0.5 ** np.arange(1, 10)
                w_meix = np.zeros((2, 1, 9, 3))
                for ki, kind in enumerate(("e", "m")):
                    ev = np.array(TAIL_E[kind] + (2.0,))
                    wv = r_st[:, None] ** ev[None, :]
                    wv[-1] /= 1.0 - 0.5 ** ev
                    w_meix[ki, 0] = wv
                w_meix = np.ascontiguousarray(w_meix)

                for idx_e, (mid_e, th, bh, H_e, is_v, has_tail) in enumerate(edge_specs):
                    if not has_tail:
                        continue
                    mid_e, th, bh = np.array(mid_e), np.array(th), np.array(bh)
                    _, _, _, lam, _, n, F, _, _, r_act, _ = last_by_edge[idx_e]
                    mids = np.ascontiguousarray(mid_e[None, :] + (d_min * r_st)[:, None] * bh[None, :])
                    for r_sub, m_sub in zip(r_st, mids):
                        t_max = H_e * np.tanh(0.5 * lam)
                        cl_s = float(clearance_L(np.vstack((m_sub - t_max * th, m_sub + t_max * th))).min())
                        assert cl_s >= 0.7 * d_min * r_sub * 0.999
                    self.lines.append((H_e, th, mids, lam, d_min, n, F, w_meix, True, r_act, is_v))

            self.v_tails = []
            if tails:
                d_min = float(self.depths[-1])
                r_st = 0.5 ** np.arange(1, 10)
                # 5 salient cube corners: 9-point radial tails weighted by NU_CUBE + 2, TAU_CUBE + 2
                ee_c, em_c = np.array([NU_CUBE + 2.0]), np.array([TAU_CUBE + 2.0])
                we_c = r_st[None, :] ** ee_c[:, None]; we_c[:, -1] /= 1.0 - 0.5 ** ee_c
                wm_c = r_st[None, :] ** em_c[:, None]; wm_c[:, -1] /= 1.0 - 0.5 ** em_c
                rpts_sal = np.ascontiguousarray([
                    np.array(vc)[None, :] + (d_min * s2 * r_st)[:, None] * np.array(vd)[None, :]
                    for vc, vd in salient_corners
                ])
                assert np.all(clearance_L(rpts_sal.reshape(-1, 3)) > 0.0)
                self.v_tails.append((rpts_sal, np.ascontiguousarray(we_c), np.ascontiguousarray(wm_c)))

                if sad_tail:
                    # 10-ray conical radial tail inside the 3 conducting octants at V_SAD = (0, 0, 0.25):
                    # 1 central spine ray (-1,-1,-1)/sqrt(3) + 3 skeleton-wall rays along each of the 3 edges meeting at V_SAD
                    sad_rays = [cone.unit(np.array([-1.0, -1.0, -1.0]))]
                    sad_edges = [
                        (np.array([1.0, 0.0, 0.0]), np.array([0.0, -s2, -s2])),   # +x notch edge -> B10
                        (np.array([0.0, 1.0, 0.0]), np.array([-s2, 0.0, -s2])),   # +y notch edge -> B01
                        (np.array([0.0, 0.0, -1.0]), np.array([-s2, -s2, 0.0])),  # -z notch edge -> B00
                    ]
                    for t_out, b_in in sad_edges:
                        for fr in (0.8, 0.4, 0.2):
                            sad_rays.append(cone.unit(t_out + fr * b_in))
                    assert len(sad_rays) == 10
                    sad_rays = np.array(sad_rays)
                    assert np.all((sad_rays[:, 2] < -1e-6) & ((sad_rays[:, 0] < -1e-6) | (sad_rays[:, 1] < -1e-6)))
                    ee_s, em_s = np.array([self.nu_sad + 2.0]), np.array([self.tau_sad + 2.0])
                    we_s = r_st[None, :] ** ee_s[:, None]; we_s[:, -1] /= 1.0 - 0.5 ** ee_s
                    wm_s = r_st[None, :] ** em_s[:, None]; wm_s[:, -1] /= 1.0 - 0.5 ** em_s
                    rpts_sad = np.ascontiguousarray([
                        V_SAD[None, :] + (1.2 * d_min * r_st)[:, None] * rd[None, :] for rd in sad_rays
                    ])
                    assert np.all(clearance_L(rpts_sad.reshape(-1, 3)) > 0.0)
                    self.v_tails.append((rpts_sad, np.ascontiguousarray(we_s), np.ascontiguousarray(wm_s)))

        # Column map in class CLS = (1, 1, -1)
        col = 3 * (len(self.pe) + len(self.pm))
        self.line_cmaps = []
        for H_e, th, mids, lam, dp, n, F, w_pair, is_tail, r_act, is_v in self.lines:
            cm = np.full((1, 6 * n), -1, np.int64)
            for kind in range(2):
                sc = np.array([-1.0, 1.0, 1.0]) * (-1.0 if kind == 1 else 1.0)
                for kk in range(n):
                    par = (-1.0) ** kk
                    for cc in range(2 if is_tail else 3):
                        if (not is_v) or (CLS[2] * par * sc[cc] > 0):
                            cm[0, (kind * n + kk) * 3 + cc] = col
                            col += 1
            self.line_cmaps.append(cm)
        self.vt_col0 = col
        for rpts, we, wm in self.v_tails:
            col += len(rpts) * (len(we) + len(wm)) * 3
        self.n_unknowns = col

        self.x, self.T, self.rw = build_collocation(d_coll, self.nu_sad, ng=4, hmax=0.25)
        self.rows = 2 * len(self.x)

    def assemble_into(self, x_sub, Tw_sub, A_buf):
        A_buf[:] = 0.0
        if len(self.pe) or len(self.pm):
            lines.points_add(x_sub, Tw_sub, REFL, self.pe, self.pm, self.k, A_buf, np.zeros(1, np.int64), SGN)
        for (H_e, th, mids, lam, dp, n, F, w_pair, is_tail, r_act, _), cm in zip(self.lines, self.line_cmaps):
            lines.line_add(x_sub, Tw_sub, REFL, r_act, th, mids, lam, dp, n, self.k, H_e, A_buf, cm, SGN,
                           rot=F, w=w_pair if is_tail else None)
        c_vt = self.vt_col0
        for rpts, we, wm in self.v_tails:
            nc_e = len(rpts) * len(we) * 3
            nc_m = len(rpts) * len(wm) * 3
            lines.radial_tails_add(x_sub, Tw_sub, REFL, rpts, we, wm, self.k, A_buf,
                                   np.array([c_vt], np.int64), np.array([c_vt + nc_e], np.int64), SGN)
            c_vt += nc_e + nc_m

    def assemble(self):
        Tw = np.ascontiguousarray(self.T * self.rw[:, None, None])
        A = np.moveaxis(np.zeros((self.rows, self.n_unknowns, 1), np.complex128, order="F"), 2, 0)
        self.assemble_into(self.x, Tw, A)
        return A[0]

    def _unpack_points(self, coef):
        ne, nm = len(self.pe), len(self.pm)
        y_list = [self.pe] if ne else []
        a_list = [coef[:3 * ne].reshape(1, ne, 3)] if ne else []
        b_pt = np.zeros((1, ne, 3), np.complex128)
        if nm:
            b_pt[0, ne - nm:] = coef[3 * ne:3 * (ne + nm)].reshape(nm, 3)
        b_list = [b_pt] if ne else []
        c_vt = self.vt_col0
        for rpts, we, wm in self.v_tails:
            R_cnt, nstep = rpts.shape[0], rpts.shape[1]
            nc_e, nc_m = R_cnt * len(we) * 3, R_cnt * len(wm) * 3
            y_list.append(rpts.reshape(-1, 3))
            ce = coef[c_vt:c_vt + nc_e].reshape(1, R_cnt, len(we), 3)
            cm = coef[c_vt + nc_e:c_vt + nc_e + nc_m].reshape(1, R_cnt, len(wm), 3)
            a_list.append(np.einsum("lrec,et->lrtc", ce, we).reshape(1, R_cnt * nstep, 3))
            b_list.append(np.einsum("lrmc,mt->lrtc", cm, wm).reshape(1, R_cnt * nstep, 3))
            c_vt += nc_e + nc_m
        if not y_list:
            return np.zeros((0, 3)), np.zeros((1, 0, 3), complex), np.zeros((1, 0, 3), complex)
        return (np.ascontiguousarray(np.vstack(y_list)),
                np.ascontiguousarray(np.concatenate(a_list, axis=1)),
                np.ascontiguousarray(np.concatenate(b_list, axis=1)))

    def field(self, coef, pts):
        pts = np.ascontiguousarray(pts, dtype=float)
        if not len(pts):
            return np.zeros((0, 3), np.complex128)
        out = np.zeros((1, len(pts), 3), np.complex128)
        y_pt, al, bt = self._unpack_points(coef)
        if len(y_pt):
            lines.eval_points_field(pts, REFL, SGN, y_pt, al, bt, self.k, out)
        for (H_e, th, mids, lam, dp, n, F, w_pair, is_tail, r_act, _), cm in zip(self.lines, self.line_cmaps):
            c_pair = np.zeros((1, 6 * n), np.complex128)
            act = cm[0] >= 0
            c_pair[0, act] = coef[cm[0, act]]
            lines.eval_line_field(pts, REFL, r_act, SGN, th, mids, lam, dp, n, self.k, H_e,
                                  c_pair.reshape(1, 2, 1, n, 3), out, rot=F, w=w_pair if is_tail else None)
        return out[0]

    def far(self, coef, dirs):
        dirs = np.ascontiguousarray(dirs, dtype=float)
        if not len(dirs):
            return np.zeros((0, 3), np.complex128)
        out = np.zeros((1, len(dirs), 3), np.complex128)
        y_pt, al, bt = self._unpack_points(coef)
        if len(y_pt):
            lines.eval_far_dipoles(dirs, REFL, np.array([0, 1], np.int64), SGN, y_pt, al, bt, self.k, out)
        for is_v_grp, r_act_grp in ((True, np.array([0], np.int64)), (False, np.array([0, 1], np.int64))):
            ys_l, al_l, bt_l = [], [], []
            for (H_e, th, mids, lam, dp, n, F, w_pair, _, _, is_v), cm in zip(self.lines, self.line_cmaps):
                if is_v != is_v_grp:
                    continue
                c_pair = np.zeros((1, 6 * n), np.complex128)
                act = cm[0] >= 0
                c_pair[0, act] = coef[cm[0, act]]
                tq, W, off, *_ = lines.build_tree(lam, dp, n, H_e, self.k)
                pm_n = int(off[1])
                c_cart = 2.0 * np.einsum("ksic,lksmc,cq->lkimq", w_pair, c_pair.reshape(1, 2, 1, n, 3), F)
                M0 = np.einsum("pm,lkimq->lkpiq", W[:pm_n], c_cart).reshape(1, 2, -1, 3)
                ys_l.append((mids[None, :, :] + tq[:pm_n, None, None] * th[None, None, :]).reshape(-1, 3))
                al_l.append(M0[:, 0])
                bt_l.append(M0[:, 1])
            if ys_l:
                lines.eval_far_dipoles(dirs, REFL, r_act_grp, SGN, np.ascontiguousarray(np.vstack(ys_l)),
                                       np.ascontiguousarray(np.concatenate(al_l, axis=1)),
                                       np.ascontiguousarray(np.concatenate(bt_l, axis=1)), self.k, out)
        return out[0]


def solve_mfs(family, J, nu_sad, tau_sad, d=D_W, p=P_W, pts=PTS, k=K):
    """Assemble, factorise, and solve one MFS configuration on Omega_L in class CLS = (1, 1, -1)."""
    t0 = time.time()
    prob = LShapeMFS(family, J, nu_sad, tau_sad, k=k)
    assert core.classes_of(k, d, p, REFL) == [CLS]
    rw2 = np.repeat(prob.rw, 2)
    b = -np.einsum("mab,mb->ma", prob.T, core.class_part(prob.x, k, d, p, REFL, CLS)).reshape(-1) * rw2
    b2 = float(np.linalg.norm(b) ** 2)
    M_pts, N = len(prob.x), prob.n_unknowns

    if 16.0 * prob.rows * N <= 1.8 * 2 ** 30:
        ta0 = time.time()
        A = prob.assemble()
        t_asm = time.time() - ta0
        tf0 = time.time()
        fac = core.factor(A, overwrite=True)
        coef, res = fac.solve(b, residual=True)
        r_norm = float(res)
        t_fac = time.time() - tf0
        del A, fac
    else:
        # Row-chunked Gram accumulation (matches mfs.tetrahedron.solve_scalene; peak RSS < 5 GiB at J = 6)
        Tw = np.ascontiguousarray(prob.T * prob.rw[:, None, None])
        G = np.zeros((N, N), np.complex128, order="F")
        rhs = np.zeros(N, np.complex128)
        n_chk = max(2, int(np.ceil(M_pts / 4500.0)))
        splits = np.linspace(0, M_pts, n_chk + 1, dtype=int)
        A_buf = np.moveaxis(np.zeros((2 * int(np.max(np.diff(splits))), N, 1), np.complex128, order="F"), 2, 0)
        ta0, tf_acc = time.time(), 0.0
        for ic in range(n_chk):
            m0, m1 = splits[ic], splits[ic + 1]
            sr = 2 * (m1 - m0)
            A_sub = A_buf[:, :sr, :] if sr == A_buf.shape[1] else np.moveaxis(
                np.zeros((sr, N, 1), np.complex128, order="F"), 2, 0)
            prob.assemble_into(np.ascontiguousarray(prob.x[m0:m1]), np.ascontiguousarray(Tw[m0:m1]), A_sub)
            tf1 = time.time()
            blas.zherk(1.0, A_sub[0], beta=1.0, c=G, trans=2, lower=0, overwrite_c=1)
            rhs += blas.zgemv(1.0, A_sub[0], b[2 * m0:2 * m1], trans=2)
            tf_acc += time.time() - tf1
        del A_buf
        gc.collect()
        t_asm, tf2 = (time.time() - ta0) - tf_acc, time.time()
        scale = np.sqrt(np.real(np.diag(G)))
        G /= scale[:, None]
        G /= scale[None, :]
        R = sla.cho_factor(G, overwrite_a=True, check_finite=False)
        coef = sla.cho_solve(R, rhs / scale, check_finite=False) / scale
        r_norm = float(np.sqrt(max(0.0, b2 - float(np.real(np.vdot(coef, rhs))))))
        del R, G
        t_fac = tf_acc + (time.time() - tf2)
    gc.collect()

    E = None if pts is None else prob.field(coef, pts)
    t_tot = time.time() - t0
    ext, sca = core.cross_sections(lambda dirs: prob.far(coef, dirs), k, d, p, RULE)
    info = dict(
        unknowns=prob.n_unknowns, rows=prob.rows, t_assemble=t_asm, t_factor=t_fac,
        residual=r_norm / np.sqrt(b2), ext=ext, sca=sca, optical=abs(ext - sca) / sca,
    )
    return E, t_tot, info, prob, coef


# ---------------------------------------------------------------------------
# 4. Watertight Outward-Oriented Graded Bempp-cl Mesh & z-Mirror Solver
# ---------------------------------------------------------------------------

def l_shape_grid(m, g=1.5):
    """
    Watertight, outward-oriented edge-graded triangle mesh of Omega_L symmetric under z -> -z,
    built from the z >= 0 half-surface (m intervals per 0.25 quarter-unit, grading exponent g)
    via mfs.bem.mirror_grid(..., REFL).
      Full mesh counts: F = 112 m^2, E = N_u = 168 m^2, V = 56 m^2 + 2, V - E + F = 2, vol = +0.375.
    """
    from mfs import bem
    m = int(m)
    sh = 0.25 * (1.0 - (1.0 - np.arange(m + 1) / m) ** float(g))          # [0, 0.25], graded towards 0.25
    s2 = np.concatenate((0.25 - sh[::-1], 0.25 + sh[1:]))                  # [0, 0.50], graded towards 0 and 0.5
    s_neg, s_pos = -0.5 + s2, s2

    v_list, t_list = [], []
    v_off = 0

    def add_grid_patch(u_nodes, v_nodes, map_fn, n_out):
        nonlocal v_off
        nu, nv = len(u_nodes) - 1, len(v_nodes) - 1
        U, V = np.meshgrid(u_nodes, v_nodes, indexing="ij")
        pts = map_fn(U.ravel(), V.ravel())
        idx = lambda i, j: i * (nv + 1) + j
        tris = []
        for i in range(nu):
            for j in range(nv):
                p00, p10, p11, p01 = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
                if np.dot(np.cross(pts[p10] - pts[p00], pts[p11] - pts[p00]), n_out) > 0:
                    tris += [(p00, p10, p11), (p00, p11, p01)]
                else:
                    tris += [(p00, p11, p10), (p00, p01, p11)]
        v_list.append(pts)
        t_list.append(np.array(tris, dtype=np.int64) + v_off)
        v_off += len(pts)

    # 3 top 0.5 x 0.5 squares at z = +0.25
    add_grid_patch(s_neg, s_neg, lambda u, v: np.column_stack((u, v, np.full_like(u, 0.25))), [0., 0., 1.])
    add_grid_patch(s_neg, s_pos, lambda u, v: np.column_stack((u, v, np.full_like(u, 0.25))), [0., 0., 1.])
    add_grid_patch(s_pos, s_neg, lambda u, v: np.column_stack((u, v, np.full_like(u, 0.25))), [0., 0., 1.])

    # 8 vertical 0.5 x 0.25 half-faces on z in [0, 0.25]
    add_grid_patch(s_neg, sh, lambda y, z: np.column_stack((np.full_like(y, -0.5), y, z)), [-1., 0., 0.])
    add_grid_patch(s_pos, sh, lambda y, z: np.column_stack((np.full_like(y, -0.5), y, z)), [-1., 0., 0.])
    add_grid_patch(s_neg, sh, lambda x, z: np.column_stack((x, np.full_like(x, -0.5), z)), [0., -1., 0.])
    add_grid_patch(s_pos, sh, lambda x, z: np.column_stack((x, np.full_like(x, -0.5), z)), [0., -1., 0.])
    add_grid_patch(s_neg, sh, lambda y, z: np.column_stack((np.full_like(y, 0.5), y, z)), [1., 0., 0.])
    add_grid_patch(s_neg, sh, lambda x, z: np.column_stack((x, np.full_like(x, 0.5), z)), [0., 1., 0.])
    add_grid_patch(s_pos, sh, lambda y, z: np.column_stack((np.full_like(y, 0.0), y, z)), [1., 0., 0.])
    add_grid_patch(s_pos, sh, lambda x, z: np.column_stack((x, np.full_like(x, 0.0), z)), [0., 1., 0.])

    grid, vmap, n_piece = bem.mirror_grid(np.vstack(v_list), np.vstack(t_list), REFL)
    return grid, vmap, n_piece


def verify_mesh(grid):
    """Verify watertightness (0 boundary edges, paired opposite half-edges), V - E + F = 2, and vol = +0.375."""
    vt, el = grid.vertices.T, grid.elements.T
    V, F = len(vt), len(el)
    he = {}
    for a, b, c in el.tolist():
        for u, v in ((a, b), (b, c), (c, a)):
            assert (u, v) not in he, f"duplicate directed half-edge {(u, v)}"
            he[(u, v)] = 1
    n_bdy = sum(1 for (u, v) in he if (v, u) not in he)
    E = len(he) // 2
    euler = V - E + F
    p0, p1, p2 = vt[el[:, 0]], vt[el[:, 1]], vt[el[:, 2]]
    vol = float(np.sum(np.einsum("ij,ij->i", p0, np.cross(p1, p2))) / 6.0)
    assert n_bdy == 0 and euler == 2 and abs(vol - 0.375) < 1e-13, (V, E, F, euler, n_bdy, vol)
    return V, E, F, euler, n_bdy, vol


def solve_bempp_l(m, g=1.5, reg=6, sing=8, d=D_W, p=P_W, pts=PTS, k=K):
    """Solve EFIE on Omega_L with exact z -> -z mirror symmetry (mfs.bem.Mirror, single class CLS)."""
    from mfs import bem
    bempp = bem.bempp
    bempp.GLOBAL_PARAMETERS.quadrature.regular = int(reg)
    bempp.GLOBAL_PARAMETERS.quadrature.singular = int(sing)
    grid, vmap, n_piece = l_shape_grid(m, g)
    V, E_cnt, F_cnt, euler, n_bdy, vol = verify_mesh(grid)
    gc.collect()
    t0 = time.time()
    sym = bem.Mirror(grid, vmap, n_piece, REFL, k)
    t_asm = sym.t_assemble
    ts0 = time.time()
    J_rwg = sym.solve(d, p)
    t_lu = time.time() - ts0
    n_cls = len(sym.lu[CLS][0])
    sym.lu.clear()
    sym.A = None
    gc.collect()
    E_scat = sym.field(J_rwg, pts)
    t_tot = time.time() - t0

    # Optical theorem check via Bempp exact far-field potential operator (outside wall timer)
    lam_gf = bempp.GridFunction(sym.div, coefficients=J_rwg)
    far_fn = lambda dirs: -(bempp.operators.far_field.maxwell.electric_field(
        sym.div, np.ascontiguousarray(dirs.T), k) * lam_gf).T
    ext, sca = core.cross_sections(far_fn, k, d, p, RULE)
    opt = abs(ext - sca) / sca
    return E_scat, t_tot, dict(
        m=m, g=g, reg=reg, sing=sing, N_full=sym.N, N_cls=n_cls, V=V, E=E_cnt, F=F_cnt,
        vol=vol, t_assemble=t_asm, t_solve=t_lu, ext=ext, sca=sca, optical=opt,
    )


# ---------------------------------------------------------------------------
# 5. CLI Benchmark Driver & Summary
# ---------------------------------------------------------------------------

def load_st():
    return dict(np.load(NPZ)) if os.path.exists(NPZ) else {"pts": PTS, "d": D_W, "p": P_W, "k": np.array(K)}


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def run_part(part):
    st = load_st()

    if part in ("all", "cone"):
        t0 = time.time()
        (nu1, sn1), (tau1, st1), np1 = saddle_cone_exponents(11)
        (nu2, sn2), (tau2, st2), np2 = saddle_cone_exponents(15)
        d_nu, d_tau = abs(nu2 - nu1), abs(tau2 - tau1)
        assert 0.45 < nu2 < 1.0 and 0.45 < tau2 < 1.0 and d_nu < 5e-5 and d_tau < 1e-4
        st["cone_poles"] = np.array([[np1, nu1, sn1, tau1, st1], [np2, nu2, sn2, tau2, st2]], float)
        st["nu_sad"] = np.array(nu2)
        st["tau_sad"] = np.array(tau2)
        np.savez(NPZ, **st)
        print(f"[cone] ({time.time() - t0:.1f}s) poles={np1}: nu_sad={nu1:.6f} (sig={sn1:.2e}), tau_sad={tau1:.6f} (sig={st1:.2e})")
        print(f"[cone]          poles={np2}: nu_sad={nu2:.6f} (sig={sn2:.2e}), tau_sad={tau2:.6f} (sig={st2:.2e})")
        print(f"[cone] stability: |d_nu|={d_nu:.2e}, |d_tau|={d_tau:.2e}", flush=True)
        if part == "cone":
            return

    nu_sad, tau_sad = float(st["nu_sad"]), float(st["tau_sad"])

    if part in ("all", "sanity"):
        from mfs import bem
        # 1. Mesh topology & Bempp z-mirror equivalence against full unsymmetrised EFIE
        grid3, _, _ = l_shape_grid(3, 1.5)
        V, E_c, F_c, euler, n_bdy, vol = verify_mesh(grid3)
        pts_chk = PTS[:40]
        bem.bempp.GLOBAL_PARAMETERS.quadrature.regular, bem.bempp.GLOBAL_PARAMETERS.quadrature.singular = 4, 6
        E_b_sym, _, ib3 = solve_bempp_l(3, 1.5, 4, 6, pts=pts_chk)
        E_b_full, N_b_full = bem.efie(grid3, K, D_W, P_W, pts_chk)
        sym_vs_full = core.rel(E_b_sym, E_b_full)
        print(f"[sanity] Mesh m=3: V={V}, E={E_c}, F={F_c}, V-E+F={euler}, bdy={n_bdy}, vol={vol:.6f}")
        print(f"[sanity] Bempp Mirror (N_cls={ib3['N_cls']}) vs full EFIE (N={N_b_full}): rel_diff={sym_vs_full:.2e}, opt={ib3['optical']:.2e}")
        assert sym_vs_full < 5e-5

        # 2. MFS far-field contraction check & coarse MFS vs Bempp agreement
        E2, t2, i2, p2, c2 = solve_mfs("full", 2, nu_sad, tau_sad, pts=pts_chk)
        E3, t3, i3, p3, c3 = solve_mfs("full", 3, nu_sad, tau_sad, pts=pts_chk)
        dirs_chk = core.fibonacci(25)
        r_inf = 1e6
        far_rel = core.rel(p2.far(c2, dirs_chk), r_inf * np.exp(-1j * K * r_inf) * p2.field(c2, r_inf * dirs_chk))
        print(f"[sanity] MFS full J=2: N_u={i2['unknowns']}, M={i2['rows']}, t={t2:.2f}s, opt={i2['optical']:.2e}, far_chk={far_rel:.2e}")
        print(f"[sanity] MFS full J=3: N_u={i3['unknowns']}, M={i3['rows']}, t={t3:.2f}s, opt={i3['optical']:.2e}, "
              f"diff_vs_J2={core.rel(E2, E3):.2e}, diff_vs_bem3={core.rel(E3, E_b_sym):.2e}", flush=True)
        assert far_rel < 2e-5
        st["mesh_check"] = np.array([V, E_c, F_c, euler, n_bdy, vol, sym_vs_full], float)
        np.savez(NPZ, **st)
        if part == "sanity":
            return

    def _warmup_mfs():
        _p = LShapeMFS("full", 1, nu_sad, tau_sad)
        _a = np.moveaxis(np.zeros((4, _p.n_unknowns, 1), np.complex128, order="F"), 2, 0)
        _p.assemble_into(_p.x[:2], np.ascontiguousarray(_p.T[:2]), _a)
        _ = _p.field(np.zeros(_p.n_unknowns, np.complex128), PTS[:2])
        _ = _p.far(np.zeros(_p.n_unknowns, np.complex128), PTS[:2])
        del _p, _a
        gc.collect()

    if part in ("all", "ref"):
        _warmup_mfs()
        E_ref, t_ref, info_ref, _, _ = solve_mfs("full", 6, nu_sad, tau_sad, pts=PTS)
        st["ref_E"] = E_ref
        st["ref_info"] = np.array([
            6, info_ref["unknowns"], info_ref["rows"], t_ref,
            info_ref["t_assemble"], info_ref["t_factor"],
            info_ref["residual"], info_ref["optical"], info_ref["ext"], info_ref["sca"],
        ], float)
        np.savez(NPZ, **st)
        print(f"[ref] MFS full J=6: N_u={info_ref['unknowns']}, M={info_ref['rows']}, t={t_ref:.2f}s "
              f"(asm={info_ref['t_assemble']:.2f}s, fac={info_ref['t_factor']:.2f}s), "
              f"res={info_ref['residual']:.2e}, opt={info_ref['optical']:.2e}, ext={info_ref['ext']:.6f}, "
              f"RSS={rss_mib():.0f} MiB", flush=True)
        if part == "ref":
            return

    E_ref = st["ref_E"]

    if part in ("all", "mfs", "mfs1"):
        _warmup_mfs()

        # Family 1: eroded inner L-surface (margin 0.07), N_u from ~380 to ~1550
        inner_prev = {int(r[0]): list(r) for r in st.get("mfs_inner", np.zeros((0, 7)))}
        for n_div in (5, 7, 8, 10, 11):
            if n_div in inner_prev and f"E_inner_{n_div}" in st:
                continue
            gc.collect()
            Ei, ti, ii, _, _ = solve_mfs("inner", n_div, nu_sad, tau_sad, pts=PTS)
            err_i = core.rel(Ei, E_ref)
            inner_prev[n_div] = [n_div, ii["unknowns"], ii["rows"], ti, ii["residual"], err_i, ii["optical"]]
            st[f"E_inner_{n_div}"] = Ei
            st["mfs_inner"] = np.array([inner_prev[k] for k in sorted(inner_prev)], float)
            np.savez(NPZ, **st)
            print(f"[mfs] inner n_div={n_div:2d}: N_u={ii['unknowns']:4d}, M={ii['rows']:5d}, t={ti:5.2f}s, "
                  f"res={ii['residual']:.2e}, err={err_i:.3e}, opt={ii['optical']:.2e}", flush=True)

        # Family 2: 'notail' (J = 2..5)
        notail_prev = {int(r[0]): list(r) for r in st.get("mfs_notail", np.zeros((0, 7)))}
        for J in (2, 3, 4, 5):
            if J in notail_prev and f"E_notail_J{J}" in st:
                continue
            gc.collect()
            En, tn, inf_n, _, _ = solve_mfs("notail", J, nu_sad, tau_sad, pts=PTS)
            err_n = core.rel(En, E_ref)
            notail_prev[J] = [J, inf_n["unknowns"], inf_n["rows"], tn, inf_n["residual"], err_n, inf_n["optical"]]
            st[f"E_notail_J{J}"] = En
            st["mfs_notail"] = np.array([notail_prev[k] for k in sorted(notail_prev)], float)
            np.savez(NPZ, **st)
            print(f"[mfs] notail J={J}: N_u={inf_n['unknowns']:5d}, M={inf_n['rows']:5d}, t={tn:6.2f}s, "
                  f"res={inf_n['residual']:.2e}, err={err_n:.3e}, opt={inf_n['optical']:.2e}", flush=True)
        if part == "mfs1":
            return

    if part in ("all", "mfs", "mfs2", "mfs3"):
        if part in ("mfs2", "mfs3"):
            _warmup_mfs()
        j_list = (2, 3, 4) if part == "mfs2" else ((5,) if part == "mfs3" else (2, 3, 4, 5))
        full_prev = {int(r[0]): list(r) for r in st.get("mfs_full", np.zeros((0, 7)))}
        for J in j_list:
            if J in full_prev and f"E_full_J{J}" in st:
                continue
            gc.collect()
            Ef, tf, inf_f, _, _ = solve_mfs("full", J, nu_sad, tau_sad, pts=PTS)
            err_f = core.rel(Ef, E_ref)
            full_prev[J] = [J, inf_f["unknowns"], inf_f["rows"], tf, inf_f["residual"], err_f, inf_f["optical"]]
            st[f"E_full_J{J}"] = Ef
            st["mfs_full"] = np.array([full_prev[k] for k in sorted(full_prev)], float)
            np.savez(NPZ, **st)
            print(f"[mfs] full   J={J}: N_u={inf_f['unknowns']:5d}, M={inf_f['rows']:5d}, t={tf:6.2f}s, "
                  f"res={inf_f['residual']:.2e}, err={err_f:.3e}, opt={inf_f['optical']:.2e}", flush=True)
        if part == "mfs2":
            return

        # Saddle-tail ablation at J = 5 ('no_saddle_tail')
        if "mfs_no_saddle_tail_J5" not in st:
            gc.collect()
            Ens, tns, inf_ns, _, _ = solve_mfs("no_saddle_tail", 5, nu_sad, tau_sad, pts=PTS)
            err_ns = core.rel(Ens, E_ref)
            st["E_no_saddle_tail_J5"] = Ens
            st["mfs_no_saddle_tail_J5"] = np.array([
                5, inf_ns["unknowns"], inf_ns["rows"], tns, inf_ns["residual"], err_ns, inf_ns["optical"],
            ], float)
            np.savez(NPZ, **st)
            print(f"[mfs] no_saddle_tail J=5: N_u={inf_ns['unknowns']:5d}, M={inf_ns['rows']:5d}, t={tns:6.2f}s, "
                  f"res={inf_ns['residual']:.2e}, err={err_ns:.3e}, opt={inf_ns['optical']:.2e} "
                  f"(RSS={rss_mib():.0f} MiB)", flush=True)
        if part in ("mfs", "mfs3"):
            return

    if part in ("all", "bem", "bem1", "bem2"):
        solve_bempp_l(2, 1.5, 6, 8, pts=PTS[:10])
        E_J5 = st["E_full_J5"]
        m_list = (4, 6, 8) if part == "bem1" else ((10,) if part == "bem2" else (4, 6, 8, 10))
        prev = {int(r[0]): list(r) for r in st.get("bem", np.zeros((0, 11)))}
        for m in m_list:
            if m in prev and f"E_bem_m{m}" in st:
                continue
            gc.collect()
            Eb, tb, ib = solve_bempp_l(m, 1.5, 6, 8, pts=PTS)
            diff_J6 = core.rel(Eb, E_ref)
            diff_J5 = core.rel(Eb, E_J5)
            prev[m] = [
                m, ib["g"], ib["reg"], ib["sing"], ib["N_full"], ib["N_cls"],
                tb, ib["t_assemble"], ib["t_solve"], diff_J6, ib["optical"],
            ]
            st[f"E_bem_m{m}"] = Eb
            st["bem"] = np.array([prev[k] for k in sorted(prev)], float)
            np.savez(NPZ, **st)
            print(f"[bem] m={m:2d} (6/8, g=1.5): N_full={ib['N_full']:5d} (N_cls={ib['N_cls']:5d}), "
                  f"t={tb:6.2f}s (asm={ib['t_assemble']:.1f}s, lu={ib['t_solve']:.1f}s), "
                  f"diff_vs_J6={diff_J6:.3e}, diff_vs_J5={diff_J5:.3e}, opt={ib['optical']:.2e}, "
                  f"RSS={rss_mib():.0f} MiB", flush=True)
        if part in ("bem", "bem1", "bem2"):
            return

    if part in ("all", "bem_q8"):
        solve_bempp_l(2, 1.5, 8, 10, pts=PTS[:10])
        gc.collect()
        Eb8, tb8, ib8 = solve_bempp_l(8, 1.5, 8, 10, pts=PTS)
        diff_J6 = core.rel(Eb8, E_ref)
        diff_J5 = core.rel(Eb8, st["E_full_J5"])
        st["E_bem_m8_q8"] = Eb8
        st["bem_q8"] = np.array([
            8, ib8["g"], ib8["reg"], ib8["sing"], ib8["N_full"], ib8["N_cls"],
            tb8, ib8["t_assemble"], ib8["t_solve"], diff_J6, ib8["optical"],
        ], float)
        np.savez(NPZ, **st)
        print(f"[bem_q8] m= 8 (8/10, g=1.5): N_full={ib8['N_full']:5d} (N_cls={ib8['N_cls']:5d}), "
              f"t={tb8:6.2f}s (asm={ib8['t_assemble']:.1f}s, lu={ib8['t_solve']:.1f}s), "
              f"diff_vs_J6={diff_J6:.3e}, diff_vs_J5={diff_J5:.3e}, opt={ib8['optical']:.2e}, "
              f"RSS={rss_mib():.0f} MiB", flush=True)


def show():
    d = np.load(NPZ)
    print("=== 3D Non-Convex L-Shaped Conducting Polyhedron (Omega_L, k = 2*pi, z -> -z symmetry) ===")
    cp = d["cone_poles"]
    for r in cp:
        print(f"  Saddle cone ({int(r[0])} poles): nu_sad = {r[1]:.6f} (sig={r[2]:.2e}), tau_sad = {r[3]:.6f} (sig={r[4]:.2e})")
    if "mesh_check" in d:
        mc = d["mesh_check"]
        print(f"  Mesh check (m=3): V={int(mc[0])}, E={int(mc[1])}, F={int(mc[2])}, V-E+F={int(mc[3])}, "
              f"bdy={int(mc[4])}, vol={mc[5]:.6f}, Mirror_vs_full={mc[6]:.2e}")
    if "ref_info" in d:
        ri = d["ref_info"]
        print(f"\n  MFS full reference (J={int(ri[0])}): N_u={int(ri[1])}, M={int(ri[2])}, t={ri[3]:.2f}s "
              f"(asm={ri[4]:.2f}s, fac={ri[5]:.2f}s), res={ri[6]:.2e}, opt={ri[7]:.2e}, ext={ri[8]:.6f}")
    for key, label in (("mfs_inner", "Family 1: inner (margin 0.07)"),
                       ("mfs_notail", "Family 2: notail (J = 2..5)"),
                       ("mfs_full", "Family 3: full (J = 2..5)")):
        if key in d:
            print(f"\n  --- {label} ---")
            for r in d[key]:
                tag = "n_div" if key == "mfs_inner" else "J"
                print(f"    {tag}={int(r[0]):2d}: N_u={int(r[1]):5d}, M={int(r[2]):5d}, t={r[3]:6.2f}s, "
                      f"res={r[4]:.2e}, err_vs_J6={r[5]:.3e}, opt={r[6]:.2e}")
            if key == "mfs_notail" and len(d[key]) >= 3:
                c = d[key]
                beta = -np.polyfit(c[:, 0], np.log2(c[:, 5]), 1)[0]
                print(f"    Algebraic fit (J=2..5): err ~ 2^(-{beta:.4f} J)  [Meixner theory: gamma = 2/3 = 0.6667]")
            if key == "mfs_full" and len(d[key]) >= 3:
                c = d[key]
                p = np.polyfit(np.sqrt(c[:, 1]), np.log(c[:, 5]), 1)
                print(f"    Root-exp fit (J=2..5): err = {np.exp(p[1]):.3g} * exp({p[0]:.4f} * sqrt(N_u))")
    if "mfs_no_saddle_tail_J5" in d:
        r = d["mfs_no_saddle_tail_J5"]
        d_sad = core.rel(d["E_no_saddle_tail_J5"], d["E_full_J5"]) if "E_full_J5" in d else np.nan
        print(f"\n  --- Ablation: no_saddle_tail (J = {int(r[0])}) ---")
        print(f"    J={int(r[0])}: N_u={int(r[1]):5d}, M={int(r[2]):5d}, t={r[3]:6.2f}s, "
              f"res={r[4]:.2e}, err_vs_J6={r[5]:.3e}, diff_vs_full_J5={d_sad:.3e}, opt={r[6]:.2e}")
    if "bem" in d:
        print("\n  --- Bempp-cl z-Mirror EFIE (g=1.5, quadrature 6/8 & 8/10 check) ---")
        E5 = d.get("E_full_J5", None)
        ms = [int(r[0]) for r in d["bem"]]
        for idx_r, r in enumerate(d["bem"]):
            m = int(r[0])
            d5_str = f", diff_vs_J5={core.rel(d[f'E_bem_m{m}'], E5):.3e}" if E5 is not None and f"E_bem_m{m}" in d else ""
            step_str = f", step={core.rel(d[f'E_bem_m{ms[idx_r-1]}'], d[f'E_bem_m{m}']):.3e}" if idx_r > 0 else ""
            print(f"    m={m:2d} (6/8): N_full={int(r[4]):5d} (N_cls={int(r[5]):5d}), t={r[6]:6.2f}s "
                  f"(asm={r[7]:.1f}s, lu={r[8]:.1f}s), diff_vs_J6={r[9]:.3e}{d5_str}{step_str}, opt={r[10]:.2e}")
    if "bem_q8" in d:
        r = d["bem_q8"]
        E5 = d.get("E_full_J5", None)
        d5_str = f", diff_vs_J5={core.rel(d['E_bem_m8_q8'], E5):.3e}" if E5 is not None and "E_bem_m8_q8" in d else ""
        q_shift = core.rel(d["E_bem_m8_q8"], d["E_bem_m8"]) if "E_bem_m8" in d else np.nan
        print(f"    m={int(r[0]):2d} (8/10): N_full={int(r[4]):5d} (N_cls={int(r[5]):5d}), t={r[6]:6.2f}s "
              f"(asm={r[7]:.1f}s, lu={r[8]:.1f}s), diff_vs_J6={r[9]:.3e}{d5_str}, q8_vs_q6={q_shift:.3e}, opt={r[10]:.2e}")


if __name__ == "__main__":
    if "--part" in sys.argv:
        run_part(sys.argv[sys.argv.index("--part") + 1])
    elif "--run" in sys.argv:
        run_part("all")
    else:
        show()
