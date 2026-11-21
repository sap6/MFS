"""
Scattering by the NASA almond (Section 6 of the paper): length L = 0.2524 m,
long axis along x with the tip at +x, mid-plane z = 0, cross-sections 3:1
ellipses. The continued field is singular near the focal curve, which runs
at the fraction ETA_F of the half-width in the mid-plane, so the sources are
electric and magnetic dipoles on the mid-plane in rows that approach it.
Near the tip each row is continued by a tail whose weights follow from the
exponents of the tip cone (Appendix A).

The almond is symmetric under y → -y and z → -z, so a wave splits into four
parity classes, each solved on a quarter of the body.
"""
import math
import numpy as np
from numba import njit, prange

from .core import class_part, cross_sections, factor, sphere_rule, tangent_basis

L = 0.2524                                                  # length [m]
T0, T1 = -0.416667, 0.58333                                 # rear end and tip, in units of L
A_, B_, C_ = 0.416667 * L, 0.193333 * L, 0.064444 * L       # semi-axes of the rear half-ellipsoid
ETA_F = np.sqrt(1 - (C_ / B_) ** 2)                         # focal curve, as a fraction of the half-width (0.943)
GAP = 1 - ETA_F                                             # from the focal curve to the rim, as a fraction of the half-width
C0 = 299792458.0                                            # speed of light [m/s]
REFL = np.array([[1.0, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1]])   # the images under y → -y and z → -z

# Exponents of the tip cone (half-angles 34.09 and 12.71 degrees) by the parities (y, z) of the
# potential: electric type E = ∇(r^nu Y), magnetic type H = ∇(r^tau Z) (mfs.cone, Table A.1)
CONE_NU = {(1, 1): (0.29959, 1.46427), (-1, 1): (1.10372, 2.24730), (1, -1): (1.04120, 2.10576), (-1, -1): (2.01295, 3.04849)}
CONE_TAU = {(1, 1): (1.05604, 1.98580), (-1, 1): (0.96803, 1.96545), (1, -1): (0.87905, 1.75454), (-1, -1): (1.97538, 2.90094)}

_COLLOC_CACHE = {}


def semi_axes(x):
    """Semi-axes (a, b) of the cross-section at x (metres)"""
    t = np.asarray(x, float) / L
    f1 = np.sqrt(np.clip(1 - (t / 0.416667) ** 2, 0, None))
    f2 = np.sqrt(np.clip(1 - (t / 2.08335) ** 2, 0, None)) - 0.96
    a = np.where(t < 0, 0.193333 * f1, 4.83345 * f2) * L
    b = np.where(t < 0, 0.064444 * f1, 1.61115 * f2) * L
    return np.clip(a, 0, None), np.clip(b, 0, None)


def midplane_comps(sz):
    """Dipole components that survive on z = 0 under z → -z (0-2: electric x, y, z; 3-5: magnetic x, y, z)"""
    return np.array([0, 1, 5]) if sz > 0 else np.array([2, 3, 4])


@njit(parallel=True, fastmath=True, cache=True)
def _block(x, T, y, comps, sy, k, out):
    """out[m * nc + a, 3 * n + i]: row a of T[m] against dipole comps[i] at y[n], assembled columnwise in F-order"""
    M, nc, N = x.shape[0], T.shape[1], y.shape[0]
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    c0, c1, c2 = comps[0], comps[1], comps[2]
    for n in prange(N):
        yn0, yn1, yn2 = y[n, 0], y[n, 1], y[n, 2]
        col = 3 * n
        for m in range(M):
            mt = m * tmul
            row = m * nc
            dx0_p = x[m, 0] - yn0
            dx1_p = x[m, 1] - yn1
            dx2_p = x[m, 2] - yn2
            r2_p = dx0_p * dx0_p + dx1_p * dx1_p + dx2_p * dx2_p
            ir_p = 1.0 / np.sqrt(r2_p)
            kr_p = k * (r2_p * ir_p)
            ik_p = ir_p * inv_k
            ik2_p = ik_p * ik_p
            ph_p = (np.cos(kr_p) + 1j * np.sin(kr_p)) * ir_p
            fI_p = (1.0 + 1j * ik_p - ik2_p) * ph_p
            fR_p = (1.0 + 3j * ik_p - 3.0 * ik2_p) * ph_p * (ir_p * ir_p)
            fm_p = (1j - ik_p) * ph_p * ir_p

            dx1_m = -x[m, 1] - yn1
            r2_m = dx0_p * dx0_p + dx1_m * dx1_m + dx2_p * dx2_p
            ir_m = 1.0 / np.sqrt(r2_m)
            kr_m = k * (r2_m * ir_m)
            ik_m = ir_m * inv_k
            ik2_m = ik_m * ik_m
            ph_m = (np.cos(kr_m) + 1j * np.sin(kr_m)) * ir_m
            fI_m = (1.0 + 1j * ik_m - ik2_m) * ph_m
            fR_m = (1.0 + 3j * ik_m - 3.0 * ik2_m) * ph_m * (ir_m * ir_m)
            fm_m = (1j - ik_m) * ph_m * ir_m
            sg_m = 2.0 * sy

            for a in range(nc):
                t0, t1, t2 = T[mt, a, 0], T[mt, a, 1], T[mt, a, 2]
                td_p = fR_p * (t0 * dx0_p + t1 * dx1_p + t2 * dx2_p)
                ge0 = 2.0 * (fI_p * t0 - td_p * dx0_p)
                ge1 = 2.0 * (fI_p * t1 - td_p * dx1_p)
                ge2 = 2.0 * (fI_p * t2 - td_p * dx2_p)
                gm0 = 2.0 * fm_p * (t1 * dx2_p - t2 * dx1_p)
                gm1 = 2.0 * fm_p * (t2 * dx0_p - t0 * dx2_p)
                gm2 = 2.0 * fm_p * (t0 * dx1_p - t1 * dx0_p)

                t1m = -t1
                td_m = fR_m * (t0 * dx0_p + t1m * dx1_m + t2 * dx2_p)
                ge0 += sg_m * (fI_m * t0 - td_m * dx0_p)
                ge1 += sg_m * (fI_m * t1m - td_m * dx1_m)
                ge2 += sg_m * (fI_m * t2 - td_m * dx2_p)
                gm0 += sg_m * fm_m * (t1m * dx2_p - t2 * dx1_m)
                gm1 += sg_m * fm_m * (t2 * dx0_p - t0 * dx2_p)
                gm2 += sg_m * fm_m * (t0 * dx1_m - t1m * dx0_p)

                out[row + a, col] = ge0 if c0 == 0 else (ge1 if c0 == 1 else (ge2 if c0 == 2 else (gm0 if c0 == 3 else (gm1 if c0 == 4 else gm2))))
                out[row + a, col + 1] = ge0 if c1 == 0 else (ge1 if c1 == 1 else (ge2 if c1 == 2 else (gm0 if c1 == 3 else (gm1 if c1 == 4 else gm2))))
                out[row + a, col + 2] = ge0 if c2 == 0 else (ge1 if c2 == 1 else (ge2 if c2 == 2 else (gm0 if c2 == 3 else (gm1 if c2 == 4 else gm2))))


@njit(parallel=True, fastmath=True, cache=True)
def _eval_midplane(pts, y, c_vec, comps, sy, k, out):
    """Direct evaluation of symmetrised mid-plane dipoles at pts (P, 3) without forming matrix A"""
    P, N = pts.shape[0], y.shape[0]
    inv_k = 1.0 / k
    c0, c1, c2 = comps[0], comps[1], comps[2]
    for m in prange(P):
        s0, s1, s2 = 0.0j, 0.0j, 0.0j
        for g in range(2):
            ry = 1.0 if g == 0 else -1.0
            sg = 2.0 if g == 0 else 2.0 * sy
            px0, px1, px2 = pts[m, 0], pts[m, 1] * ry, pts[m, 2]
            e0, e1, e2 = 0.0j, 0.0j, 0.0j
            for n in range(N):
                dx0 = px0 - y[n, 0]
                dx1 = px1 - y[n, 1]
                dx2 = px2 - y[n, 2]
                r2 = dx0 * dx0 + dx1 * dx1 + dx2 * dx2
                inv_r = 1.0 / np.sqrt(r2)
                kr = k * (r2 * inv_r)
                ik = inv_r * inv_k
                ik2 = ik * ik
                ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
                fI = (1.0 + 1j * ik - ik2) * ph
                fR = (1.0 + 3j * ik - 3.0 * ik2) * ph * (inv_r * inv_r)
                fm = (1j - ik) * ph * inv_r
                ax, ay, az, bx, by, bz = 0.0j, 0.0j, 0.0j, 0.0j, 0.0j, 0.0j
                for idx, cid in ((0, c0), (1, c1), (2, c2)):
                    v = c_vec[n, idx]
                    if cid == 0:
                        ax = v
                    elif cid == 1:
                        ay = v
                    elif cid == 2:
                        az = v
                    elif cid == 3:
                        bx = v
                    elif cid == 4:
                        by = v
                    else:
                        bz = v
                da = fR * (dx0 * ax + dx1 * ay + dx2 * az)
                e0 += fI * ax - da * dx0 + fm * (dx1 * bz - dx2 * by)
                e1 += fI * ay - da * dx1 + fm * (dx2 * bx - dx0 * bz)
                e2 += fI * az - da * dx2 + fm * (dx0 * by - dx1 * bx)
            s0 += sg * e0
            s1 += sg * ry * e1
            s2 += sg * e2
        out[m, 0], out[m, 1], out[m, 2] = s0, s1, s2


@njit(parallel=True, fastmath=True, cache=True)
def _far_midplane(dirs, y, c_mat, comps, sy, k, diag, out):
    """Exact r -> inf far-field pattern for W waves; if diag=True, dirs has shape (W, 3) and out has shape (W, 3)"""
    D, N, W = dirs.shape[0], y.shape[0], c_mat.shape[2]
    c0, c1, c2 = comps[0], comps[1], comps[2]
    for m in prange(D):
        w0 = m if diag else 0
        w1 = m + 1 if diag else W
        for w in range(w0, w1):
            s0, s1, s2 = 0.0j, 0.0j, 0.0j
            for g in range(2):
                ry = 1.0 if g == 0 else -1.0
                sg = 2.0 if g == 0 else 2.0 * sy
                d0, d1, d2 = dirs[m, 0], dirs[m, 1] * ry, dirs[m, 2]
                e0, e1, e2 = 0.0j, 0.0j, 0.0j
                for n in range(N):
                    ang = -k * (d0 * y[n, 0] + d1 * y[n, 1] + d2 * y[n, 2])
                    ph = np.cos(ang) + 1j * np.sin(ang)
                    ax, ay, az, bx, by, bz = 0.0j, 0.0j, 0.0j, 0.0j, 0.0j, 0.0j
                    for idx, cid in ((0, c0), (1, c1), (2, c2)):
                        v = c_mat[n, idx, w]
                        if cid == 0:
                            ax = v
                        elif cid == 1:
                            ay = v
                        elif cid == 2:
                            az = v
                        elif cid == 3:
                            bx = v
                        elif cid == 4:
                            by = v
                        else:
                            bz = v
                    da = d0 * ax + d1 * ay + d2 * az
                    fm = 1j * ph
                    e0 += ph * (ax - da * d0) + fm * (d1 * bz - d2 * by)
                    e1 += ph * (ay - da * d1) + fm * (d2 * bx - d0 * bz)
                    e2 += ph * (az - da * d2) + fm * (d0 * by - d1 * bx)
                s0 += sg * e0
                s1 += sg * ry * e1
                s2 += sg * e2
            if diag:
                out[m, 0, 0], out[m, 1, 0], out[m, 2, 0] = s0, s1, s2
            else:
                out[m, 0, w], out[m, 1, w], out[m, 2, w] = s0, s1, s2


def tail_weights(tail, cls, n_e=2, n_m=1):
    """
    Weight matrix W (3 N_t, n_tail) tying the tail dipoles on each row with (r / r_cut)^e,
    e = nu + 2 (electric) or tau + 2 (magnetic), using the n_e smallest electric and n_m
    smallest magnetic exponents of the class (H is axial, so tau uses opposite parities)
    """
    sy, sz = cls
    comps = midplane_comps(sz)
    pts, row, rr = tail
    if not len(pts):
        return np.zeros((0, 0))
    cols = []
    for sel, exps in ((comps < 3, CONE_NU[(sy, sz)][:n_e]), (comps >= 3, CONE_TAU[(-sy, -sz)][:n_m])):
        idx = np.nonzero(sel)[0]
        for i in np.unique(row):
            p = np.nonzero(row == i)[0]
            for e in exps:
                W = np.zeros((len(pts), 3, len(idx)))
                W[p[:, None], idx, np.arange(len(idx))] = (rr[p] ** (e + 2))[:, None]
                cols.append(W.reshape(3 * len(pts), len(idx)))
    return np.hstack(cols)


def _expand_sources(y, c, cls, tail):
    nf = 3 * len(y)
    c2 = c[:, None] if c.ndim == 1 else c
    if tail is None or not len(tail[0]):
        return np.ascontiguousarray(y, dtype=float), np.ascontiguousarray(c2.reshape(len(y), 3, -1))
    W = tail_weights(tail, cls)
    ct = W @ c2[nf:]
    y_all = np.ascontiguousarray(np.vstack((y, tail[0])), dtype=float)
    c_all = np.ascontiguousarray(np.concatenate((c2[:nf], ct), axis=0).reshape(len(y_all), 3, -1))
    return y_all, c_all


def matrix(x, y, k, T=None, cls=(-1, 1), tail=None, batch=1024):
    """
    Fields of the dipoles at y on the mid-plane (and the tip tail, if given), symmetrised for cls
      x: points (M, 3), T: tangents (M, 2, 3) or None for the Cartesian components
      Returns an F-ordered array (M nc, unknowns)
    """
    x, y = np.ascontiguousarray(x), np.ascontiguousarray(y)
    comps = midplane_comps(cls[1])
    T_ = np.eye(3)[None] if T is None else np.ascontiguousarray(T)
    nc, nf = T_.shape[1], len(y) * len(comps)
    W = np.zeros((0, 0)) if tail is None else tail_weights(tail, cls)
    A = np.empty((len(x) * nc, nf + W.shape[1]), complex, order="F")
    _block(x, T_, y, comps, float(cls[0]), float(k), A[:, :nf])
    if W.shape[1]:
        tpts = np.ascontiguousarray(tail[0])
        buf = np.empty((min(len(x), batch) * nc, 3 * len(tpts)), complex, order="F")
        for i in range(0, len(x), batch):
            j = min(i + batch, len(x))
            gb = buf if (j - i) * nc == buf.shape[0] else np.empty(((j - i) * nc, 3 * len(tpts)), complex, order="F")
            _block(x[i:j], T_[i:j] if T_.shape[0] > 1 else T_, tpts, comps, float(cls[0]), float(k), gb)
            A[i * nc:j * nc, nf:] = gb @ W
    return A


# Sources on the mid-plane

def surf_dist(x, y):
    """Approximate distance from the mid-plane point (x, y, 0) to the surface"""
    a, b = semi_axes(x)
    eta = np.minimum(np.abs(y) / np.maximum(a, 1e-300), 1.0)
    return np.minimum(a - np.abs(y), b * np.sqrt(1 - eta ** 2))


def row_etas(n_mid, n_F):
    """Fractions of the half-width of the rows: n_mid in the middle, n_F + 1 approaching the focal curve"""
    mid = (ETA_F - 3 * GAP) * np.cos((np.arange(n_mid) + 0.5) * np.pi / (2 * n_mid))
    return list(mid) + [ETA_F - 2 * GAP * 0.5 ** j for j in range(n_F + 1)]


def source_rows(n_mid, n_F, s=0.7, tip_cut=0.004, hmax=np.inf, tip_end=4e-7):
    """
    Source points of the layout n_mid:n_F in the quarter y >= 0 of the mid-plane, and the tail at the tip
      s: spacing along a row, as a fraction of the distance to the surface
      Rows follow the outline in the front half and confocal ellipses in the rear half. The free
      points stop tip_cut · L before the tip; the tail continues each row to tip_end · L from it
      Returns the points (N, 3), and the tail: points (N_t, 3), row index (N_t,) and the distance
      to the tip relative to the cut (N_t,)
    """
    x0, xt = L * (T1 - tip_cut), L * T1
    pts, tail = [], []
    for i, eta in enumerate(row_etas(n_mid, n_F)):
        # Front half: y = eta a(x), up to the cut, then the tail
        x = 0.5 * min(s * surf_dist(0.0, eta * semi_axes(0.0)[0]), hmax)
        while x < L * (T1 - tip_end):
            y = eta * semi_axes(x)[0]
            if x < x0:
                pts.append((x, y, 0.0))
            else:
                tail.append((x, y, 0.0, i, (xt - x) / (xt - x0)))
            x += min(s * surf_dist(x, y), hmax)

        # Rear half: the ellipse through the same point of the widest section, confocal with the focal ellipse
        Amu = np.sqrt(A_ ** 2 - B_ ** 2 * (1 - eta ** 2))
        psi = 0.5 * min(s * surf_dist(0.0, eta * B_), hmax) / Amu
        while psi < np.pi / 2:
            x, y = -Amu * np.sin(psi), eta * B_ * np.cos(psi)
            ds = min(s * surf_dist(x, y), hmax)
            if y > 0.3 * ds:
                pts.append((x, y, 0.0))
            psi += ds / np.sqrt((Amu * np.cos(psi)) ** 2 + (eta * B_ * np.sin(psi)) ** 2)
    tail = np.array(tail, float).reshape(-1, 5)
    return np.array(pts), (tail[:, :3], tail[:, 3].astype(int), tail[:, 4])


# Boundary points on the quarter y, z >= 0 of the surface

def surface_point(u, phi, rear):
    """Front half: u = x. Rear half: x = -A sin(u)"""
    if rear:
        x, a, b = -A_ * np.sin(u), B_ * np.cos(u), C_ * np.cos(u)
    else:
        x, (a, b) = u, semi_axes(u)
    return np.stack((x * np.ones_like(phi), a * np.cos(phi), b * np.sin(phi)), axis=-1)


def _count_collocation(s, hmax, tip_cut=0.002, end_cut=0.02):
    """Fast scalar count of collocation points without forming surface arrays or finite-difference normals"""
    rings, phi = [], 0.0
    while phi < math.pi / 2:
        sp, cp = math.sin(phi), math.cos(phi)
        sp9, cp9 = math.sin(phi + 1e-9), math.cos(phi + 1e-9)
        rd = math.sqrt((sp9 / 3.0) ** 2 + max(0.0, cp9 - ETA_F) ** 2)
        dphi = min(s * rd / math.sqrt(sp * sp + cp * cp / 9.0), hmax / B_, math.pi / 2 - phi)
        rings.append((phi + 0.5 * dphi, dphi))
        phi += dphi
    cnt = 0
    for phi, _ in rings:
        sp, cp = math.sin(phi), math.cos(phi)
        rd = math.sqrt((sp / 3.0) ** 2 + max(0.0, cp - ETA_F) ** 2)
        cf = cp * cp + sp * sp / 9.0
        for rear in (False, True):
            u, u_end = 0.0, ((math.pi / 2 - end_cut) if rear else L * (T1 - tip_cut))
            while u < u_end:
                if rear:
                    cu, su = math.cos(u), math.sin(u)
                    a = B_ * cu
                    metric = math.sqrt((A_ * cu) ** 2 + (B_ * su) ** 2 * cf)
                else:
                    t = u / L
                    a = 4.83345 * (math.sqrt(max(0.0, 1.0 - (t / 2.08335) ** 2)) - 0.96) * L
                    metric = 1.0
                du = min(min(s * a * rd, hmax) / metric, u_end - u)
                cnt += 1
                u += du
    return cnt


def collocation(s=0.35, hmax=0.012, tip_cut=0.002, end_cut=0.02):
    """
    Boundary points that refine towards the rim, as the sources do
      Returns points (M, 3), area weights (M,) and outward normals (M, 3)
    """
    key = (round(float(s), 12), round(float(hmax), 12), tip_cut, end_cut)
    if key in _COLLOC_CACHE:
        return _COLLOC_CACHE[key]

    def rel_dist(phi):                                          # distance to the focal curve, relative to the half-width
        return math.sqrt((math.sin(phi) / 3.0) ** 2 + max(0.0, math.cos(phi) - ETA_F) ** 2)

    rings, phi = [], 0.0
    while phi < math.pi / 2:
        sp, cp = math.sin(phi), math.cos(phi)
        dphi = min(s * rel_dist(phi + 1e-9) / math.sqrt(sp * sp + cp * cp / 9.0), hmax / B_, math.pi / 2 - phi)
        rings.append((phi + 0.5 * dphi, dphi))
        phi += dphi

    pts, wts, nrm, e = [], [], [], 1e-7
    for phi, dphi in rings:
        rd, cf = rel_dist(phi), math.cos(phi) ** 2 + math.sin(phi) ** 2 / 9.0
        for rear in (False, True):
            u, u_end, ums, dus = 0.0, ((math.pi / 2 - end_cut) if rear else L * (T1 - tip_cut)), [], []
            while u < u_end:
                if rear:
                    cu, su = math.cos(u), math.sin(u)
                    a = B_ * cu
                    metric = math.sqrt((A_ * cu) ** 2 + (B_ * su) ** 2 * cf)
                else:
                    t = u / L
                    a = max(0.0, 4.83345 * (math.sqrt(max(0.0, 1.0 - (t / 2.08335) ** 2)) - 0.96) * L)
                    metric = 1.0
                du = min(min(s * a * rd, hmax) / metric, u_end - u)
                ums.append(u + 0.5 * du)
                dus.append(du)
                u += du
            um, du = np.array(ums), np.array(dus)
            ph = np.full_like(um, phi)
            pu = (surface_point(um + e, ph, rear) - surface_point(um - e, ph, rear)) / (2 * e)
            pp = (surface_point(um, ph + e, rear) - surface_point(um, ph - e, rear)) / (2 * e)
            n = np.cross(pu, pp)
            ln = np.linalg.norm(n, axis=1)
            pts.append(surface_point(um, ph, rear))
            wts.append(ln * du * dphi)
            nrm.append((n if rear else -n) / ln[:, None])       # u runs towards -x on the rear and +x on the front
    out = (np.vstack(pts), np.concatenate(wts), np.vstack(nrm))
    if len(_COLLOC_CACHE) >= 8:
        _COLLOC_CACHE.clear()
    _COLLOC_CACHE[key] = out
    return out


def solve_class(k, y, cls, d, p, tail=None, s_col=0.35, ratio=2.5):
    """
    Least-squares MFS for one parity class and any number of waves (d, p), shapes (W, 3)
      Returns the coefficients (unknowns, W), the relative boundary residual of each wave, and the sizes
    """
    hmax = min(0.012, 2 * np.pi / k / 5)
    m0 = _count_collocation(s_col, hmax)
    if 2 * m0 < ratio * 3 * len(y):
        s_col *= np.sqrt(2 * m0 / (ratio * 3 * len(y)))
    x, w, nrm = collocation(s_col, hmax)
    T = tangent_basis(nrm)
    sw = np.sqrt(w / w.sum())

    # Form F-ordered matrix A and the right-hand sides, rows weighted by sqrt(area)
    Tw = T * sw[:, None, None]
    A = matrix(x, y, k, Tw, cls, tail)
    b = -np.einsum("mab,mwb->maw", Tw, class_part(x, k, d, p, REFL, (1, *cls))).reshape(-1, len(d))

    # Solve for all waves at once
    c, res = factor(A, overwrite=True).solve(b, residual=True)
    return c, res / np.maximum(np.linalg.norm(b, axis=0), 1e-300), A.shape


def far(k, y, c, dirs, cls, tail=None):
    """Exact far-field pattern F (E ~ F e^{ikr}/r) in the directions dirs, for the coefficients c of one class"""
    dirs = np.ascontiguousarray(dirs, dtype=float)
    y_all, c_all = _expand_sources(y, c, cls, tail)
    out = np.empty((len(dirs), 3, c_all.shape[2]), complex)
    _far_midplane(dirs, y_all, c_all, midplane_comps(cls[1]), float(cls[0]), float(k), False, out)
    return out[:, :, 0] if c.ndim == 1 else out


def nose_on(k, n_mid, n_F, pts=None, s=0.7, s_col=0.35, hmax=np.inf):
    """
    The wave along -x (meeting the tip first), polarised along y, which excites the single
    parity class (-1, +1)
      Returns (y, tail), the coefficients, the field at pts and the shape of the matrix
    """
    d, p = np.array([[-1.0, 0, 0]]), np.array([[0.0, 1, 0]])
    y, tail = source_rows(n_mid, n_F, s=s, hmax=hmax)
    c, _, shape = solve_class(k, y, (-1, 1), d, p, tail, s_col=s_col)
    if pts is None:
        E = None
    else:
        y_all, c_all = _expand_sources(y, c[:, 0], (-1, 1), tail)
        E = np.empty((len(pts), 3), complex)
        _eval_midplane(np.ascontiguousarray(pts, dtype=float), y_all, c_all[:, :, 0],
                       midplane_comps(1), -1.0, float(k), E)
    return (y, tail), c[:, 0], E, shape


def optical(k, src, c):
    """Extinction and scattering cross-sections (ext, sca) of the nose-on solution: equal for the exact field"""
    y, tail = src if isinstance(src, tuple) else (src, None)
    d, p = np.array([-1.0, 0, 0]), np.array([0.0, 1, 0])
    return cross_sections(lambda dirs: far(k, y, c, dirs, (-1, 1), tail), k, d, p, sphere_rule(48, 96))


def rcs(freq, n_mid=6, n_F=3, step=0.5, s=0.5, s_col=0.3):
    """
    Monostatic radar cross-section against the azimuth in the plane of the almond (0 = nose-on)
      freq: frequency [GHz]; s: spacing of the sources along a row
      HH (electric field in the plane) needs the two classes even in z, VV the two odd in z;
      the factorisation of a class serves all look angles
      Returns the angles [deg] and sigma_HH, sigma_VV [m²]
    """
    k = 2 * np.pi * freq * 1e9 / C0
    hmax = 2 * np.pi / k / 4
    y, tail = source_rows(n_mid, n_F, s=s, hmax=hmax)
    phi = np.radians(np.arange(0.0, 180.0 + 1e-9, step))
    d = -np.column_stack((np.cos(phi), np.sin(phi), 0 * phi))
    back = np.ascontiguousarray(-d)
    out = {"phi_deg": np.degrees(phi)}
    p_hh = np.column_stack((-np.sin(phi), np.cos(phi), 0 * phi))     # electric field in the plane of the almond
    p_vv = np.tile([0.0, 0, 1], (len(phi), 1))                        # electric field along z
    F_diag = np.empty((len(phi), 3, 1), complex)
    for pol, p, sz in (("HH", p_hh, 1), ("VV", p_vv, -1)):
        amp = 0
        for sy in (1, -1):
            c, _, _ = solve_class(k, y, (sy, sz), d, p, tail, s_col)
            y_all, c_all = _expand_sources(y, c, (sy, sz), tail)
            _far_midplane(back, y_all, c_all, midplane_comps(sz), float(sy), float(k), True, F_diag)
            amp = amp + np.einsum("na,na->n", p, F_diag[:, :, 0])
        out["sigma_" + pol] = 4 * np.pi * np.abs(amp) ** 2
    return out
