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
CENTRE = np.array([0.02, 0, 0])                             # a point inside the almond, to orient the normals

# Exponents of the tip cone (half-angles 34.09 and 12.71 degrees) by the parities (y, z) of the
# potential: electric type E = ∇(r^nu Y), magnetic type H = ∇(r^tau Z) (mfs.cone, Table A.1)
CONE_NU = {(1, 1): (0.29959, 1.46427), (-1, 1): (1.10372, 2.24730), (1, -1): (1.04120, 2.10576), (-1, -1): (2.01295, 3.04849)}
CONE_TAU = {(1, 1): (1.05604, 1.98580), (-1, 1): (0.96803, 1.96545), (1, -1): (0.87905, 1.75454), (-1, -1): (1.97538, 2.90094)}


def semi_axes(x):
    """Semi-axes (a, b) of the cross-section at x (metres)"""
    t = np.asarray(x, float) / L
    f1 = np.sqrt(np.clip(1 - (t / 0.416667) ** 2, 0, None))
    f2 = np.sqrt(np.clip(1 - (t / 2.08335) ** 2, 0, None)) - 0.96
    a = np.where(t < 0, 0.193333 * f1, 4.83345 * f2) * L
    b = np.where(t < 0, 0.064444 * f1, 1.61115 * f2) * L
    return np.clip(a, 0, None), np.clip(b, 0, None)


def parity_class(sy, sz):
    """
    Dipole components that survive on the mid-plane for the class (sy, sz), and the signs
    of the four images. Components 0-2: electric x, y, z; 3-5: magnetic x, y, z
    """
    return (np.array([0, 1, 5]) if sz > 0 else np.array([2, 3, 4])), np.array([1.0, sy, sz, sy * sz])


@njit(parallel=True, fastmath=True, cache=True)
def _block(x, T, use_T, y, comps, sig, k, out):
    """out[m, a, n, i]: component a (or tangent a) at x[m] of dipole comps[i] at y[n], summed over the images"""
    for m in prange(x.shape[0]):
        v = np.empty(3, np.complex128)
        out[m] = 0.0
        for g in range(4):
            r0, r1, r2 = x[m, 0], x[m, 1] * (1 - 2 * (g % 2)), x[m, 2] * (1 - 2 * (g // 2))
            for n in range(y.shape[0]):
                d0, d1, d2 = r0 - y[n, 0], r1 - y[n, 1], r2 - y[n, 2]
                r = np.sqrt(d0 * d0 + d1 * d1 + d2 * d2)
                kr = k * r
                ph = (np.cos(kr) + 1j * np.sin(kr)) / r
                e0, e1, e2 = d0 / r, d1 / r, d2 / r
                fI = (1.0 + 1j / kr - 1.0 / (kr * kr)) * ph
                fR = (1.0 + 3j / kr - 3.0 / (kr * kr)) * ph
                fm = (1j - 1.0 / kr) * ph
                for i in range(comps.shape[0]):
                    c = comps[i]
                    if c < 3:                                   # electric dipole along e_c
                        ec = e0 if c == 0 else (e1 if c == 1 else e2)
                        v[0], v[1], v[2] = -fR * e0 * ec, -fR * e1 * ec, -fR * e2 * ec
                        v[c] += fI
                    elif c == 3:                                # magnetic dipoles: fm e × e_c
                        v[0], v[1], v[2] = 0.0, fm * e2, -fm * e1
                    elif c == 4:
                        v[0], v[1], v[2] = -fm * e2, 0.0, fm * e0
                    else:
                        v[0], v[1], v[2] = fm * e1, -fm * e0, 0.0
                    s = sig[g]
                    v1, v2 = v[1] * (1 - 2 * (g % 2)), v[2] * (1 - 2 * (g // 2))
                    if use_T:
                        for a in range(2):
                            out[m, a, n, i] += s * (T[m, a, 0] * v[0] + T[m, a, 1] * v1 + T[m, a, 2] * v2)
                    else:
                        out[m, 0, n, i] += s * v[0]
                        out[m, 1, n, i] += s * v1
                        out[m, 2, n, i] += s * v2


def matrix(x, y, k, T=None, cls=(-1, 1)):
    """
    Fields of the dipoles at y on the mid-plane, symmetrised for the class cls
      x: points (M, 3), T: tangents (M, 2, 3) or None for the Cartesian components
      Returns an array (M nc, N · 3)
    """
    comps, sig = parity_class(*cls)
    nc = 3 if T is None else 2
    out = np.empty((len(x), nc, len(y), len(comps)), complex)
    _block(np.ascontiguousarray(x), np.zeros((1, 2, 3)) if T is None else np.ascontiguousarray(T),
           T is not None, np.ascontiguousarray(y), comps, sig, k, out)
    return out.reshape(len(x) * nc, -1)


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
        Amu, psi, first = np.sqrt(A_ ** 2 - B_ ** 2 * (1 - eta ** 2)), 0.0, True
        while psi < np.pi / 2:
            x, y = -Amu * np.sin(psi), eta * B_ * np.cos(psi)
            ds = min(s * surf_dist(x, y), hmax)
            dpsi = ds / np.sqrt((Amu * np.cos(psi)) ** 2 + (eta * B_ * np.sin(psi)) ** 2)
            if not first and y > 0.3 * ds:
                pts.append((x, y, 0.0))
            psi, first = psi + (0.5 * dpsi if first else dpsi), False
    tail = np.array(tail)
    return np.array(pts), (tail[:, :3], tail[:, 3].astype(int), tail[:, 4])


def tail_matrix(x, tail, k, T, cls, n_e=2, n_m=1):
    """
    Columns of the tip tail: on each row the dipoles are tied with the weights (r / r_cut)^e,
    e = nu + 2 (electric) or tau + 2 (magnetic), with the n_e smallest electric and the
    n_m smallest magnetic exponents of the class
    """
    sy, sz = cls
    comps, _ = parity_class(sy, sz)
    pts, row, rr = tail
    nc = 3 if T is None else 2
    G = matrix(x, pts, k, T, cls).reshape(len(x) * nc, len(pts), 3)
    cols = []
    # H is a pseudovector, so the potential of the magnetic type has the opposite parities
    for sel, exps in ((comps < 3, CONE_NU[(sy, sz)][:n_e]), (comps >= 3, CONE_TAU[(-sy, -sz)][:n_m])):
        cols += [np.einsum("mnc,n->mc", G[:, row == i][:, :, sel], rr[row == i] ** (e + 2)) for i in np.unique(row) for e in exps]
    return np.concatenate(cols, axis=1)


# Boundary points on the quarter y, z >= 0 of the surface

def surface_point(u, phi, rear):
    """Front half: u = x. Rear half: x = -A sin(u)"""
    if rear:
        x, a, b = -A_ * np.sin(u), B_ * np.cos(u), C_ * np.cos(u)
    else:
        x, (a, b) = u, semi_axes(u)
    return np.stack((x * np.ones_like(phi), a * np.cos(phi), b * np.sin(phi)), axis=-1)


def collocation(s=0.35, hmax=0.012, tip_cut=0.002, end_cut=0.02):
    """
    Boundary points that refine towards the rim, as the sources do
      Returns points (M, 3), area weights (M,) and outward normals (M, 3)
    """
    def rel_dist(phi):                                          # distance to the focal curve, relative to the half-width
        return np.sqrt((np.sin(phi) / 3) ** 2 + np.maximum(0.0, np.cos(phi) - ETA_F) ** 2)

    # Angles around the cross-section, the same for every x (every section is a 3:1 ellipse)
    phis, dphis, phi = [], [], 0.0
    while phi < np.pi / 2:
        dphi = min(s * rel_dist(phi + 1e-9) / np.sqrt(np.sin(phi) ** 2 + np.cos(phi) ** 2 / 9), hmax / B_, np.pi / 2 - phi)
        phis.append(phi + 0.5 * dphi)
        dphis.append(dphi)
        phi += dphi

    # Along the body, in the front and in the rear half, with normals by finite differences
    pts, wts, nrm = [], [], []
    e = 1e-7                                                    # step of the finite differences
    for phi, dphi in zip(phis, dphis):
        for rear in (False, True):
            u, u_end = 0.0, (np.pi / 2 - end_cut) if rear else L * (T1 - tip_cut)
            while u < u_end:
                if rear:
                    a = B_ * np.cos(u)
                    metric = np.sqrt((A_ * np.cos(u)) ** 2 + (B_ * np.sin(u)) ** 2 * (np.cos(phi) ** 2 + np.sin(phi) ** 2 / 9))
                else:
                    a, metric = semi_axes(u)[0], 1.0
                du = min(min(s * a * rel_dist(phi), hmax) / metric, u_end - u)
                um = np.array(u + 0.5 * du)
                pu = (surface_point(um + e, np.array(phi), rear) - surface_point(um - e, np.array(phi), rear)) / (2 * e)
                pp = (surface_point(um, np.array(phi + e), rear) - surface_point(um, np.array(phi - e), rear)) / (2 * e)
                n = np.cross(pu, pp)
                pts.append(surface_point(um, np.array(phi), rear))
                wts.append(np.linalg.norm(n) * du * dphi)
                nrm.append(n / np.linalg.norm(n))
                u += du
    pts, nrm = np.array(pts), np.array(nrm)
    nrm *= np.sign(np.sum(nrm * (pts - CENTRE), axis=1))[:, None]
    return pts, np.array(wts), nrm


def solve_class(k, y, cls, d, p, tail=None, s_col=0.35, ratio=2.5):
    """
    Least-squares MFS for one parity class and any number of waves (d, p), shapes (W, 3)
      Returns the coefficients (unknowns, W), the relative boundary residual of each wave, and the sizes
    """
    # Boundary points, with at least ratio equations per unknown
    hmax = min(0.012, 2 * np.pi / k / 5)
    x, w, nrm = collocation(s_col, hmax)
    if 2 * len(x) < ratio * 3 * len(y):
        s_col *= np.sqrt(2 * len(x) / (ratio * 3 * len(y)))
        x, w, nrm = collocation(s_col, hmax)
    T = tangent_basis(nrm)
    sw = np.repeat(np.sqrt(w / w.sum()), 2)[:, None]

    # Form matrix A and the right-hand sides, rows weighted by sqrt(area)
    A = matrix(x, y, k, T, cls)
    if tail is not None:
        A = np.concatenate((A, tail_matrix(x, tail, k, T, cls)), axis=1)
    A *= sw
    b = -np.einsum("mab,mwb->maw", T, class_part(x, k, d, p, REFL, (1, *cls))).reshape(-1, len(d)) * sw

    # Solve for all waves at once
    c, res = factor(np.ascontiguousarray(A)).solve(b, residual=True)
    return c, res / np.maximum(np.linalg.norm(b, axis=0), 1e-300), A.shape


def far(k, y, c, dirs, cls, tail=None, rfar=1.0e5):
    """Far-field pattern F (E ~ F e^{ikr}/r) in the directions dirs, for the coefficients c of one class"""
    A = matrix(rfar * dirs, y, k, None, cls)
    if tail is not None:
        A = np.concatenate((A, tail_matrix(rfar * dirs, tail, k, None, cls)), axis=1)
    return (A @ c).reshape(len(dirs), 3, *c.shape[1:]) * rfar * np.exp(-1j * k * rfar)


def nose_on(k, n_mid, n_F, pts=None):
    """
    The wave along -x (meeting the tip first), polarised along y. It has the parities
    (-1, +1), which carry no current into the tip, so no tail is needed.
      Returns the source points, the coefficients, the field at pts, the extinction and
      scattering cross-sections (ext, sca), and the shape of the matrix
    """
    d, p = np.array([[-1.0, 0, 0]]), np.array([[0.0, 1, 0]])
    y, _ = source_rows(n_mid, n_F)
    c, _, shape = solve_class(k, y, (-1, 1), d, p)
    E = None if pts is None else (matrix(pts, y, k) @ c[:, 0]).reshape(-1, 3)
    ext, sca = cross_sections(lambda dirs: far(k, y, c[:, 0], dirs, (-1, 1)), k, d[0], p[0], sphere_rule(48, 96))
    return y, c[:, 0], E, (ext, sca), shape


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
    out = {"phi_deg": np.degrees(phi)}
    p_hh = np.column_stack((-np.sin(phi), np.cos(phi), 0 * phi))     # electric field in the plane of the almond
    p_vv = np.tile([0.0, 0, 1], (len(phi), 1))                        # electric field along z
    for pol, p, sz in (("HH", p_hh, 1), ("VV", p_vv, -1)):
        amp = 0
        for sy in (1, -1):
            c, _, _ = solve_class(k, y, (sy, sz), d, p, tail, s_col)
            F = far(k, y, c, -d, (sy, sz), tail)                   # (waves, 3, waves): backscatter of wave n is F[n, :, n]
            amp = amp + np.einsum("na,nan->n", p, F)
        out["sigma_" + pol] = 4 * np.pi * np.abs(amp) ** 2
    return out
