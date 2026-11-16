"""MFS for a perfectly conducting NASA almond with "skeleton" sources.

Every cross-section of the almond is a 3:1 ellipse, so the scattered field
continues analytically inside the body up to (roughly) the focal curve in the
mid-plane z = 0, which lies within 6% of the rim. Sources are electric and
magnetic point dipoles on the mid-plane, in rows that follow the rim and
cluster towards the focal curve. Their spacing along a row is proportional to
the local distance to the surface, as is the spacing of the collocation points.

This module is the library of the almond solvers: geometry, source rows,
collocation, the symmetrised matrices for any mirror-parity class
(`parity_class`), and the tip tail (`tip_tail`, `tail_matrix`, with exponents
from the tip cone, `tip_exponents`). The general solver, for any look angle
and both polarisations, is almond_rcs.py.

Run as a script it solves one case only: nose-on incidence (travelling in -x,
towards the tip), polarised along y. That field is odd under y -> -y and even
under z -> -z, a class that does not carry current into the tip, so the
script uses neither the tip tail nor the finer row spacing (0.5) that the
class even under both mirrors needs; do not copy its settings for other
incidences.

Usage: python pec_almond_skeleton.py [freq_GHz]
"""
import json
import sys
import time

import numpy as np
from numba import njit, prange
from scipy.linalg import lstsq, qr, solve_triangular

from pec_almond import L, T0, T1, semi_axes
from pec_sphere_convergence import fibonacci_sphere, tangent_basis

A_, B_, C_ = 0.416667 * L, 0.193333 * L, 0.064444 * L
ETA_F = np.sqrt(1 - (C_ / B_) ** 2)          # focal curve, as a fraction of the half-width
GAP = 1 - ETA_F
D_INC = np.array([-1.0, 0.0, 0.0])
P_INC = np.array([0.0, 1.0, 0.0])
REFL = np.array([[1.0, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1]])
SIG = np.array([1.0, -1.0, 1.0, -1.0])
COMPS = np.array([0, 1, 5])                  # e_x, e_y, m_z survive on the mid-plane


def ab(x):
    return semi_axes(np.asarray(x, float) / L)


@njit(parallel=True, fastmath=True, cache=True)
def sym_block(x, T, use_T, y, comps, refl, sig, k, out):
    """Symmetrised fields of point dipoles. comps index 0-2: electric x,y,z;
    3-5: magnetic x,y,z. out: (M, nc, N, len(comps))."""
    M, N, ncmp = x.shape[0], y.shape[0], comps.shape[0]
    nc = 2 if use_T else 3
    for m in prange(M):
        v = np.empty(3, np.complex128)
        for n in range(N):
            for ci in range(ncmp):
                for a in range(nc):
                    out[m, a, n, ci] = 0.0
        for gi in range(refl.shape[0]):
            xr0, xr1, xr2 = x[m, 0] * refl[gi, 0], x[m, 1] * refl[gi, 1], x[m, 2] * refl[gi, 2]
            for n in range(N):
                d0, d1, d2 = xr0 - y[n, 0], xr1 - y[n, 1], xr2 - y[n, 2]
                r = np.sqrt(d0 * d0 + d1 * d1 + d2 * d2)
                kr = k * r
                ph = (np.cos(kr) + 1j * np.sin(kr)) / r
                e0, e1, e2 = d0 / r, d1 / r, d2 / r
                fI = (1.0 + 1j / kr - 1.0 / (kr * kr)) * ph
                fR = (1.0 + 3j / kr - 3.0 / (kr * kr)) * ph
                fm = (1j - 1.0 / kr) * ph
                for ci in range(ncmp):
                    c = comps[ci]
                    if c == 0:
                        v[0] = fI - fR * e0 * e0; v[1] = -fR * e1 * e0; v[2] = -fR * e2 * e0
                    elif c == 1:
                        v[0] = -fR * e0 * e1; v[1] = fI - fR * e1 * e1; v[2] = -fR * e2 * e1
                    elif c == 2:
                        v[0] = -fR * e0 * e2; v[1] = -fR * e1 * e2; v[2] = fI - fR * e2 * e2
                    elif c == 3:
                        v[0] = 0.0; v[1] = fm * e2; v[2] = -fm * e1
                    elif c == 4:
                        v[0] = -fm * e2; v[1] = 0.0; v[2] = fm * e0
                    else:
                        v[0] = fm * e1; v[1] = -fm * e0; v[2] = 0.0
                    v0, v1, v2 = v[0] * refl[gi, 0], v[1] * refl[gi, 1], v[2] * refl[gi, 2]
                    if use_T:
                        for a in range(2):
                            out[m, a, n, ci] += sig[gi] * (T[m, a, 0] * v0 + T[m, a, 1] * v1 + T[m, a, 2] * v2)
                    else:
                        out[m, 0, n, ci] += sig[gi] * v0
                        out[m, 1, n, ci] += sig[gi] * v1
                        out[m, 2, n, ci] += sig[gi] * v2


def parity_class(sy, sz):
    """Dipole components that survive on the mid-plane, and the image signs,
    for fields with parity sy under y -> -y and sz under z -> -z (parity of
    the tangential-to-mirror electric field). The nose-on case is (-1, +1)."""
    comps = np.array([0, 1, 5]) if sz > 0 else np.array([2, 3, 4])
    return comps, np.array([1.0, sy, sz, sy * sz])


def matrix(x, y, k, T=None, comps=COMPS, sig=SIG):
    nc = 2 if T is not None else 3
    out = np.empty((len(x), nc, len(y), len(comps)), complex)
    sym_block(np.ascontiguousarray(x), np.ascontiguousarray(T) if T is not None else np.zeros((1, 2, 3)),
              T is not None, np.ascontiguousarray(y), comps, REFL, sig, k, out)
    return out.reshape(len(x) * nc, len(y) * len(comps))


# ---------- sources on the mid-plane -------------------------------------------

def surf_dist(x, y):
    """Approximate distance from the mid-plane point (x, y, 0) to the surface."""
    a, b = ab(x)
    eta = np.minimum(np.abs(y) / np.maximum(a, 1e-300), 1.0)
    return np.minimum(a - np.abs(y), b * np.sqrt(1 - eta ** 2))


def row_etas(n_mid, J, q, beyond=()):
    etas = list((ETA_F - 3 * GAP) * np.cos((np.arange(n_mid) + 0.5) * np.pi / (2 * n_mid)))
    etas += [ETA_F - 2 * GAP * q ** j for j in range(J + 1)]
    return etas + [ETA_F + f * GAP for f in beyond]


TIP_EXPS = (2.2, 2.4, 2.7, 3.0)          # the fitted set used before the cone exponents were computed

# Exponents of the tip cone (half-angles 34.09 and 12.71 degrees), from cone_exponents.py:
# electric type E = grad(r^nu Y), magnetic type H = grad(r^tau Z); keyed by the mirror
# parities (y, z) of the potential Y or Z.
CONE_NU = {(1, 1): (0.29959, 1.46427, 2.02232), (-1, 1): (1.10372, 2.24730), (1, -1): (1.04120, 2.10576),
           (-1, -1): (2.01295, 3.04849)}
CONE_TAU = {(1, 1): (1.05604, 1.98580), (-1, 1): (0.96803, 1.96545), (1, -1): (0.87905, 1.75454),
            (-1, -1): (1.97538, 2.90094)}


def tip_exponents(sy, sz, n_e=2, n_m=1):
    """Dipole exponents of the tip tail for the class (sy, sz), from the cone.
    A field E ~ r^(nu - 1) is produced by electric dipoles of strength r^(nu + 2)
    (per source, sources spaced in proportion to r); its potential has the
    parities of the class. The magnetic-type fields H ~ r^(tau - 1) come from
    magnetic dipoles r^(tau + 2); H is a pseudovector, so its potential has the
    opposite parities. Returns (electric exponents, magnetic exponents)."""
    return (tuple(v + 2 for v in CONE_NU[(sy, sz)][:n_e]), tuple(v + 2 for v in CONE_TAU[(-sy, -sz)][:n_m]))


def tip_tail(n_mid, J, q, beyond=(), s=0.7, tip_cut=0.004, hmax=np.inf, tip_end=4e-7):
    """Continuation of every source row from where it stops (tip_cut * L from
    the tip) towards the tip, by the same spacing rule. Near the tip the body
    is an elliptic cone and the field behaves like r^(nu - 1) with an exponent
    nu that is small when the wave drives current along the axis; free sources
    cut off at a fixed distance then leave an error that decays only like
    cut^nu. The continued points are not free: on each row they are tied
    together with weights (r / r_cut)^e, dipole exponent e = nu + 2, for a
    small set of exponents (TIP_EXPS) among which the fit chooses.
    Returns points (N, 3), row index (N,), r / r_cut (N,)."""
    x0, xt = L * (T1 - tip_cut), L * T1
    pts, row, rr = [], [], []
    for i, eta in enumerate(row_etas(n_mid, J, q, beyond)):
        x = 0.5 * min(s * surf_dist(0.0, eta * ab(0.0)[0]), hmax)
        while x < x0:
            x += min(s * surf_dist(x, eta * ab(x)[0]), hmax)
        while x < L * (T1 - tip_end):
            yv = eta * ab(x)[0]
            pts.append((x, yv, 0.0))
            row.append(i)
            rr.append((xt - x) / (xt - x0))
            x += s * surf_dist(x, yv)
    return np.array(pts), np.array(row), np.array(rr)


def tail_matrix(x, tail, k, T=None, comps=COMPS, sig=SIG, exps=TIP_EXPS, mexps=()):
    """Columns of the tip tail: per row, per electric component of the class,
    per exponent (and per magnetic component and exponent in mexps)."""
    pts, row, rr = tail
    nc = 2 if T is not None else 3
    cols = []
    for cc, ee in ((comps[comps < 3], exps), (comps[comps >= 3], mexps)):
        if len(ee) == 0 or len(cc) == 0:
            continue
        M = matrix(x, pts, k, T, cc, sig).reshape(len(x) * nc, len(pts), len(cc))
        cols += [np.einsum("mnc,n->mc", M[:, row == i, :], rr[row == i] ** e) for i in np.unique(row) for e in ee]
    return np.concatenate(cols, axis=1) if cols else np.zeros((len(x) * nc, 0), complex)


def source_rows(n_mid, J, q, beyond=(), s=0.7, tip_cut=0.004, hmax=np.inf):
    etas = row_etas(n_mid, J, q, beyond)
    pts = []
    x_end = L * (T1 - tip_cut)
    for eta in etas:
        # front part: y = eta * a(x)
        x = 0.0
        step = min(s * surf_dist(0.0, eta * ab(0.0)[0]), hmax)
        x += 0.5 * step
        while x < x_end:
            yv = eta * ab(x)[0]
            pts.append((x, yv))
            x += min(s * surf_dist(x, yv), hmax)
        # rear part: confocal ellipse in the mid-plane
        Amu = np.sqrt(A_ ** 2 - B_ ** 2 * (1 - eta ** 2))
        psi = 0.0
        first = True
        while psi < np.pi / 2:
            xv, yv = -Amu * np.sin(psi), eta * B_ * np.cos(psi)
            dist = surf_dist(xv, yv)
            ds = min(s * dist, hmax)
            dpsi = ds / np.sqrt((Amu * np.cos(psi)) ** 2 + (eta * B_ * np.sin(psi)) ** 2)
            if first:
                psi += 0.5 * dpsi
                first = False
                continue
            if yv > 0.3 * ds:
                pts.append((xv, yv))
            psi += dpsi
    pts = np.array(pts)
    return np.column_stack((pts, np.zeros(len(pts))))


# ---------- collocation points on the surface ----------------------------------

def surface_pt(u, phi, rear):
    """Surface parametrisation. Front: u = x. Rear: u = psi, x = -A sin(psi)."""
    if rear:
        x = -A_ * np.sin(u)
        a, b = B_ * np.cos(u), C_ * np.cos(u)
    else:
        x = u
        a, b = ab(u)
    return np.stack((x * np.ones_like(phi), a * np.cos(phi), b * np.sin(phi)), axis=-1)


def collocation(eta_max, s=0.35, hmax=0.012, tip_cut=0.002, end_cut=0.02):
    # angular sequence (independent of x because every section is a 3:1 ellipse)
    def rel_dist(phi):
        return np.sqrt((np.sin(phi) / 3) ** 2 + np.maximum(0.0, np.cos(phi) - eta_max) ** 2)
    phis, dphis = [], []
    phi = 0.0
    while phi < np.pi / 2:
        dphi = min(s * rel_dist(phi + 1e-9) / np.sqrt(np.sin(phi) ** 2 + np.cos(phi) ** 2 / 9), hmax / B_)
        dphi = min(dphi, np.pi / 2 - phi)
        phis.append(phi + 0.5 * dphi)
        dphis.append(dphi)
        phi += dphi
    pts, wts, nrm = [], [], []
    e = 1e-7
    for phi, dphi in zip(phis, dphis):
        for rear in (False, True):
            u, u_end = 0.0, (np.pi / 2 - end_cut) if rear else L * (T1 - tip_cut)
            while u < u_end:
                if rear:
                    a = B_ * np.cos(u)
                    metric = np.sqrt((A_ * np.cos(u)) ** 2 + (B_ * np.sin(u)) ** 2 * (np.cos(phi) ** 2 + np.sin(phi) ** 2 / 9))
                else:
                    a = ab(u)[0]
                    metric = 1.0
                du = min(s * a * rel_dist(phi), hmax) / metric
                du = min(du, u_end - u)
                um = u + 0.5 * du
                p0 = surface_pt(np.array(um), np.array(phi), rear)
                pu = (surface_pt(np.array(um + e), np.array(phi), rear) - surface_pt(np.array(um - e), np.array(phi), rear)) / (2 * e)
                pp = (surface_pt(np.array(um), np.array(phi + e), rear) - surface_pt(np.array(um), np.array(phi - e), rear)) / (2 * e)
                n = np.cross(pu, pp)
                jac = np.linalg.norm(n)
                pts.append(p0)
                wts.append(jac * du * dphi)
                nrm.append(n / jac)
                u += du
    pts, nrm = np.array(pts), np.array(nrm)
    nrm *= np.sign(np.sum(nrm * (pts - np.array([0.02, 0, 0])), axis=1))[:, None]
    return pts, np.array(wts), nrm


def incident(x, k, d=D_INC, p=P_INC):
    return p[None, :] * np.exp(1j * k * (x @ d))[:, None]


def cross_sections(k, y, c, nt=48, nphi=96, rfar=1.0e5):
    """Extinction and scattering cross-sections of the computed field; equal
    for a perfect conductor (optical theorem). Their mismatch is an error
    measure that needs no reference solution."""
    g, w = np.polynomial.legendre.leggauss(nt)
    ph = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    ct, phh = np.meshgrid(g, ph, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    dirs = np.column_stack(((st * np.cos(phh)).ravel(), (st * np.sin(phh)).ravel(), ct.ravel()))
    wts = np.repeat(w, nphi) * 2 * np.pi / nphi
    dirs = np.vstack((dirs, D_INC[None, :]))
    F = (matrix(rfar * dirs, y, k) @ c).reshape(-1, 3) * rfar * np.exp(-1j * k * rfar)
    ext = 4 * np.pi / k * float(np.imag(np.conj(P_INC) @ F[-1]))
    sca = float(np.sum(wts * np.sum(np.abs(F[:-1]) ** 2, axis=1)))
    return ext, sca


def col_hmax(k):
    """Largest collocation spacing: 0.012 m, or a fifth of the wavelength."""
    return min(0.012, 2 * np.pi / k / 5)


def solve(k, y, s_col=0.35, cond=1e-13, eta_max=None, ratio=2.5):
    eta_max = eta_max if eta_max is not None else ETA_F
    x, w, nrm = collocation(eta_max, s_col, col_hmax(k))
    if 2 * len(x) < ratio * len(COMPS) * len(y):      # keep rows / unknowns >= ratio
        s_col *= np.sqrt(2 * len(x) / (ratio * len(COMPS) * len(y)))
        x, w, nrm = collocation(eta_max, s_col, col_hmax(k))
    T = tangent_basis(nrm)
    sw = np.sqrt(w / w.sum())
    t0 = time.time()
    A = matrix(x, y, k, T) * np.repeat(sw, 2)[:, None]
    b = (-np.einsum("mab,mb->ma", T, incident(x, k)) * sw[:, None]).reshape(-1)
    scale = np.linalg.norm(A, axis=0)
    A /= scale
    t1 = time.time()
    Q, R = qr(A, mode="economic", overwrite_a=True)       # plain QR; rank-revealing only if needed
    dR = np.abs(np.diag(R))
    if dR.min() > cond * dR.max():
        c, rank = solve_triangular(R, Q.conj().T @ b), len(dR)
    else:
        c, _, rank, _ = lstsq(Q @ R, b, lapack_driver="gelsy", cond=cond)
    return c / scale, {"rows": len(b), "unknowns": len(scale), "rank": int(rank), "N": len(y), "s_col": float(s_col),
                       "t_assemble": t1 - t0, "t_solve": time.time() - t1, "coef_norm": float(np.linalg.norm(c))}


def boundary_residual(k, y, c, eta_max, s=0.27):
    x, w, nrm = collocation(eta_max, s, col_hmax(k))
    T = tangent_basis(nrm)
    rhs = -np.einsum("mab,mb->ma", T, incident(x, k))
    res = (matrix(x, y, k, T) @ c).reshape(-1, 2) - rhs
    return float(np.sqrt(np.sum(w * np.sum(np.abs(res) ** 2, 1)) / np.sum(w * np.sum(np.abs(rhs) ** 2, 1))))


def bem_reference(k, pts, sizes=((40, 36), (60, 54))):
    import bempp_cl.api as bempp
    import bem_reference as br
    from pec_almond import mesh
    br.D, br.P = D_INC, P_INC
    out = []
    for n_t, n_s in sizes:
        v, e = mesh(n_t, n_s)
        t0 = time.time()
        E, ndof = br.scattered_field(bempp.Grid(v.T, e.T.astype("uint32")), k, pts)
        out.append((ndof, E, time.time() - t0))
    return out


if __name__ == "__main__":
    freq = float(sys.argv[1]) if len(sys.argv) > 1 else 1.19
    k = 2 * np.pi * freq * 1e9 / 299792458.0
    pts = 0.3 * fibonacci_sphere(400)
    print(f"f = {freq} GHz, k = {k:.3f}, almond length = {L * k / (2 * np.pi):.2f} wavelengths", flush=True)
    configs = [dict(n_mid=4, J=1, q=0.5), dict(n_mid=6, J=3, q=0.5), dict(n_mid=8, J=5, q=0.5),
               dict(n_mid=10, J=7, q=0.5)]
    for a in sys.argv[2:]:
        if a.startswith("--configs="):          # e.g. --configs=6:3,8:5  (n_mid:J)
            configs = [dict(n_mid=int(p.split(":")[0]), J=int(p.split(":")[1]), q=0.5) for p in a[10:].split(",")]
    results, fields = [], []
    for cfg in configs:
        y = source_rows(**cfg, hmax=2 * np.pi / k / 4)
        c, info = solve(k, y)
        info.update(cfg, bc=boundary_residual(k, y, c, ETA_F, 0.77 * info["s_col"]))
        info["sigma_ext"], info["sigma_sca"] = cross_sections(k, y, c)
        info["optical"] = abs(info["sigma_ext"] - info["sigma_sca"]) / info["sigma_sca"]
        E = (matrix(pts, y, k) @ c).reshape(-1, 3)
        if fields:
            info["change"] = float(np.linalg.norm(E - fields[-1]) / np.linalg.norm(E))
        fields.append(E)
        results.append(info)
        print(f"{cfg}: N={info['N']} unk={info['unknowns']} rows={info['rows']} rank={info['rank']} "
              f"bc={info['bc']:.2e} change={info.get('change', float('nan')):.2e} |c|={info['coef_norm']:.1e} "
              f"asm={info['t_assemble']:.0f}s solve={info['t_solve']:.0f}s\n      sigma_ext={info['sigma_ext']:.6e} "
              f"sigma_sca={info['sigma_sca']:.6e} optical={info['optical']:.1e}", flush=True)
        np.save(f"almond_skeleton_field_{freq}GHz.npy", np.array(fields))
        json.dump(results, open(f"results_pec_almond_skeleton_{freq}GHz.json", "w"), indent=1)
    bem_sizes = [a[6:] for a in sys.argv if a.startswith("--bem=")]      # e.g. --bem=40x36,60x54,70x63
    if "--bem" in sys.argv or bem_sizes:
        sizes = (tuple(tuple(int(v) for v in p.split("x")) for p in bem_sizes[0].split(","))
                 if bem_sizes else ((40, 36), (60, 54)))
        for ndof, Eb, t in bem_reference(k, pts, sizes):
            line = f"BEM dofs={ndof} ({t:.0f}s):"
            for r, E in zip(results, fields):
                line += f" {np.linalg.norm(E - Eb) / np.linalg.norm(Eb):.2e}"
            print(line + "   (difference from each MFS run)", flush=True)
            np.save(f"almond_bem_field_{freq}GHz_{ndof}.npy", Eb)
