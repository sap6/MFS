"""MFS for a perfectly conducting NASA almond with "skeleton" sources.

Every cross-section of the almond is a 3:1 ellipse, so the scattered field
continues analytically inside the body up to (roughly) the focal curve in the
mid-plane z = 0, which lies within 6% of the rim. Sources are electric and
magnetic point dipoles on the mid-plane, in rows that follow the rim and
cluster towards the focal curve. Their spacing along a row is proportional to
the local distance to the surface, as is the spacing of the collocation points.

Incidence is nose-on (travelling in -x, towards the tip), polarised along y,
so the field is odd under y -> -y and even under z -> -z; the unknowns and
collocation points are reduced to one quarter accordingly.

Usage: python pec_almond_skeleton.py [freq_GHz]
"""
import json
import sys
import time

import numpy as np
from numba import njit, prange
from scipy.linalg import lstsq

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


def matrix(x, y, k, T=None):
    nc = 2 if T is not None else 3
    out = np.empty((len(x), nc, len(y), len(COMPS)), complex)
    sym_block(np.ascontiguousarray(x), np.ascontiguousarray(T) if T is not None else np.zeros((1, 2, 3)),
              T is not None, np.ascontiguousarray(y), COMPS, REFL, SIG, k, out)
    return out.reshape(len(x) * nc, len(y) * len(COMPS))


# ---------- sources on the mid-plane -------------------------------------------

def surf_dist(x, y):
    """Approximate distance from the mid-plane point (x, y, 0) to the surface."""
    a, b = ab(x)
    eta = np.minimum(np.abs(y) / np.maximum(a, 1e-300), 1.0)
    return np.minimum(a - np.abs(y), b * np.sqrt(1 - eta ** 2))


def source_rows(n_mid, J, q, beyond=(), s=0.7, tip_cut=0.004):
    etas = list((ETA_F - 3 * GAP) * np.cos((np.arange(n_mid) + 0.5) * np.pi / (2 * n_mid)))
    etas += [ETA_F - 2 * GAP * q ** j for j in range(J + 1)]
    etas += [ETA_F + f * GAP for f in beyond]
    pts = []
    x_end = L * (T1 - tip_cut)
    for eta in etas:
        # front part: y = eta * a(x)
        x = 0.0
        step = s * surf_dist(0.0, eta * ab(0.0)[0])
        x += 0.5 * step
        while x < x_end:
            yv = eta * ab(x)[0]
            pts.append((x, yv))
            x += s * surf_dist(x, yv)
        # rear part: confocal ellipse in the mid-plane
        Amu = np.sqrt(A_ ** 2 - B_ ** 2 * (1 - eta ** 2))
        psi = 0.0
        first = True
        while psi < np.pi / 2:
            xv, yv = -Amu * np.sin(psi), eta * B_ * np.cos(psi)
            dist = surf_dist(xv, yv)
            ds = s * dist
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
        dphi = min(s * rel_dist(phi + 1e-9) / np.sqrt(np.sin(phi) ** 2 + np.cos(phi) ** 2 / 9), 0.25)
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


def incident(x, k):
    return P_INC[None, :] * np.exp(1j * k * (x @ D_INC))[:, None]


def solve(k, y, s_col=0.35, cond=1e-13, eta_max=None):
    eta_max = eta_max if eta_max is not None else ETA_F
    x, w, nrm = collocation(eta_max, s_col)
    T = tangent_basis(nrm)
    sw = np.sqrt(w / w.sum())
    t0 = time.time()
    A = matrix(x, y, k, T) * np.repeat(sw, 2)[:, None]
    b = (-np.einsum("mab,mb->ma", T, incident(x, k)) * sw[:, None]).reshape(-1)
    scale = np.linalg.norm(A, axis=0)
    A /= scale
    t1 = time.time()
    c, _, rank, _ = lstsq(A, b, lapack_driver="gelsy", cond=cond, overwrite_a=True, overwrite_b=True)
    return c / scale, {"rows": len(b), "unknowns": len(scale), "rank": int(rank), "N": len(y),
                       "t_assemble": t1 - t0, "t_solve": time.time() - t1, "coef_norm": float(np.linalg.norm(c))}


def boundary_residual(k, y, c, eta_max, s=0.27):
    x, w, nrm = collocation(eta_max, s)
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
    results, fields = [], []
    for cfg in configs:
        y = source_rows(**cfg)
        c, info = solve(k, y)
        info.update(cfg, bc=boundary_residual(k, y, c, ETA_F))
        E = (matrix(pts, y, k) @ c).reshape(-1, 3)
        if fields:
            info["change"] = float(np.linalg.norm(E - fields[-1]) / np.linalg.norm(E))
        fields.append(E)
        results.append(info)
        print(f"{cfg}: N={info['N']} unk={info['unknowns']} rows={info['rows']} rank={info['rank']} "
              f"bc={info['bc']:.2e} change={info.get('change', float('nan')):.2e} |c|={info['coef_norm']:.1e} "
              f"asm={info['t_assemble']:.0f}s solve={info['t_solve']:.0f}s", flush=True)
        np.save(f"almond_skeleton_field_{freq}GHz.npy", np.array(fields))
        json.dump(results, open(f"results_pec_almond_skeleton_{freq}GHz.json", "w"), indent=1)
    if "--bem" in sys.argv:
        for ndof, Eb, t in bem_reference(k, pts):
            line = f"BEM dofs={ndof} ({t:.0f}s):"
            for r, E in zip(results, fields):
                line += f" {np.linalg.norm(E - Eb) / np.linalg.norm(Eb):.2e}"
            print(line + "   (difference from each MFS run)", flush=True)
            np.save(f"almond_bem_field_{freq}GHz_{ndof}.npy", Eb)
