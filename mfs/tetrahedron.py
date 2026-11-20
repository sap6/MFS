"""
Scattering by the perfectly conducting regular tetrahedron of unit edge length
(Section 5.6 of the paper), with vertices (a, a, a), (a, -a, -a), (-a, a, -a)
and (-a, -a, a), a = 1 / (2 sqrt(2)).

By Theorem 1 the sources lie on the skeleton, the six triangular bisector walls
joining the edges to the origin:
  * edge lines parallel to each edge at the depths d_j = 0.1 q^j down to d_min = 0.1 · 2^{-J},
    each with an electric and a magnetic dipole density (Legendre modes in a stretched variable);
  * vertex points on the four rays from the origin to the vertices, and a uniform grid of
    electric dipoles on the faces of an inner tetrahedron;
  * singular tails: nine sub-lines below d_min weighted by the wedge exponents gamma = pi / Phi
    (Phi = 2 pi - arccos(1/3), irrational), and radial nine-point tails on the central ray and
    three skeleton-wall rays per edge at each vertex, weighted by the spherical corner exponents
    nu = (0.32221, 1.09227) and tau = 0.907 (Section 6.2 and Appendix A).

The three half-turns about the coordinate axes form D_2 (four parity classes, each solved on one
equilateral face), and the one-third turn about (1, 1, 1) leaves at most two matrices to factorise.
"""
import os
import time

import numpy as np
from numba import njit, prange

from . import lines
from .core import class_part, classes_of, cross_sections, factor, sign, sphere_rule, tangent_basis
from .cube import ROT, graded_nodes, orbits

H = 0.5                                                              # half edge length
A_TET = 0.5 / np.sqrt(2.0)                                           # vertex coordinate for unit edge length
CLASSES = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]         # parities under the three half-turns of D_2
REFL = np.array(CLASSES, float)                                      # proper 180-degree rotations (det = +1)
VERTS = A_TET * REFL                                                 # v_0 = (a, a, a) and its three images on F_0
FACE0 = np.array([VERTS[1], VERTS[3], VERTS[2]])                     # fundamental face, outward normal -(1, 1, 1)/sqrt(3)
PHI = 2 * np.pi - np.arccos(1.0 / 3.0)                               # exterior wedge angle (289.47 degrees)
GAMMA = np.pi / PHI                                                  # wedge exponent pi / Phi = 0.621825
TAIL_E = {"e": (GAMMA, GAMMA + 1.0), "m": (2 * GAMMA, 2 * GAMMA + 1.0)}
NU, TAU = (0.32221, 1.09227), (0.907,)                               # corner exponents from spherical MPS (mfs.cone)
RULE = sphere_rule(40, 80)


@njit(parallel=True, fastmath=True, cache=True)
def _line_kernel(x, T, refl, r_act, mid_ax, F, w, tq, W, off, tlo, thi, ltab, g, gw, lam, d, n, k, A, cmap, sgn):
    """Fields of an electric and magnetic skeleton line (or 9-line tail) in local frame F = [t_hat, b_hat, n_hat]"""
    M, npos, ns, nc = x.shape[0], w.shape[2], w.shape[1], T.shape[1]
    tmul = 1 if T.shape[0] > 1 else 0
    a_tet = 0.5 / np.sqrt(2.0)
    r_vec = np.ones(npos)
    if npos > 1:
        for i in range(npos):
            r_vec[i] = 0.5 ** (i + 1)
    nch = min(16, max(1, M))
    for ch in prange(nch):
        accr, acci = np.empty((2, 2, ns, nc, 3, n)), np.empty((2, 2, ns, nc, 3, n))
        S = np.empty((2, ns, nc, 3), np.complex128)
        Ge, Gm = np.empty((nc, 3), np.complex128), np.empty((nc, 3), np.complex128)
        Tm, xm, e, Wq = np.empty((nc, 3)), np.empty(3), np.empty(3), np.empty(n)
        stack_l, stack_i = np.empty(88, np.int64), np.empty(88, np.int64)
        for m in range((M * ch) // nch, (M * (ch + 1)) // nch):
            accr[:] = 0.0
            acci[:] = 0.0
            for ia in range(2):
                r = r_act[ia]
                for j in range(3):
                    xm[j] = x[m, j] * refl[r, j]
                    for a in range(nc):
                        Tm[a, j] = T[m * tmul, a, j] * refl[r, j]
                dx0 = xm[0] - (a_tet if mid_ax == 0 else 0.0)
                dx1 = xm[1] - (a_tet if mid_ax == 1 else 0.0)
                dx2 = xm[2] - (a_tet if mid_ax == 2 else 0.0)
                xt = F[0, 0] * dx0 + F[0, 1] * dx1 + F[0, 2] * dx2
                xb = F[1, 0] * dx0 + F[1, 1] * dx1 + F[1, 2] * dx2
                xn = F[2, 0] * dx0 + F[2, 1] * dx1 + F[2, 2] * dx2
                rho2 = 1e300
                for i in range(npos):
                    di = d * r_vec[i]
                    rho2 = min(rho2, (xb - di) ** 2 + xn ** 2)

                sp, stack_l[0], stack_i[0] = 1, 0, 0
                while sp > 0:
                    sp -= 1
                    lev, idx = stack_l[sp], stack_i[sp]
                    wid = 2.0 / (1 << lev)
                    p = (1 << lev) - 1 + idx if lev <= ltab else -1
                    a0 = tlo[p] if p >= 0 else 0.5 * np.tanh(0.5 * lam * (-1.0 + wid * idx))
                    b0 = thi[p] if p >= 0 else 0.5 * np.tanh(0.5 * lam * (-1.0 + wid * (idx + 1)))
                    gap = max(a0 - xt, 0.0, xt - b0)
                    if b0 - a0 > np.sqrt(rho2 + gap * gap) and lev < 40:
                        stack_l[sp], stack_i[sp], stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx, lev + 1, 2 * idx + 1
                        sp += 2
                        continue
                    for iq in range(off[p + 1] - off[p] if p >= 0 else g.shape[0]):
                        if p >= 0:
                            t = tq[off[p] + iq]
                            for kk in range(n):
                                Wq[kk] = W[off[p] + iq, kk]
                        else:
                            u = -1.0 + wid * (idx + 0.5) + 0.5 * wid * g[iq]
                            ex = np.exp(lam * u)
                            t = 0.5 * (ex - 1.0) / (ex + 1.0)
                            Wq[0] = 0.5 * wid * gw[iq] * lam * ex / ((ex + 1.0) * (ex + 1.0)) / d
                            if n > 1:
                                Wq[1] = Wq[0] * u
                            for kk in range(2, n):
                                Wq[kk] = ((2 * kk - 1) * u * Wq[kk - 1] - (kk - 1) * Wq[kk - 2]) / kk
                        S[:] = 0.0
                        for i in range(npos):
                            di = d * r_vec[i]
                            for q in range(3):
                                e[q] = xm[q] - ((a_tet if q == mid_ax else 0.0) + di * F[1, q] + t * F[0, q])
                            rr = np.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
                            e /= rr
                            lines.dipoles(Tm, e, k, rr, Ge, Gm)
                            for a in range(nc):
                                for cc in range(3):
                                    le = Ge[a, 0] * F[cc, 0] + Ge[a, 1] * F[cc, 1] + Ge[a, 2] * F[cc, 2]
                                    lm = Gm[a, 0] * F[cc, 0] + Gm[a, 1] * F[cc, 1] + Gm[a, 2] * F[cc, 2]
                                    for s in range(ns):
                                        S[0, s, a, cc] += w[0, s, i, cc] * le
                                        S[1, s, a, cc] += w[1, s, i, cc] * lm
                        for kind in range(2):
                            for s in range(ns):
                                for a in range(nc):
                                    for cc in range(3):
                                        sr, si = S[kind, s, a, cc].real, S[kind, s, a, cc].imag
                                        for kk in range(n):
                                            accr[ia, kind, s, a, cc, kk] += sr * Wq[kk]
                                            acci[ia, kind, s, a, cc, kk] += si * Wq[kk]
            for cl in range(sgn.shape[0]):
                for kind in range(2):
                    for s in range(ns):
                        for kk in range(n):
                            for cc in range(3):
                                col = cmap[cl, ((kind * ns + s) * n + kk) * 3 + cc]
                                if col >= 0:
                                    for a in range(nc):
                                        vr, vi = 0.0, 0.0
                                        for ia in range(2):
                                            r = r_act[ia]
                                            vr += sgn[cl, r] * accr[ia, kind, s, a, cc, kk]
                                            vi += sgn[cl, r] * acci[ia, kind, s, a, cc, kk]
                                        A[cl, m * nc + a, col] += 2.0 * (vr + 1j * vi)


@njit(parallel=True, fastmath=True, cache=True)
def _corner_tails_add(x, T, refl, ray_pts, we, wm, k, A, col_e, col_m, sgn):
    """Fields of the radial corner-tail rays (9 tied dipoles per ray, sharing electric/magnetic kernel calls)"""
    M, nc = x.shape[0], T.shape[1]
    ncls, nr = sgn.shape
    R_cnt, nstep = ray_pts.shape[0], ray_pts.shape[1]
    ne, nm = we.shape[0], wm.shape[0]
    tmul = 1 if T.shape[0] > 1 else 0
    for m in prange(M):
        acce = np.empty((ncls, ne, nc, 3), np.complex128)
        accm = np.empty((ncls, nm, nc, 3), np.complex128)
        Ge, Gm = np.empty((nc, 3), np.complex128), np.empty((nc, 3), np.complex128)
        Tm, e = np.empty((nc, 3)), np.empty(3)
        for ir in range(R_cnt):
            acce[:] = 0.0
            accm[:] = 0.0
            for r in range(nr):
                for a in range(nc):
                    for q in range(3):
                        Tm[a, q] = T[m * tmul, a, q] * refl[r, q]
                for t in range(nstep):
                    for q in range(3):
                        e[q] = x[m, q] * refl[r, q] - ray_pts[ir, t, q]
                    rr = np.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
                    e /= rr
                    lines.dipoles(Tm, e, k, rr, Ge, Gm)
                    for s in range(ncls):
                        sg = sgn[s, r]
                        for ie in range(ne):
                            wv = sg * we[ie, t]
                            for a in range(nc):
                                for cc in range(3):
                                    acce[s, ie, a, cc] += wv * Ge[a, cc]
                        for im in range(nm):
                            wv = sg * wm[im, t]
                            for a in range(nc):
                                for cc in range(3):
                                    accm[s, im, a, cc] += wv * Gm[a, cc]
            for s in range(ncls):
                be = col_e[s] + ir * ne * 3
                for ie in range(ne):
                    for a in range(nc):
                        for cc in range(3):
                            A[s, m * nc + a, be + ie * 3 + cc] += acce[s, ie, a, cc]
                bm = col_m[s] + ir * nm * 3
                for im in range(nm):
                    for a in range(nc):
                        for cc in range(3):
                            A[s, m * nc + a, bm + im * 3 + cc] += accm[s, im, a, cc]


class Line:
    """Dipole density on the edge bisector through +a e_ax at depth d (or 9-line wedge tail with tail=True)"""
    def __init__(self, ax, d, n, kind, tail=False):
        self.ax, self.d, self.n, self.kind, self.tail = ax, d, n, kind, tail
        t_hat = np.zeros(3)
        t_hat[(ax + 1) % 3] = t_hat[(ax + 2) % 3] = 1.0 / np.sqrt(2.0)
        b_hat = np.zeros(3)
        b_hat[ax] = -1.0
        self.F = np.ascontiguousarray(np.vstack((t_hat, b_hat, np.cross(t_hat, b_hat))))
        self.r_act = np.array([r for r in range(4) if REFL[r, (ax + 1) % 3] > 0], np.int64)
        self.lam = np.log((1.0 - 1.5 * d) / (1.5 * d))          # |t| <= 1/2 - 1.5 d stays inside the skeleton wall (1.5 > sqrt(2))
        if tail:
            r = 0.5 ** np.arange(1, 10)
            e = np.array(TAIL_E[kind] + (2.0,))
            self.w = r[None, :, None] ** e[None, None, :]
            self.w[:, -1, :] /= 1.0 - 0.5 ** e
        else:
            self.w = np.ones((1, 1, 3))

    def keep(self, cls):
        """Columns kept by parity class cls under the half-turn REFL[1 + ax] that maps the edge to itself"""
        chi = sign(cls, REFL[1 + self.ax])
        par = (-1.0) ** np.arange(self.n)[:, None]
        kept = chi * par * np.array([[-1.0, 1.0, -1.0]]) > 0
        if self.tail:
            kept[:, 2] = False
        return kept


def face_kites(sc, w):
    """Three kites covering F_0, graded towards all three edges and vertices (sc -> 0.5 is xi -> 0)"""
    c0 = FACE0.mean(axis=0)
    xi = 1.0 - 2.0 * sc
    U, V = [g.ravel() for g in np.meshgrid(xi, xi, indexing="ij")]
    W2 = 4.0 * np.outer(w, w).ravel()
    pts, wts = [], []
    for i in range(3):
        A, M1, M2 = FACE0[i], 0.5 * (FACE0[i] + FACE0[(i + 1) % 3]), 0.5 * (FACE0[i] + FACE0[(i + 2) % 3])
        pts.append(((1 - U) * (1 - V))[:, None] * A + ((1 - U) * V)[:, None] * M1 + (U * (1 - V))[:, None] * M2 + (U * V)[:, None] * c0)
        du = -(1 - V)[:, None] * A - V[:, None] * M1 + (1 - V)[:, None] * M2 + V[:, None] * c0
        dv = -(1 - U)[:, None] * A + (1 - U)[:, None] * M1 - U[:, None] * M2 + U[:, None] * c0
        wts.append(np.linalg.norm(np.cross(du, dv), axis=1) * W2)
    pts = np.vstack(pts)
    dists = [np.linalg.norm(pts - (FACE0[i] + np.clip((pts - FACE0[i]) @ (FACE0[(i + 1) % 3] - FACE0[i]), 0.0, 1.0)[:, None]
                                   * (FACE0[(i + 1) % 3] - FACE0[i])), axis=1) for i in range(3)]
    return pts, np.concatenate(wts), np.min(dists, axis=0)


class Tetrahedron:
    """
    MFS for the regular tetrahedron of unit edge length at level J (smallest depth 0.1 · 2^{-J}).
      q      : depth ratio, 2^{-1/2} (two per octave) or 1/2
      tails  : with or without the edge and corner singular tails
      J=None : plain inner-tetrahedron layout (n_base barycentric order per face) that stalls
    """
    def __init__(self, J, k=2 * np.pi, q=2 ** -0.5, extra_modes=None, tails=True, classes=CLASSES, n_base=12):
        self.k, self.classes, self.tails = k, [tuple(c) for c in classes], tails
        s = max(0.0, k / (2 * np.pi) - 1.0)
        m = max(1, int(round(np.log(0.5) / np.log(q))))
        if J is None:
            self.depths, self.tails, ex = np.zeros(0), False, 0.0
        else:
            self.depths = 0.1 * q ** np.arange(m * J + 1)
            n0 = max(4, 2 * J - 4)
            ex = 2 * s * (s + 1) if extra_modes is None else extra_modes

        # Uniform interior barycentric grid on an inner tetrahedron (scale 0.56), plus free vertex points
        bary = np.array([(i + 1, j + 1, n_base - i - j) for i in range(n_base) for j in range(n_base - i)], float) / (n_base + 2)
        base = 0.56 * (bary @ FACE0)
        v_hat = np.ones(3) / np.sqrt(3.0)
        vertex = VERTS[0] - self.depths[:, None] * v_hat[None, :]
        self.pe, self.pm = np.concatenate((base, vertex)), vertex
        self.lines = [Line(ax, d, n0 + int(round(ex * (1.0 - 3.0 * d))) + int(np.ceil(2 * np.log(1.0 / d))), kind)
                      for d in self.depths for ax in range(3) for kind in "em"]
        n_ray_slots = 0
        if self.tails:
            self.lines += [Line(ln.ax, ln.d, ln.n, ln.kind, tail=True) for ln in self.lines[-6:]]
            ray_dirs = [v_hat]
            for ax in range(3):                                 # three radial rays per skeleton wall, clustering towards the edge
                t_hat = np.zeros(3)
                t_hat[(ax + 1) % 3] = t_hat[(ax + 2) % 3] = 1.0 / np.sqrt(2.0)
                b_hat = np.zeros(3)
                b_hat[ax] = -1.0
                for fr in (0.8, 0.4, 0.2):
                    vec = t_hat - (fr / np.sqrt(2.0)) * b_hat
                    ray_dirs.append(vec / np.linalg.norm(vec))
            r_steps = 0.5 ** np.arange(1, 10)
            ee, em = np.array(NU) + 2.0, np.array(TAU) + 2.0
            self.we = r_steps[None, :] ** ee[:, None]
            self.we[:, -1] /= 1.0 - 0.5 ** ee
            self.wm = r_steps[None, :] ** em[:, None]
            self.wm[:, -1] /= 1.0 - 0.5 ** em
            self.ray_pts = np.array([VERTS[0][None, :] - (1.2 * self.depths[-1] * r_steps)[:, None] * rd[None, :] for rd in ray_dirs])
            n_ray_slots = len(ray_dirs) * (len(ee) + len(em))

        sizes = [len(self.pe) + len(self.pm)] + [ln.n for ln in self.lines] + ([n_ray_slots] if n_ray_slots else [])
        self.slot0 = np.cumsum([0] + sizes)
        self.keep = {c: np.concatenate([np.ones((sizes[0], 3), bool)] + [ln.keep(c) for ln in self.lines] +
                                       ([np.ones((n_ray_slots, 3), bool)] if n_ray_slots else [])).ravel() for c in CLASSES}
        self.unknowns = {c: int(kp.sum()) for c, kp in self.keep.items()}

        # Graded boundary points on F_0 (xi = 1 - 2 sc, so H/8 in sc gives four panels per kite side), weighted by sqrt(area · dist)
        ng = 2 + max(m, 1 + int(np.ceil(0.5 * s)))
        d_min = self.depths[-1] if len(self.depths) else 0.1
        levels = max(2, int(np.ceil(np.log2(H / d_min))))
        sc, w = graded_nodes(levels, ng, min(H / 8, np.pi / (4 * k)))
        self.x, W, dist = face_kites(sc, w)
        self.T = tangent_basis(np.tile(-v_hat, (len(self.x), 1)))
        self.rw = np.sqrt(W * dist)
        self.rows = 2 * len(self.x)

        self.fac, self.t_assemble, self.t_factor = {}, 0.0, 0.0
        if not self.classes:
            return
        t0 = time.time()
        A = self.matrix(self.x, self.T * self.rw[:, None, None], self.classes)
        self.t_assemble, t0 = time.time() - t0, time.time()
        self.fac = {c: factor(A[i, :, :self.unknowns[c]], overwrite=True) for i, c in enumerate(self.classes)}
        self.t_factor = time.time() - t0

    def matrix(self, x, T, classes):
        nmax = max(self.unknowns[c] for c in classes)
        cmap = np.full((len(classes), 3 * self.slot0[-1]), -1, np.int64)
        for i, c in enumerate(classes):
            cmap[i, self.keep[c]] = np.arange(self.unknowns[c])
        T = np.eye(3)[None] if T is None else np.ascontiguousarray(T)
        A = np.moveaxis(np.zeros((len(x) * T.shape[1], nmax, len(classes)), complex, order="F"), 2, 0)
        sgn = np.array([[sign(c, R) for R in REFL] for c in classes])
        x = np.ascontiguousarray(x)
        lines.points_add(x, T, REFL, self.pe, self.pm, self.k, A, np.ascontiguousarray(cmap[:, 0]), sgn)
        for i, ln in enumerate(self.lines):
            if ln.kind == "e":
                lo = 3 * self.slot0[1 + i]
                w = np.ascontiguousarray(np.stack((ln.w, self.lines[i + 1].w)))
                tq, W, off, tlo, thi, ltab, g, gw = lines.build_tree(ln.lam, ln.d, ln.n, H, self.k)
                _line_kernel(x, T, REFL, ln.r_act, ln.ax, ln.F, w, tq, W, off, tlo, thi, ltab, g, gw,
                             ln.lam, ln.d, ln.n, self.k, A, np.ascontiguousarray(cmap[:, lo:lo + 6 * ln.n]), sgn)
        if self.tails:
            lo = 3 * self.slot0[-2]
            ne_cols = len(self.ray_pts) * len(self.we) * 3
            _corner_tails_add(x, T, REFL, self.ray_pts, self.we, self.wm, self.k, A,
                              np.ascontiguousarray(cmap[:, lo]), np.ascontiguousarray(cmap[:, lo + ne_cols]), sgn)
        return A

    def solve(self, d, p, only=None):
        coef, r2, b2 = {}, 0.0, 0.0
        for c in classes_of(self.k, d, p, REFL):
            if only is not None and c not in only:
                continue
            b = -np.einsum("mab,mb->ma", self.T, class_part(self.x, self.k, d, p, REFL, c)).reshape(-1) * np.repeat(self.rw, 2)
            coef[c], res = self.fac[c].solve(b, residual=True)
            r2, b2 = r2 + float(res) ** 2, b2 + np.linalg.norm(b) ** 2
        return coef, np.sqrt(r2 / b2), b2

    def field(self, coef, pts, batch=1000):
        pts = np.asarray(pts, float)
        if not len(pts):
            return np.empty((0, 3), complex)
        cl, out = list(coef), []
        for i in range(0, len(pts), batch):
            A = self.matrix(pts[i:i + batch], None, cl)
            out.append(sum(A[j][:, :self.unknowns[c]] @ coef[c] for j, c in enumerate(cl)).reshape(-1, 3))
        return np.concatenate(out)

    def far(self, coef, dirs, rfar=1.0e6):
        return self.field(coef, rfar * dirs) * rfar * np.exp(-1j * self.k * rfar)


def solve(J, d=(0, 0, -1.0), p=(0, 1.0, 0), pts=None, k=2 * np.pi, max_gb=None, **layout):
    """One plane wave at level J using D_2 and the one-third turn about (1, 1, 1) (at most two matrices)"""
    if max_gb is None:
        pages = os.sysconf("SC_PHYS_PAGES") if "SC_PHYS_PAGES" in getattr(os, "sysconf_names", ()) else -1
        ram = pages * os.sysconf("SC_PAGE_SIZE") / 2 ** 30 if pages > 0 else 12.0
        try:
            cg = open("/sys/fs/cgroup/memory.max").read().strip()
            ram = min(ram, int(cg) / 2 ** 30) if cg.isdigit() else ram
        except OSError:
            pass
        max_gb = min(6.0, 0.55 * ram)
    d, p = np.asarray(d, float), np.asarray(p, float)
    pts = None if pts is None else np.asarray(pts, float)
    t0 = time.time()
    orbit = orbits(classes_of(k, d, p, REFL))
    todo, E, sols = list(orbit), (None if pts is None else np.zeros((len(pts), 3), complex)), []
    r2, b2, info = 0.0, 0.0, dict(unknowns=[], t_assemble=0.0, t_factor=0.0)
    probe = Tetrahedron(J, k, classes=[], **layout)
    nmax = max(probe.unknowns.values())
    g = max(1, int(max_gb * 2 ** 30 / (16.0 * (probe.rows * nmax + nmax ** 2))))
    for i in range(0, len(todo), g):
        tet = Tetrahedron(J, k, classes=todo[i:i + g], **layout)
        bsols = []
        for j in range(3):
            only = [c for c in tet.classes if j in orbit[c]]
            if only:
                Q = np.linalg.matrix_power(ROT, j)
                coef, res, nb2 = tet.solve(Q.T @ d, Q.T @ p, only)
                bsols.append((Q, coef))
                r2, b2 = r2 + res ** 2 * nb2, b2 + nb2
        tet.fac.clear()
        if pts is not None and bsols:
            for p0 in range(0, len(pts), 1000):
                pb = pts[p0:p0 + 1000]
                Ap = tet.matrix(np.vstack([pb @ Q for Q, _ in bsols]), None, tet.classes)
                for m, (Q, coef) in enumerate(bsols):
                    sl = slice(3 * m * len(pb), 3 * (m + 1) * len(pb))
                    E[p0:p0 + 1000] += sum(Ap[s, sl, :tet.unknowns[c]] @ coef[c] for s, c in enumerate(tet.classes) if c in coef).reshape(-1, 3) @ Q.T
        sols.extend(bsols)
        info["unknowns"] += [tet.unknowns[c] for c in tet.classes]
        info["t_assemble"] += tet.t_assemble
        info["t_factor"] += tet.t_factor
        del tet
    t = time.time() - t0

    far = lambda dirs: sum(probe.far(coef, dirs @ Q) @ Q.T for Q, coef in sols)
    info["ext"], info["sca"] = cross_sections(far, k, d, p, RULE)
    info["residual"] = float(np.sqrt(r2 / b2))
    info["optical"] = abs(info["ext"] - info["sca"]) / info["sca"]
    info["rows"] = probe.rows
    return E, t, info
