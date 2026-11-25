"""
Scattering by the perfectly conducting regular tetrahedron of unit edge length
(Section 5.6 of the paper), with vertices (a, a, a), (a, -a, -a), (-a, a, -a)
and (-a, -a, a), a = 1 / (2 sqrt(2)), together with zero-symmetry convex tetrahedra.

By Theorem 1 the sources lie on the skeleton, the six triangular bisector walls
joining the edges to the insphere centre (the origin for the regular tetrahedron):
  * edge lines parallel to each edge at the depths d_j = 0.14 q^j down to d_min = 0.14 · 2^{-J},
    meeting the 0.56-scaled inner tetrahedron (edge depth (1 - 0.56)/(2 sqrt(2)) = 0.156),
    each with an electric and a magnetic dipole density (Legendre modes in a stretched variable);
  * vertex points on the four rays from the centre to the vertices, and a uniform grid of
    electric dipoles on the faces of an inner tetrahedron;
  * singular tails: nine sub-lines below d_min weighted by the wedge exponents gamma = pi / Phi
    (Phi = 2 pi - arccos(1/3), irrational), and radial nine-point tails on the central ray and
    three skeleton-wall rays per edge at each vertex, weighted by the spherical corner exponents
    nu = (0.32221, 1.09227) and tau = 0.907 (Section 5.6 and Appendix A).

The three half-turns about the coordinate axes form D_2 (four parity classes, each solved on one
equilateral face), and the one-third turn about (1, 1, 1) leaves at most two matrices to factorise.
"""
import itertools
import os
import time

import numpy as np
import scipy.linalg as sla
from scipy.linalg import blas

from . import cone, lines
from .core import class_part, classes_of, cross_sections, factor, plane_wave, sign, sphere_rule, tangent_basis
from .cube import ROT, Skeleton, graded_uv, orbits

H = 0.5                                                              # half edge length
A_TET = 0.5 / np.sqrt(2.0)                                           # vertex coordinate for unit edge length
CLASSES = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]         # parities under the three half-turns of D_2
REFL = np.array(CLASSES, float)                                      # proper 180-degree rotations (det = +1)
VERTS = A_TET * REFL                                                 # v_0 = (a, a, a) and its three images on F_0
FACE0 = np.array([VERTS[1], VERTS[3], VERTS[2]])                     # fundamental face, outward normal -(1, 1, 1)/sqrt(3)
PHI = 2 * np.pi - np.arccos(1.0 / 3.0)                               # exterior wedge angle (289.47 degrees)
GAMMA = np.pi / PHI                                                  # wedge exponent pi / Phi = 0.621823
TAIL_E = {"e": (GAMMA, GAMMA + 1.0), "m": (2 * GAMMA, 2 * GAMMA + 1.0)}
NU, TAU = (0.32221, 1.09227), (0.907,)                               # corner exponents from spherical MPS (mfs.cone)
RULE = sphere_rule(40, 80)


class Line:
    """Dipole density on the edge bisector through +a e_ax at depth d (or 9-line wedge tail with tail=True)"""
    def __init__(self, ax, d, n, kind, tail=False):
        self.ax, self.d, self.n, self.kind, self.tail = ax, d, n, kind, tail
        t_hat = np.zeros(3)
        t_hat[(ax + 1) % 3] = t_hat[(ax + 2) % 3] = 1.0 / np.sqrt(2.0)
        b_hat = np.zeros(3)
        b_hat[ax] = -1.0
        self.t_hat = t_hat
        self.rot = np.ascontiguousarray(np.vstack((t_hat, b_hat, np.cross(t_hat, b_hat))))
        self.r_act = np.array([r for r in range(4) if REFL[r, (ax + 1) % 3] > 0], np.int64)
        self.lam = np.log((1.0 - 1.5 * d) / (1.5 * d))          # |t| <= 1/2 - 1.5 d stays inside the skeleton wall (1.5 > sqrt(2))
        r = 0.5 ** np.arange(1, 10) if tail else np.ones(1)
        mid0 = np.zeros(3)
        mid0[ax] = A_TET
        self.mids = np.ascontiguousarray(mid0[None, :] + (d * r)[:, None] * b_hat[None, :])
        self.w = lines.geom_weights(TAIL_E[kind] + (2.0,)).T[None, :, :] if tail else None

    def keep(self, cls):
        """Columns kept by parity class cls under the half-turn REFL[1 + self.ax] that maps the edge to itself"""
        chi = sign(cls, REFL[1 + self.ax])
        par = (-1.0) ** np.arange(self.n)[:, None]
        kept = chi * par * np.array([[-1.0, 1.0, -1.0]]) > 0
        if self.tail:
            kept[:, 2] = False
        return kept


def _tri_kites(fverts, su, sv, wuv, nu1s):
    """Three graded kites covering triangle fverts (3, 3), with corner weights r_v^{1 - 2 nu1s[i]}"""
    c0, U, V, W2 = fverts.mean(axis=0), 1.0 - 2.0 * su, 1.0 - 2.0 * sv, 4.0 * wuv
    pts, wts, rfac = [], [], []
    for i in range(3):
        A, M1, M2 = fverts[i], 0.5 * (fverts[i] + fverts[(i + 1) % 3]), 0.5 * (fverts[i] + fverts[(i + 2) % 3])
        pk = ((1 - U) * (1 - V))[:, None] * A + ((1 - U) * V)[:, None] * M1 + (U * (1 - V))[:, None] * M2 + (U * V)[:, None] * c0
        du = -(1 - V)[:, None] * A - V[:, None] * M1 + (1 - V)[:, None] * M2 + V[:, None] * c0
        dv = -(1 - U)[:, None] * A + (1 - U)[:, None] * M1 - U[:, None] * M2 + U[:, None] * c0
        pts.append(pk)
        wts.append(np.linalg.norm(np.cross(du, dv), axis=1) * W2)
        rfac.append(np.linalg.norm(pk - A, axis=1) ** (1.0 - 2.0 * nu1s[i]))
    pts = np.vstack(pts)
    dists = [np.linalg.norm(pts - (fverts[i] + np.clip(((pts - fverts[i]) @ (fverts[(i + 1) % 3] - fverts[i]))
             / float(np.sum((fverts[(i + 1) % 3] - fverts[i]) ** 2)), 0.0, 1.0)[:, None]
             * (fverts[(i + 1) % 3] - fverts[i])), axis=1) for i in range(3)]
    return pts, np.concatenate(wts), np.min(dists, axis=0), np.concatenate(rfac)


class Tetrahedron(Skeleton):
    """
    MFS for the regular tetrahedron of unit edge length at level J (smallest depth 0.14 · 2^{-J}).
      q      : depth ratio, 2^{-1/2} (two per octave) or 1/2
      tails  : with or without the edge and corner singular tails
      J=None : plain inner-tetrahedron layout (n_base barycentric order per face) that stalls
    """
    refl = REFL
    r_far = A_TET * np.sqrt(3.0) + 2.0 * H

    def __init__(self, J, k=2 * np.pi, q=2 ** -0.5, extra_modes=None, tails=True, classes=CLASSES, n_base=12):
        self.k, self.classes, self.tails = k, [tuple(c) for c in classes], bool(tails and J is not None)
        s = max(0.0, k / (2 * np.pi) - 1.0)
        m = max(1, int(round(np.log(0.5) / np.log(q))))
        if J is None:
            self.depths, depths_v, ex = np.zeros(0), np.zeros(0), 0.0
        else:
            self.depths = 0.14 * q ** np.arange(m * J + 1)
            mv = max(2, m)                                       # vertex ray uses at least two dipoles per octave
            depths_v = 0.14 * 2.0 ** (-np.arange(mv * J + 1) / mv)
            n0 = max(4, 2 * J - 4)
            ex = 2 * s * (s + 1) if extra_modes is None else extra_modes

        # Uniform interior barycentric grid on an inner tetrahedron (scale 0.56), plus free vertex points
        bary = np.array([(i + 1, j + 1, n_base - i - j) for i in range(n_base) for j in range(n_base - i)], float) / (n_base + 2)
        base = 0.56 * (bary @ FACE0)
        v_hat = np.ones(3) / np.sqrt(3.0)
        vertex = VERTS[0] - depths_v[:, None] * v_hat[None, :]
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
                    ray_dirs.append(cone.unit(t_hat - (fr / np.sqrt(2.0)) * b_hat))
            r_steps = 0.5 ** np.arange(1, 10)
            self.we = lines.geom_weights(np.array(NU) + 2.0)
            self.wm = lines.geom_weights(np.array(TAU) + 2.0)
            self.ray_pts = np.array([VERTS[0][None, :] - (1.2 * self.depths[-1] * r_steps)[:, None] * rd[None, :] for rd in ray_dirs])
            n_ray_slots = len(ray_dirs) * (len(self.we) + len(self.wm))

        sizes = [len(self.pe) + len(self.pm)] + [ln.n for ln in self.lines] + ([n_ray_slots] if n_ray_slots else [])
        self.slot0 = np.cumsum([0] + sizes)
        self.keep = {c: np.concatenate([np.ones((sizes[0], 3), bool)] + [ln.keep(c) for ln in self.lines] +
                                       ([np.ones((n_ray_slots, 3), bool)] if n_ray_slots else [])).ravel() for c in CLASSES}
        self.unknowns = {c: int(kp.sum()) for c, kp in self.keep.items()}

        # Anisotropic graded kites on F_0; since nu_1 = 0.32221 < 1/2, weighting by
        # dist · r_v^{1 - 2 nu_1} keeps |E|^2 integrable at both edges and vertices
        ng = 2 + max(m, 1 + int(np.ceil(0.5 * s)))
        d_min = self.depths[-1] if len(self.depths) else 0.14
        levels = max(2, int(np.ceil(np.log2(2 * H / d_min))))
        su, sv, wuv = graded_uv(levels, ng, min(H / 4, np.pi / (2 * k)))
        self.x, W, dist, rfac = _tri_kites(FACE0, su, sv, wuv, [NU[0]] * 3)
        self.T = tangent_basis(np.tile(-v_hat, (len(self.x), 1)))
        self.rw = np.sqrt(W * dist * rfac)
        self.rows = 2 * len(self.x)

        self.fac, self.t_assemble, self.t_factor = {}, 0.0, 0.0
        if not self.classes:
            return
        t0 = time.time()
        A = self.matrix(self.x, self.T * self.rw[:, None, None], self.classes)
        self.t_assemble, t0 = time.time() - t0, time.time()
        self.fac = {c: factor(A[i, :, :self.unknowns[c]], overwrite=True) for i, c in enumerate(self.classes)}
        self.t_factor = time.time() - t0

    def _add_vertex_tails(self, x, T, A, cmap, sgn):
        lo = 3 * self.slot0[-2]
        ne_cols = len(self.ray_pts) * len(self.we) * 3
        lines.radial_tails_add(x, T, REFL, self.ray_pts, self.we, self.wm, self.k, A,
                               np.ascontiguousarray(cmap[:, lo]), np.ascontiguousarray(cmap[:, lo + ne_cols]), sgn)

    def _unpack_points(self, cl, c_all):
        ncls, ne, nm = len(cl), len(self.pe), len(self.pm)
        n_ray = self.ray_pts.shape[0] * self.ray_pts.shape[1] if self.tails else 0
        y_pt = np.vstack((self.pe, self.ray_pts.reshape(-1, 3))) if n_ray else self.pe
        alpha = np.zeros((ncls, ne + n_ray, 3), complex)
        beta = np.zeros((ncls, ne + n_ray, 3), complex)
        alpha[:, :ne] = c_all[:, :3 * ne].reshape(ncls, ne, 3)
        if nm:
            beta[:, ne - nm:ne] = c_all[:, 3 * ne:3 * (ne + nm)].reshape(ncls, nm, 3)
        if n_ray:
            R_cnt, ne_exp, nm_exp = self.ray_pts.shape[0], len(self.we), len(self.wm)
            lo = 3 * self.slot0[-2]
            ne_cols, nm_cols = R_cnt * ne_exp * 3, R_cnt * nm_exp * 3
            ce_ray = c_all[:, lo:lo + ne_cols].reshape(ncls, R_cnt, ne_exp, 3)
            cm_ray = c_all[:, lo + ne_cols:lo + ne_cols + nm_cols].reshape(ncls, R_cnt, nm_exp, 3)
            alpha[:, ne:] = np.einsum("lrec,et->lrtc", ce_ray, self.we).reshape(ncls, n_ray, 3)
            beta[:, ne:] = np.einsum("lrmc,mt->lrtc", cm_ray, self.wm).reshape(ncls, n_ray, 3)
        return np.ascontiguousarray(y_pt), np.ascontiguousarray(alpha), np.ascontiguousarray(beta)


def solve(J, d=(0, 0, -1.0), p=(0, 1.0, 0), pts=None, k=2 * np.pi, max_gb=None, **layout):
    """
    One plane wave on the regular tetrahedron at level J, using D_2 and the one-third turn about (1, 1, 1)
      pts: evaluation points, shape (P, 3), or None
      Returns the scattered field at pts, the wall time, and a dictionary of checks
    """
    if max_gb is None:
        pages = os.sysconf("SC_PHYS_PAGES") if "SC_PHYS_PAGES" in getattr(os, "sysconf_names", ()) else -1
        ram = pages * os.sysconf("SC_PAGE_SIZE") / 2 ** 30 if pages > 0 else 12.0
        try:
            cg = open("/sys/fs/cgroup/memory.max").read().strip()
            ram = min(ram, int(cg) / 2 ** 30) if cg.isdigit() else ram
        except OSError:
            pass
        max_gb = min(6.0, 0.60 * ram)
    d, p = np.asarray(d, float), np.asarray(p, float)
    pts = None if pts is None else np.asarray(pts, float)
    t0 = time.time()
    orbit = orbits(classes_of(k, d, p, REFL))
    todo, E, sols = list(orbit), (None if pts is None else np.zeros((len(pts), 3), complex)), []
    r2, b2, info = 0.0, 0.0, dict(unknowns=[], t_assemble=0.0, t_factor=0.0)
    probe = Tetrahedron(J, k, classes=[], **layout)
    nmax = max(probe.unknowns.values())
    rw2 = np.repeat(probe.rw, 2)
    Tw = probe.T * probe.rw[:, None, None]
    g = max(1, int((max_gb * 2 ** 30 - 16.0 * nmax ** 2) / (16.0 * probe.rows * nmax)))
    for i in range(0, len(todo), g):
        batch_cls = todo[i:i + g]
        ta0 = time.time()
        A = probe.matrix(probe.x, Tw, batch_cls)
        info["t_assemble"] += time.time() - ta0
        by_rot = {j: {} for j in range(3)}
        tf0 = time.time()
        for ic, c in enumerate(batch_cls):
            fac_c = factor(A[ic, :, :probe.unknowns[c]], overwrite=True)
            for j in orbit[c]:
                Q = np.linalg.matrix_power(ROT, j)
                b = -np.einsum("mab,mb->ma", probe.T, class_part(probe.x, k, Q.T @ d, Q.T @ p, REFL, c)).reshape(-1) * rw2
                by_rot[j][c], res = fac_c.solve(b, residual=True)
                r2, b2 = r2 + float(res) ** 2, b2 + float(np.linalg.norm(b) ** 2)
            del fac_c
        info["t_factor"] += time.time() - tf0
        del A
        bsols = [(np.linalg.matrix_power(ROT, j), coef) for j, coef in by_rot.items() if coef]
        if pts is not None and bsols:
            for Q, coef in bsols:
                E += probe.field(coef, pts @ Q) @ Q.T
        sols.extend(bsols)
        info["unknowns"] += [probe.unknowns[c] for c in batch_cls]
    t = time.time() - t0

    far = lambda dirs: sum(probe.far(coef, dirs @ Q) @ Q.T for Q, coef in sols)
    info["ext"], info["sca"] = cross_sections(far, k, d, p, RULE)
    info["residual"] = float(np.sqrt(r2 / b2))
    info["optical"] = abs(info["ext"] - info["sca"]) / info["sca"]
    info["rows"] = probe.rows
    return E, t, info


# Zero-symmetry convex tetrahedron (six distinct dihedral angles, four distinct corner cones)

SCALENE_VERTS = np.array([
    [0.42, 0.35, 0.30], [0.36, -0.40, -0.31], [-0.44, 0.31, -0.36], [-0.34, -0.26, 0.37],
])
SCALENE_EXP = [                                                      # 5-digit spherical MPS exponents for SCALENE_VERTS
    ((0.31823, 1.08018), (0.92206,)), ((0.31917, 1.08018), (0.92247,)),
    ((0.30948, 1.07668), (0.92297,)), ((0.34176, 1.12037), (0.90798,)),
]
_REFL1, _RACT1, _SGN1, _CONE_CACHE = np.ones((1, 3)), np.zeros(1, np.int64), np.ones((1, 1)), {}


class ScaleneTetrahedron:
    """
    Zero-symmetry MFS on an arbitrary convex tetrahedron verts (4, 3): the six edge-bisector planes
    meet at the insphere centre c_in, forming six flat triangular walls with distinct wedge exponents
    gamma_e = pi / Phi_e and four trihedral corners with distinct spherical MPS exponents (nu_v, tau_v).
    """
    def __init__(self, J, k=2 * np.pi, q=2 ** -0.5, tails=True, verts=SCALENE_VERTS, n_base=12):
        self.k, self.tails, self.verts = float(k), bool(tails and J is not None), np.asarray(verts, float)
        V = self.verts
        fnorms, areas = [], []
        for f in range(4):
            oth = [V[j] for j in range(4) if j != f]
            cr = np.cross(oth[1] - oth[0], oth[2] - oth[0])
            nf = cr / np.linalg.norm(cr)
            fnorms.append(-nf if nf @ (V[f] - oth[0]) > 0 else nf)
            areas.append(0.5 * np.linalg.norm(cr))
        self.fnorms, areas = np.array(fnorms), np.array(areas)
        self.c_in = (areas[:, None] * V).sum(axis=0) / areas.sum()
        self.r_in = float((V[1] - self.c_in) @ self.fnorms[0])
        d0 = 0.68 * self.r_in

        edges_info = []
        for a, b in itertools.combinations(range(4), 2):
            fc, fd = [f for f in range(4) if f not in (a, b)]
            L = float(np.linalg.norm(V[b] - V[a]))
            t_hat = (V[b] - V[a]) / L
            b_hat = cone.unit(-(self.fnorms[fc] + self.fnorms[fd]))
            gamma = np.pi / (2 * np.pi - float(np.arccos(np.clip(-self.fnorms[fc] @ self.fnorms[fd], -1.0, 1.0))))
            D_e = float((self.c_in - V[a]) @ b_hat)
            tan_a, tan_b = D_e / float((self.c_in - V[a]) @ t_hat), D_e / float((V[b] - self.c_in) @ t_hat)
            F = np.ascontiguousarray(np.vstack((t_hat, b_hat, np.cross(t_hat, b_hat))))
            edges_info.append((a, b, L, 0.5 * (V[a] + V[b]), t_hat, b_hat, F, gamma, D_e, tan_a, tan_b))

        L_mean = float(np.mean([e[2] for e in edges_info]))
        s = max(0.0, k * L_mean / (2 * np.pi) - 1.0)
        m = max(1, int(round(np.log(0.5) / np.log(q))))
        if J is None:
            self.depths, depths_v, ex = np.zeros(0), np.zeros(0), 0.0
        else:
            self.depths = d0 * q ** np.arange(m * J + 1)
            mv = max(2, m)
            depths_v = d0 * 2.0 ** (-np.arange(mv * J + 1) / mv)
            n0, ex = max(4, 2 * J - 4), 2 * s * (s + 1)

        bary = np.array([(i + 1, j + 1, n_base - i - j) for i in range(n_base) for j in range(n_base - i)], float) / (n_base + 2)
        base = np.vstack([self.c_in + 0.56 * (bary @ V[[j for j in range(4) if j != f]] - self.c_in) for f in range(4)])
        u_verts = [cone.unit(V[v] - self.c_in) for v in range(4)]
        v_pts = np.vstack([V[v][None, :] - depths_v[:, None] * u_verts[v][None, :] for v in range(4)]) if len(depths_v) else np.zeros((0, 3))
        self.pe, self.pm = np.vstack((base, v_pts)), v_pts

        self.lines = []
        for d in self.depths:
            for a, b, L, mid0, t_hat, b_hat, F, gamma, D_e, tan_a, tan_b in edges_info:
                al = 0.53 * d / D_e
                lam = np.log((1.0 - al) / al)
                n = n0 + int(round(ex * L * (1.0 - 2.0 * al))) + int(np.ceil(2 * np.log(L / d)))
                mids = np.ascontiguousarray((mid0 + (d / D_e) * (self.c_in - mid0))[None, :])
                self.lines.append((0.5 * L, t_hat, mids, lam, d, n, F, np.ones((2, 1, 1, 3)), False))
        if self.tails:
            d_min, r_st = self.depths[-1], 0.5 ** np.arange(1, 10)
            for (a, b, L, mid0, t_hat, b_hat, F, gamma, D_e, tan_a, tan_b), ln_last in zip(edges_info, self.lines[-6:]):
                n, lam = ln_last[5], ln_last[3]
                q_e = ln_last[2][0] - d_min * b_hat                  # shift purely along b_hat so Meixner modes stay unsheared
                mids = np.ascontiguousarray(q_e[None, :] + (d_min * r_st)[:, None] * b_hat[None, :])
                w_pair = np.stack([lines.geom_weights((e0, e1, 2.0)).T[None, :, :]
                                   for e0, e1 in ((gamma, gamma + 1.0), (2 * gamma, 2 * gamma + 1.0))])
                self.lines.append((0.5 * L, t_hat, mids, lam, d_min, n, F, np.ascontiguousarray(w_pair), True))

        if np.allclose(V, SCALENE_VERTS, atol=1e-12):
            self.v_exp = SCALENE_EXP
        else:
            self.v_exp = []
            for v in range(4):
                Vv = cone.unit([V[j] - V[v] for j in range(4) if j != v])
                key = tuple(np.round(Vv.ravel(), 8))
                if key not in _CONE_CACHE:
                    cv = cone.CubeVertex(n=12, ni=150, V=Vv)
                    _CONE_CACHE[key] = (tuple(cv.exponents(False, hi=1.35)[:2]), (cv.exponents(True, hi=1.1)[0],))
                self.v_exp.append(_CONE_CACHE[key])

        self.v_tails, n_vt_cols = [], 0
        if self.tails:
            r_st = 0.5 ** np.arange(1, 10)
            for v in range(4):
                nu_v, tau_v = self.v_exp[v]
                ray_dirs = [u_verts[v]]
                for a, b, L, mid0, t_hat, b_hat, F, gamma, D_e, tan_a, tan_b in edges_info:
                    if v in (a, b):
                        t_out, tan_ev = (-t_hat, tan_a) if v == a else (t_hat, tan_b)
                        for fr in (0.8, 0.4, 0.2):
                            ray_dirs.append(cone.unit(t_out - (fr * tan_ev) * b_hat))
                we = lines.geom_weights(np.array(nu_v) + 2.0)
                wm = lines.geom_weights(np.array(tau_v) + 2.0)
                rpts = np.ascontiguousarray([V[v][None, :] - (1.2 * self.depths[-1] * r_st)[:, None] * rd[None, :] for rd in ray_dirs])
                nc_e, nc_m = len(ray_dirs) * len(we) * 3, len(ray_dirs) * len(wm) * 3
                self.v_tails.append((rpts, np.ascontiguousarray(we), np.ascontiguousarray(wm), nc_e, nc_m))
                n_vt_cols += nc_e + nc_m

        col = 3 * (len(self.pe) + len(self.pm))
        self.line_cmaps = []
        for *_, n, _, _, is_tail in self.lines:
            cm = np.full((1, 6 * n), -1, np.int64)
            for kind in range(2):
                for kk in range(n):
                    for cc in range(2 if is_tail else 3):
                        cm[0, (kind * n + kk) * 3 + cc] = col
                        col += 1
            self.line_cmaps.append(cm)
        self.vt_col0 = col
        self.n_unknowns = col + n_vt_cols

        ng = 3 + int(np.ceil(0.5 * s))
        d_min = self.depths[-1] if len(self.depths) else d0
        levels = max(2, int(np.ceil(np.log2(L_mean / d_min))))
        su, sv, wuv = graded_uv(levels, ng, min(0.125, np.pi / (2 * k * L_mean)))
        xs, Ts, rws = [], [], []
        for f in range(4):
            idx_v = [j for j in range(4) if j != f]
            xf, wf, df, rf = _tri_kites(V[idx_v], su, sv, wuv, [self.v_exp[j][0][0] for j in idx_v])
            xs.append(xf)
            Ts.append(tangent_basis(np.tile(self.fnorms[f], (len(xf), 1))))
            rws.append(np.sqrt(wf * df * rf))
        self.x, self.T, self.rw = np.ascontiguousarray(np.vstack(xs)), np.ascontiguousarray(np.vstack(Ts)), np.concatenate(rws)
        self.rows = 2 * len(self.x)

    def _assemble_into(self, x_sub, Tw_sub, A_buf):
        A_buf[:] = 0.0
        lines.points_add(x_sub, Tw_sub, _REFL1, self.pe, self.pm, self.k, A_buf, np.zeros(1, np.int64), _SGN1)
        for (H_e, t_hat, mids, lam, dp, n, F, w_pair, is_tail), cm in zip(self.lines, self.line_cmaps):
            lines.line_add(x_sub, Tw_sub, _REFL1, _RACT1, t_hat, mids, lam, dp, n, self.k, H_e, A_buf, cm, _SGN1,
                           rot=F, w=w_pair if is_tail else None)
        c_vt = self.vt_col0
        for rpts, we, wm, nc_e, nc_m in self.v_tails:
            lines.radial_tails_add(x_sub, Tw_sub, _REFL1, rpts, we, wm, self.k, A_buf,
                                   np.array([c_vt], np.int64), np.array([c_vt + nc_e], np.int64), _SGN1)
            c_vt += nc_e + nc_m

    def _unpack(self, coef):
        ne, nm = len(self.pe), len(self.pm)
        y_list, a_list = [self.pe], [coef[:3 * ne].reshape(1, ne, 3)]
        b_pt = np.zeros((1, ne, 3), complex)
        if nm:
            b_pt[0, ne - nm:] = coef[3 * ne:3 * (ne + nm)].reshape(nm, 3)
        b_list = [b_pt]
        c_vt = self.vt_col0
        for rpts, we, wm, nc_e, nc_m in self.v_tails:
            R_cnt, nstep = rpts.shape[0], rpts.shape[1]
            y_list.append(rpts.reshape(-1, 3))
            ce = coef[c_vt:c_vt + nc_e].reshape(1, R_cnt, len(we), 3)
            cm = coef[c_vt + nc_e:c_vt + nc_e + nc_m].reshape(1, R_cnt, len(wm), 3)
            a_list.append(np.einsum("lrec,et->lrtc", ce, we).reshape(1, R_cnt * nstep, 3))
            b_list.append(np.einsum("lrmc,mt->lrtc", cm, wm).reshape(1, R_cnt * nstep, 3))
            c_vt += nc_e + nc_m
        return np.ascontiguousarray(np.vstack(y_list)), np.ascontiguousarray(np.concatenate(a_list, axis=1)), np.ascontiguousarray(np.concatenate(b_list, axis=1))

    def field(self, coef, pts):
        """Scattered field at pts (P, 3), evaluated directly from mode coefficients without forming A"""
        pts = np.ascontiguousarray(pts, dtype=float)
        if not len(pts):
            return np.zeros((0, 3), complex)
        y_pt, al, bt = self._unpack(coef)
        out = np.zeros((1, len(pts), 3), complex)
        lines.eval_points_field(pts, _REFL1, _SGN1, y_pt, al, bt, self.k, out)
        for (H_e, t_hat, mids, lam, dp, n, F, w_pair, is_tail), cm in zip(self.lines, self.line_cmaps):
            c_pair = np.zeros((1, 6 * n), complex)
            act = cm[0] >= 0
            c_pair[0, act] = coef[cm[0, act]]
            lines.eval_line_field(pts, _REFL1, _RACT1, _SGN1, t_hat, mids, lam, dp, n, self.k, H_e,
                                  c_pair, out, rot=F, w=w_pair if is_tail else None)
        return out[0]

    def far(self, coef, dirs):
        """Exact r -> infinity far-field pattern F, with E ~ F e^{ikr}/r, in directions dirs (D, 3)"""
        dirs = np.ascontiguousarray(dirs, dtype=float)
        if not len(dirs):
            return np.zeros((0, 3), complex)
        y_pt, al, bt = self._unpack(coef)
        out = np.zeros((1, len(dirs), 3), complex)
        lines.eval_far_dipoles(dirs, _REFL1, _RACT1, _SGN1, y_pt, al, bt, self.k, out)
        ys_l, al_l, bt_l = [], [], []
        for (H_e, t_hat, mids, lam, dp, n, F, w_pair, _), cm in zip(self.lines, self.line_cmaps):
            c_pair = np.zeros((1, 6 * n), complex)
            act = cm[0] >= 0
            c_pair[0, act] = coef[cm[0, act]]
            tq, W, off, *_ = lines.build_tree(lam, dp, n, H_e, self.k)
            pm_n = int(off[1])
            c_cart = 2.0 * np.einsum("kic,lkmc,cq->lkimq", w_pair.reshape(2, len(mids), 3), c_pair.reshape(1, 2, n, 3), F)
            M0 = np.einsum("pm,lkimq->lkpiq", W[:pm_n], c_cart).reshape(1, 2, -1, 3)
            ys_l.append((mids[None, :, :] + tq[:pm_n, None, None] * t_hat[None, None, :]).reshape(-1, 3))
            al_l.append(M0[:, 0])
            bt_l.append(M0[:, 1])
        if ys_l:
            lines.eval_far_dipoles(dirs, _REFL1, _RACT1, _SGN1, np.ascontiguousarray(np.vstack(ys_l)),
                                   np.ascontiguousarray(np.concatenate(al_l, axis=1)),
                                   np.ascontiguousarray(np.concatenate(bt_l, axis=1)), self.k, out)
        return out[0]


def solve_scalene(J, d=(0, 0, -1.0), p=(0, 1.0, 0), pts=None, k=2 * np.pi, verts=SCALENE_VERTS, max_gb=1.2, **layout):
    """
    One plane wave on a zero-symmetry convex tetrahedron (chunking A^H A face-by-face above max_gb)
      pts  : evaluation points, shape (P, 3), or None
      verts: vertex coordinates, shape (4, 3)
      Returns the scattered field at pts, the wall time, and a dictionary of checks
    """
    d, p = np.asarray(d, float), np.asarray(p, float)
    pts = None if pts is None else np.asarray(pts, float)
    t0 = time.time()
    probe = ScaleneTetrahedron(J, k, verts=verts, **layout)
    M_tot, N = len(probe.x), probe.n_unknowns
    Tw = np.ascontiguousarray(probe.T * probe.rw[:, None, None])
    b = -np.einsum("mab,mb->ma", probe.T, plane_wave(probe.x, k, d, p)).reshape(-1) * np.repeat(probe.rw, 2)
    b2 = float(np.linalg.norm(b) ** 2)

    if 16.0 * probe.rows * N <= max_gb * 2 ** 30:
        ta0 = time.time()
        A = np.moveaxis(np.zeros((probe.rows, N, 1), complex, order="F"), 2, 0)
        probe._assemble_into(probe.x, Tw, A)
        t_asm, tf0 = time.time() - ta0, time.time()
        fac = factor(A[0], overwrite=True)
        coef, res = fac.solve(b, residual=True)
        t_fac = time.time() - tf0
        r_norm = float(res)
        del A, fac
    else:
        G = np.zeros((N, N), np.complex128, order="F")
        rhs = np.zeros(N, np.complex128)
        splits = np.linspace(0, M_tot, 5, dtype=int)
        A_buf = np.moveaxis(np.zeros((2 * int(np.max(np.diff(splits))), N, 1), np.complex128, order="F"), 2, 0)
        ta0, tf_acc = time.time(), 0.0
        for ic in range(4):
            m0, m1 = splits[ic], splits[ic + 1]
            sr = 2 * (m1 - m0)
            A_sub = A_buf[:, :sr, :] if sr == A_buf.shape[1] else np.moveaxis(np.zeros((sr, N, 1), np.complex128, order="F"), 2, 0)
            probe._assemble_into(np.ascontiguousarray(probe.x[m0:m1]), np.ascontiguousarray(Tw[m0:m1]), A_sub)
            tf1 = time.time()
            blas.zherk(1.0, A_sub[0], beta=1.0, c=G, trans=2, lower=0, overwrite_c=1)
            rhs += blas.zgemv(1.0, A_sub[0], b[2 * m0:2 * m1], trans=2)
            tf_acc += time.time() - tf1
        del A_buf, A_sub
        t_asm, tf2 = (time.time() - ta0) - tf_acc, time.time()
        scale = np.sqrt(np.real(np.diag(G)))
        G /= scale[:, None]
        G /= scale[None, :]
        R = sla.cho_factor(G, overwrite_a=True, check_finite=False)
        coef = sla.cho_solve(R, rhs / scale, check_finite=False) / scale
        r_norm = float(np.sqrt(max(0.0, b2 - float(np.real(np.vdot(coef, rhs))))))
        del R, G
        t_fac = tf_acc + (time.time() - tf2)

    E = None if pts is None else probe.field(coef, pts)
    t = time.time() - t0
    ext, sca = cross_sections(lambda dirs: probe.far(coef, dirs), k, d, p, RULE)
    info = dict(unknowns=[N], rows=probe.rows, t_assemble=t_asm, t_factor=t_fac,
                residual=r_norm / np.sqrt(b2), ext=ext, sca=sca, optical=abs(ext - sca) / sca)
    return E, t, info
