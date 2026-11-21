"""
Scattering by the perfectly conducting regular tetrahedron of unit edge length
(Section 5.6 of the paper), with vertices (a, a, a), (a, -a, -a), (-a, a, -a)
and (-a, -a, a), a = 1 / (2 sqrt(2)).

By Theorem 1 the sources lie on the skeleton, the six triangular bisector walls
joining the edges to the origin:
  * edge lines parallel to each edge at the depths d_j = 0.14 q^j down to d_min = 0.14 · 2^{-J},
    meeting the 0.56-scaled inner tetrahedron (edge depth (1 - 0.56)/(2 sqrt(2)) = 0.156),
    each with an electric and a magnetic dipole density (Legendre modes in a stretched variable);
  * vertex points on the four rays from the origin to the vertices, and a uniform grid of
    electric dipoles on the faces of an inner tetrahedron;
  * singular tails: nine sub-lines below d_min weighted by the wedge exponents gamma = pi / Phi
    (Phi = 2 pi - arccos(1/3), irrational), and radial nine-point tails on the central ray and
    three skeleton-wall rays per edge at each vertex, weighted by the spherical corner exponents
    nu = (0.32221, 1.09227) and tau = 0.907 (Section 5.6 and Appendix A).

The three half-turns about the coordinate axes form D_2 (four parity classes, each solved on one
equilateral face), and the one-third turn about (1, 1, 1) leaves at most two matrices to factorise.
"""
import os
import time

import numpy as np

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
        self.F = np.ascontiguousarray(np.vstack((t_hat, b_hat, np.cross(t_hat, b_hat))))
        self.r_act = np.array([r for r in range(4) if REFL[r, (ax + 1) % 3] > 0], np.int64)
        self.lam = np.log((1.0 - 1.5 * d) / (1.5 * d))          # |t| <= 1/2 - 1.5 d stays inside the skeleton wall (1.5 > sqrt(2))
        r = 0.5 ** np.arange(1, 10) if tail else np.ones(1)
        mid0 = np.zeros(3)
        mid0[ax] = A_TET
        self.mids = np.ascontiguousarray(mid0[None, :] + (d * r)[:, None] * b_hat[None, :])
        if tail:
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
    MFS for the regular tetrahedron of unit edge length at level J (smallest depth 0.14 · 2^{-J}).
      q      : depth ratio, 2^{-1/2} (two per octave) or 1/2
      tails  : with or without the edge and corner singular tails
      J=None : plain inner-tetrahedron layout (n_base barycentric order per face) that stalls
    """
    def __init__(self, J, k=2 * np.pi, q=2 ** -0.5, extra_modes=None, tails=True, classes=CLASSES, n_base=12):
        self.k, self.classes, self.tails = k, [tuple(c) for c in classes], bool(tails and J is not None)
        s = max(0.0, k / (2 * np.pi) - 1.0)
        m = max(1, int(round(np.log(0.5) / np.log(q))))
        if J is None:
            self.depths, ex = np.zeros(0), 0.0
        else:
            self.depths = 0.14 * q ** np.arange(m * J + 1)
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

        # Graded boundary points on F_0 (|dx/d(sc)| = 1 along each half-edge, so H/4 matches the cube);
        # since nu_1 = 0.32221 < 1/2, weighting by dist · r_v^{1 - 2 nu_1} keeps |E|^2 bounded at vertices too
        ng = max(4, 2 + max(m, 1 + int(np.ceil(0.5 * s))))
        d_min = self.depths[-1] if len(self.depths) else 0.14
        levels = max(2, int(np.ceil(np.log2(2 * H / d_min))))
        sc, w = graded_nodes(levels, ng, min(H / 4, np.pi / (2 * k)))
        self.x, W, dist = face_kites(sc, w)
        self.T = tangent_basis(np.tile(-v_hat, (len(self.x), 1)))
        rv = np.min([np.linalg.norm(self.x - v, axis=1) for v in FACE0], axis=0)
        self.rw = np.sqrt(W * dist * (rv ** (1.0 - 2.0 * NU[0])))
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
                lines.line_add(x, T, REFL, ln.r_act, ln.t_hat, ln.mids, ln.lam, ln.d, ln.n, self.k, H, A,
                               cmap[:, lo:lo + 6 * ln.n], sgn, rot=ln.F, w=w)
        if self.tails:
            lo = 3 * self.slot0[-2]
            ne_cols = len(self.ray_pts) * len(self.we) * 3
            lines.radial_tails_add(x, T, REFL, self.ray_pts, self.we, self.wm, self.k, A,
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

    def field(self, coef, pts, batch=1000):
        pts = np.ascontiguousarray(pts, dtype=float)
        if not len(pts) or not coef:
            return np.zeros((len(pts), 3), complex)
        cl = list(coef)
        ncls = len(cl)
        c_all = np.zeros((ncls, 3 * self.slot0[-1]), complex)
        for j, c in enumerate(cl):
            c_all[j, self.keep[c]] = coef[c]
        sgn = np.array([[sign(c, R) for R in REFL] for c in cl])
        y_pt, alpha, beta = self._unpack_points(cl, c_all)
        out = np.zeros((ncls, len(pts), 3), complex)
        lines.eval_points_field(pts, REFL, sgn, y_pt, alpha, beta, float(self.k), out)
        for i, ln in enumerate(self.lines):
            if ln.kind == "e":
                lo = 3 * self.slot0[1 + i]
                c_pair = c_all[:, lo:lo + 6 * ln.n].reshape(ncls, 2, 1, ln.n, 3)
                w = np.stack((ln.w, self.lines[i + 1].w))
                lines.eval_line_field(pts, REFL, ln.r_act, sgn, ln.t_hat, ln.mids, ln.lam, ln.d, ln.n,
                                      float(self.k), H, c_pair, out, rot=ln.F, w=w)
        return out.sum(axis=0)

    def far(self, coef, dirs):
        dirs = np.ascontiguousarray(dirs, dtype=float)
        if not len(dirs) or not coef:
            return np.zeros((len(dirs), 3), complex)
        cl = list(coef)
        ncls = len(cl)
        c_all = np.zeros((ncls, 3 * self.slot0[-1]), complex)
        for j, c in enumerate(cl):
            c_all[j, self.keep[c]] = coef[c]
        sgn = np.array([[sign(c, R) for R in REFL] for c in cl])
        y_pt, alpha, beta = self._unpack_points(cl, c_all)
        out = np.zeros((ncls, len(dirs), 3), complex)
        lines.eval_far_dipoles(dirs, REFL, np.arange(4, dtype=np.int64), sgn, y_pt, alpha, beta, float(self.k), out)
        by_ax = {0: ([], [], []), 1: ([], [], []), 2: ([], [], [])}
        for i, ln in enumerate(self.lines):
            if ln.kind == "e":
                lo = 3 * self.slot0[1 + i]
                c_pair = c_all[:, lo:lo + 6 * ln.n].reshape(ncls, 2, 1, ln.n, 3)
                w = np.stack((ln.w, self.lines[i + 1].w))
                tq, W, off, *_ = lines.build_tree(ln.lam, ln.d, ln.n, H, self.k)
                pm_n = int(off[1])
                c_cart = 2.0 * np.einsum("ksic,lksmc,cq->lkimq", w, c_pair, ln.F)
                M0 = np.einsum("pm,lkimq->lkpiq", W[:pm_n], c_cart).reshape(ncls, 2, -1, 3)
                ys = (ln.mids[None, :, :] + tq[:pm_n, None, None] * ln.t_hat[None, None, :]).reshape(-1, 3)
                by_ax[ln.ax][0].append(ys)
                by_ax[ln.ax][1].append(M0[:, 0])
                by_ax[ln.ax][2].append(M0[:, 1])
        for ax, (yl, al, bl) in by_ax.items():
            if yl:
                r_act = np.array([ri for ri in range(4) if REFL[ri, (ax + 1) % 3] > 0], np.int64)
                lines.eval_far_dipoles(dirs, REFL, r_act, sgn, np.ascontiguousarray(np.vstack(yl)),
                                       np.ascontiguousarray(np.concatenate(al, axis=1)),
                                       np.ascontiguousarray(np.concatenate(bl, axis=1)), float(self.k), out)
        return out.sum(axis=0)


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
            for Q, coef in bsols:
                E += probe.field(coef, pts @ Q) @ Q.T
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
