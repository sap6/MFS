"""
Scattering by the perfectly conducting unit cube [-1/2, 1/2]³ (Section 5 of the paper).
By Theorem 1 the sources lie on the skeleton, the twelve bisector planes of the edges:
  * edge lines parallel to each edge at the depths d_j = 0.2 q^j down to d_min = 0.2 · 2^{-J},
    each with an electric and a magnetic dipole density (Legendre modes in a stretched variable);
  * vertex points: electric and magnetic dipoles on the body diagonals at the same depths;
  * base points: an 8 × 8 grid of electric dipoles on each face of an inner cube;
  * singular tails: nine more lines below d_min (and nine points below the deepest vertex
    point), tied with weights from the edge and vertex exponents.
The three mirror planes split a wave into eight parity classes, each solved on one octant,
and the rotation about the body diagonal means that only four class matrices are ever needed.
"""
import itertools
import os
import time

import numpy as np

from . import lines
from .core import class_part, classes_of, cross_sections, dyadic, factor, magnetic, sign, sphere_rule

H = 0.5                                                          # half side length
CLASSES = list(itertools.product((1, -1), repeat=3))             # parities (sx, sy, sz) under the three mirrors
REFL = np.array(CLASSES, float)                                  # the eight reflections, in the same order
ROT = np.eye(3)[[1, 2, 0]]                                       # a third of a turn about the body diagonal
NU, TAU = 0.45417, 0.840                                         # smallest electric and magnetic exponents of the corner
TAIL_E = {"e": (2 / 3, 5 / 3), "m": (4 / 3, 7 / 3)}              # edge tail exponents: (along the edge, along the bisector)
RULE = sphere_rule(40, 80)                                       # quadrature on the unit sphere, for the optical theorem


def orbits(classes):
    """
    Turning a solution by ROT turns its parities (c_x, c_y, c_z) into (c_y, c_z, c_x),
    so the part of parity c turned j times comes from the matrix of class c with
    the wave turned back j times. Returns {class to assemble: [j, ...]}
    """
    orbit, todo = {}, list(classes)
    while todo:
        c = todo.pop(0)
        orbit[c] = [0]
        for j in (1, 2):
            cj = c[j:] + c[:j]
            if cj in todo:
                todo.remove(cj)
                orbit[c].append(j)
    return orbit


def graded_nodes(levels, ng, hmax=H / 4):
    """Gauss-Legendre nodes on (0, H): ng per panel, panels halving towards H, none longer than hmax"""
    dist = np.concatenate(([0.0], H * 0.5 ** np.arange(levels, -1, -1)))   # distance from the edge
    brk = [0.0]
    for a, b in zip(dist[:-1], dist[1:]):
        n = int(np.ceil((b - a) / hmax))
        brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
    ends = H - np.array(brk)[::-1]
    a, b = ends[:-1, None], ends[1:, None]
    g, w = np.polynomial.legendre.leggauss(ng)
    return (0.5 * (a + b) + 0.5 * (b - a) * g).ravel(), (0.5 * (b - a) * w).ravel()


def graded_uv(levels, ng, n_max=0, hmax=H / 4):
    """
    Anisotropic Gauss-Legendre rule on (0, H)²: panels halve towards u = H and v = H,
    with ng × ng nodes near the corner (|p - q| <= 1) and ng across the nearer edge times
    max(3, ng - 1) along it away from the corner (|p - q| >= 2)
    """
    dist = np.concatenate(([0.0], H * 0.5 ** np.arange(levels, -1, -1)))
    brk = [0.0]
    for a, b in zip(dist[:-1], dist[1:]):
        n = int(np.ceil((b - a) / hmax))
        brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
    ends = H - np.array(brk)[::-1]
    a, b = ends[:-1], ends[1:]
    P = len(a)
    ng_t = max(3, ng - 1)
    rules = {n: np.polynomial.legendre.leggauss(n) for n in {ng, ng_t}}
    us, vs, ws = [], [], []
    for p in range(P):
        for q in range(P):
            if abs(p - q) <= 1:
                nu, nv = ng, ng
            elif q > p:
                nu, nv = ng_t, ng
            else:
                nu, nv = ng, ng_t
            gu, wu = rules[nu]
            gv, wv = rules[nv]
            up = 0.5 * (a[p] + b[p]) + 0.5 * (b[p] - a[p]) * gu
            vq = 0.5 * (a[q] + b[q]) + 0.5 * (b[q] - a[q]) * gv
            wp = 0.5 * (b[p] - a[p]) * wu
            wq = 0.5 * (b[q] - a[q]) * wv
            U, V = np.meshgrid(up, vq, indexing="ij")
            us.append(U.ravel())
            vs.append(V.ravel())
            ws.append(np.outer(wp, wq).ravel())
    return np.concatenate(us), np.concatenate(vs), np.concatenate(ws)


class Line:
    """
    Dipole density on a line parallel to the axis ax at depth d below the edge
    (b, c) = (+H, +H), or, with tail=True, the nine tied lines below it
    """
    def __init__(self, ax, d, n, kind, tail=False):
        self.ax, self.d, self.n, self.kind, self.tail = ax, d, n, kind, tail
        self.b, self.c = [i for i in range(3) if i != ax]
        self.lam = np.log((2 * H - d) / d)
        r = 0.5 ** np.arange(1, 10) if tail else np.ones(1)      # depths relative to d
        self.pos = np.column_stack((H - d * r / np.sqrt(2), H - d * r / np.sqrt(2)))
        self.mids = np.zeros((len(r), 3))
        self.mids[:, self.b], self.mids[:, self.c] = self.pos[:, 0], self.pos[:, 1]
        self.t_hat = np.zeros(3)
        self.t_hat[ax] = 1.0
        self.r_act = np.array([ri for ri in range(8) if REFL[ri, ax] > 0], np.int64)
        if tail:
            # Local frame (along the edge, along the bisector, normal to it) and weights r^e
            self.rot = np.zeros((3, 3))
            self.rot[0, ax] = 1.0
            self.rot[1:, [self.b, self.c]] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
            e = np.array(TAIL_E[kind] + (2.0,))                  # the normal component is not used
            self.w = r[None, :, None] ** e[None, None, :]
            self.w[:, -1, :] /= 1 - 0.5 ** e                     # the last line carries the rest of the sheet
        else:
            self.rot, self.w = None, None

    def keep(self, cls):
        """Which (mode, component) columns a parity class keeps, shape (n, 3)"""
        par = (-1.0) ** np.arange(self.n)[:, None]              # parity of the mode along the edge
        axis = 0 if self.tail else self.ax                       # index of the component along the edge
        sc = np.where(np.arange(3) == axis, -1.0, 1.0)[None, :] * (-1.0 if self.kind == "m" else 1.0)
        kept = cls[self.ax] * par * sc > 0
        if self.tail:
            kept[:, 2] = False                                   # a surface current has no normal component
        return kept


class VertexTail:
    """Nine dipoles on the body diagonal below the deepest vertex point, tied with weights r^{nu + 2} (or r^{tau + 2})"""
    def __init__(self, d, kind):
        self.kind, self.n = kind, 1
        r = 0.5 ** np.arange(1, 10)
        e = (TAU if kind == "m" else NU) + 2.0
        self.w = r ** e
        self.w[-1] /= 1 - 0.5 ** e
        self.pts = np.ones((1, 3)) * (H - d * r[:, None] / np.sqrt(3.0))

    def block(self, k, x, T):
        """Field of the tied dipoles at x, tested with T: (M nc, 3)"""
        G = (magnetic if self.kind == "m" else dyadic)(x, self.pts, k)  # (M, 9, 3, 3)
        return np.einsum("mab,mtbc,t->mac", T, G, self.w).reshape(len(x) * T.shape[1], 3)


class Cube:
    """
    MFS for the cube at level J (smallest depth 0.2 · 2^{-J}). The given parity
    classes are assembled and factorised once; solve() then takes any wave
    whose classes are among them.
      q          : ratio of the depths, 2^{-1/2} (two per octave, the paper's layout) or 1/2
      extra_modes: extra modes per unit line length (default: 2 s (s + 1) for s = 2 H k/(2π) - 1 > 0)
      tails      : with or without the singular tails
      J = None gives the usual layout that fails: only the base grid, n_base × n_base per face
    """
    def __init__(self, J, k=2 * np.pi, q=2 ** -0.5, extra_modes=None, tails=True, classes=CLASSES, n_base=8):
        self.k, self.classes, self.tails = k, [tuple(c) for c in classes], bool(tails and J is not None)
        s = max(0.0, k * 2 * H / (2 * np.pi) - 1.0)              # wavelengths per side, minus one
        m = max(1, int(round(np.log(0.5) / np.log(q))))          # depths per octave (2 for q = 2^{-1/2}, 1 for q = 1/2)
        if J is None:
            self.depths, depths_v, ex = np.zeros(0), np.zeros(0), 0.0
        else:
            self.depths = 0.2 * q ** np.arange(m * J + 1)
            mv = max(2, m)                                       # vertex ray uses at least two dipoles per octave
            depths_v = 0.2 * 2.0 ** (-np.arange(mv * J + 1) / mv)
            n0 = max(4, 2 * J - 4)
            ex = 2 * s * (s + 1) if extra_modes is None else extra_modes

        # Sources in the octant x, y, z > 0
        g0 = (np.arange(n_base) + 0.5) / n_base * 0.6 - 0.3      # base grid on the faces of the inner cube, side 0.6
        g0 = g0[g0 > 0]
        U, V = [g.ravel() for g in np.meshgrid(g0, g0, indexing="ij")]
        base = np.concatenate([np.insert(np.column_stack((U, V)), ax, 0.3, axis=1) for ax in range(3)])
        vertex = H - depths_v[:, None] / np.sqrt(3) * np.ones((1, 3))
        self.pe, self.pm = np.concatenate((base, vertex)), vertex
        self.lines = [Line(ax, d, n0 + int(round(ex * 2 * (H - d))) + int(np.ceil(2 * np.log(2 * H / d))), kind)
                      for d in self.depths for ax in range(3) for kind in "em"]
        if self.tails:                                           # the six lines at d_min (three axes, e and m) continue as tails
            self.lines += [Line(ln.ax, ln.d, ln.n, ln.kind, tail=True) for ln in self.lines[-6:]]
            self.lines += [VertexTail(self.depths[-1], kind) for kind in "em"]

        # Columns: three per slot (point, or mode of a line); a class keeps some of them
        sizes = [len(self.pe) + len(self.pm)] + [ln.n for ln in self.lines]
        self.slot0 = np.cumsum([0] + sizes)
        self.keep = {c: np.concatenate([np.ones((sizes[0], 3), bool)] + [
            ln.keep(c) if isinstance(ln, Line) else np.ones((1, 3), bool) for ln in self.lines]).ravel()
            for c in CLASSES}
        self.unknowns = {c: int(kp.sum()) for c, kp in self.keep.items()}

        # Boundary points in the octant: anisotropic graded Gauss panels halving towards the edges
        # from min(H/4, λ/4) down to d_min; rows weighted by sqrt(area · dist)
        ng = 2 + max(m, 1 + int(np.ceil(0.5 * s)))
        d_min = self.depths[-1] if len(self.depths) else 0.2
        levels = max(2, int(np.ceil(np.log2(H / d_min))))
        n_max = max((ln.n for ln in self.lines if isinstance(ln, Line)), default=0)
        U, V, W = graded_uv(levels, ng, n_max, min(H / 4, np.pi / (2 * k)))
        self.x = np.concatenate([np.insert(np.column_stack((U, V)), ax, H, axis=1) for ax in range(3)])
        self.T = np.repeat([np.delete(np.eye(3), ax, axis=0) for ax in range(3)], len(W), axis=0)
        dist = H - np.sort(np.abs(self.x), axis=1)[:, 1]        # distance to the nearest edge
        self.rw = np.sqrt(np.tile(W, 3)) * dist ** 0.5
        self.rows = 2 * len(self.x)

        # Assemble and factorise the class matrices
        self.fac, self.t_assemble, self.t_factor = {}, 0.0, 0.0
        if not self.classes:
            return
        t0 = time.time()
        A = self.matrix(self.x, self.T * self.rw[:, None, None], self.classes)
        self.t_assemble, t0 = time.time() - t0, time.time()
        self.fac = {c: factor(A[i, :, :self.unknowns[c]], overwrite=True) for i, c in enumerate(self.classes)}
        self.t_factor = time.time() - t0

    def matrix(self, x, T, classes):
        """
        Matrices of the classes, summed over the eight mirror images of the points
          x: points (M, 3), T: rows to test with (M, nc, 3), or None for the Cartesian components
          Returns an array (classes, M nc, columns), padded to the largest class
        """
        nmax = max(self.unknowns[c] for c in classes)
        cmap = np.full((len(classes), 3 * self.slot0[-1]), -1, np.int64)   # column of each unknown in each class
        for i, c in enumerate(classes):
            cmap[i, self.keep[c]] = np.arange(self.unknowns[c])
        T = np.eye(3)[None] if T is None else np.ascontiguousarray(T)
        nc = T.shape[1]
        A = np.moveaxis(np.zeros((len(x) * nc, nmax, len(classes)), complex, order="F"), 2, 0)   # each A[i] is F-ordered
        sgn = np.array([[sign(c, R) for R in REFL] for c in classes])
        x = np.ascontiguousarray(x)

        # Point dipoles: electric, then magnetic
        lines.points_add(x, T, REFL, self.pe, self.pm, self.k, A, np.ascontiguousarray(cmap[:, 0]), sgn)

        # Lines, each electric line with its magnetic twin, and the vertex tails
        for i, ln in enumerate(self.lines):
            lo = 3 * self.slot0[1 + i]
            if isinstance(ln, Line) and ln.kind == "e":          # the magnetic twin follows in the next 3 n columns
                w = None if ln.w is None else np.stack((ln.w, self.lines[i + 1].w))
                lines.line_add(x, T, REFL, ln.r_act, ln.t_hat, ln.mids, ln.lam, ln.d, ln.n, self.k, H, A,
                               cmap[:, lo:lo + 6 * ln.n], sgn, rot=ln.rot, w=w)
        if self.tails:
            vt_e, vt_m = self.lines[-2], self.lines[-1]
            lo = 3 * self.slot0[-3]
            lines.radial_tails_add(x, T, REFL, np.ascontiguousarray(vt_e.pts[None]),
                                   np.ascontiguousarray(vt_e.w[None]), np.ascontiguousarray(vt_m.w[None]),
                                   self.k, A, np.ascontiguousarray(cmap[:, lo]), np.ascontiguousarray(cmap[:, lo + 3]), sgn)
        return A

    def solve(self, d, p, only=None):
        """Coefficients of each class for the wave (d, p), the relative boundary residual, and |b|²"""
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
        n_vt = 9 if self.tails else 0
        y_pt = np.vstack((self.pe, self.lines[-2].pts)) if n_vt else self.pe
        alpha = np.zeros((ncls, ne + n_vt, 3), complex)
        beta = np.zeros((ncls, ne + n_vt, 3), complex)
        alpha[:, :ne] = c_all[:, :3 * ne].reshape(ncls, ne, 3)
        if nm:
            beta[:, ne - nm:ne] = c_all[:, 3 * ne:3 * (ne + nm)].reshape(ncls, nm, 3)
        if n_vt:
            vt_e, vt_m = self.lines[-2], self.lines[-1]
            lo = 3 * self.slot0[-3]
            alpha[:, ne:] = c_all[:, lo:lo + 3][:, None, :] * vt_e.w[None, :, None]
            beta[:, ne:] = c_all[:, lo + 3:lo + 6][:, None, :] * vt_m.w[None, :, None]
        return np.ascontiguousarray(y_pt), np.ascontiguousarray(alpha), np.ascontiguousarray(beta)

    def field(self, coef, pts, batch=1000):
        """Scattered field at the points pts (P, 3), evaluated directly from mode coefficients without forming A"""
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
        line_entries = [(i, ln) for i, ln in enumerate(self.lines) if isinstance(ln, Line) and ln.kind == "e"]
        if not line_entries:
            return out.sum(axis=0)

        # Beyond H sqrt(3) + 2 H, the root panel is accepted for every line and mirror image
        far_mask = np.linalg.norm(pts, axis=1) >= H * np.sqrt(3.0) + 2.0 * H + 1e-9
        if np.any(far_mask):
            pts_f = np.ascontiguousarray(pts[far_mask])
            out_f = np.zeros((ncls, len(pts_f), 3), complex)
            by_ax = {0: ([], [], []), 1: ([], [], []), 2: ([], [], [])}
            for i, ln in line_entries:
                lo = 3 * self.slot0[1 + i]
                c_pair = c_all[:, lo:lo + 6 * ln.n].reshape(ncls, 2, 1, ln.n, 3)
                w = np.ones((2, 1, len(ln.mids), 3)) if ln.w is None else np.stack((ln.w, self.lines[i + 1].w))
                rot = np.eye(3) if ln.rot is None else ln.rot
                tq, W, off, *_ = lines.build_tree(ln.lam, ln.d, ln.n, H, self.k)
                pm_n = int(off[1])
                c_cart = 2.0 * np.einsum("ksic,lksmc,cq->lkimq", w, c_pair, rot)
                M0 = np.einsum("pm,lkimq->lkpiq", W[:pm_n], c_cart).reshape(ncls, 2, -1, 3)
                ys = (ln.mids[None, :, :] + tq[:pm_n, None, None] * ln.t_hat[None, None, :]).reshape(-1, 3)
                by_ax[ln.ax][0].append(ys)
                by_ax[ln.ax][1].append(M0[:, 0])
                by_ax[ln.ax][2].append(M0[:, 1])
            for ax, (yl, al, bl) in by_ax.items():
                if yl:
                    r_act = np.array([ri for ri in range(8) if REFL[ri, ax] > 0], np.int64)
                    lines.eval_points_field(pts_f, np.ascontiguousarray(REFL[r_act]), np.ascontiguousarray(sgn[:, r_act]),
                                            np.ascontiguousarray(np.vstack(yl)),
                                            np.ascontiguousarray(np.concatenate(al, axis=1)),
                                            np.ascontiguousarray(np.concatenate(bl, axis=1)), float(self.k), out_f)
            out[:, far_mask] += out_f

        if np.any(~far_mask):
            pts_n = np.ascontiguousarray(pts[~far_mask])
            out_n = np.zeros((ncls, len(pts_n), 3), complex)
            for i, ln in line_entries:
                lo = 3 * self.slot0[1 + i]
                c_pair = c_all[:, lo:lo + 6 * ln.n].reshape(ncls, 2, 1, ln.n, 3)
                w = None if ln.w is None else np.stack((ln.w, self.lines[i + 1].w))
                lines.eval_line_field(pts_n, REFL, ln.r_act, sgn, ln.t_hat, ln.mids, ln.lam, ln.d, ln.n,
                                      float(self.k), H, c_pair, out_n, rot=ln.rot, w=w)
            out[:, ~far_mask] += out_n
        return out.sum(axis=0)

    def far(self, coef, dirs):
        """Exact r -> infinity far-field pattern F, with E ~ F e^{ikr}/r"""
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
        lines.eval_far_dipoles(dirs, REFL, np.arange(8, dtype=np.int64), sgn, y_pt, alpha, beta, float(self.k), out)
        by_ax = {0: ([], [], []), 1: ([], [], []), 2: ([], [], [])}
        for i, ln in enumerate(self.lines):
            if isinstance(ln, Line) and ln.kind == "e":
                lo = 3 * self.slot0[1 + i]
                c_pair = c_all[:, lo:lo + 6 * ln.n].reshape(ncls, 2, 1, ln.n, 3)
                w = np.ones((2, 1, len(ln.mids), 3)) if ln.w is None else np.stack((ln.w, self.lines[i + 1].w))
                rot = np.eye(3) if ln.rot is None else ln.rot
                tq, W, off, *_ = lines.build_tree(ln.lam, ln.d, ln.n, H, self.k)
                pm_n = int(off[1])
                c_cart = 2.0 * np.einsum("ksic,lksmc,cq->lkimq", w, c_pair, rot)
                M0 = np.einsum("pm,lkimq->lkpiq", W[:pm_n], c_cart).reshape(ncls, 2, -1, 3)
                ys = (ln.mids[None, :, :] + tq[:pm_n, None, None] * ln.t_hat[None, None, :]).reshape(-1, 3)
                by_ax[ln.ax][0].append(ys)
                by_ax[ln.ax][1].append(M0[:, 0])
                by_ax[ln.ax][2].append(M0[:, 1])
        for ax, (yl, al, bl) in by_ax.items():
            if yl:
                r_act = np.array([ri for ri in range(8) if REFL[ri, ax] > 0], np.int64)
                lines.eval_far_dipoles(dirs, REFL, r_act, sgn, np.ascontiguousarray(np.vstack(yl)),
                                       np.ascontiguousarray(np.concatenate(al, axis=1)),
                                       np.ascontiguousarray(np.concatenate(bl, axis=1)), float(self.k), out)
        return out.sum(axis=0)


def solve(J, d=(0, 0, -1.0), p=(0, 1.0, 0), pts=None, k=2 * np.pi, max_gb=None, **layout):
    """
    One plane wave at level J, using the mirror planes and the rotation about the body diagonal
      pts: points where the scattered field is wanted
      max_gb: if the class matrices would exceed this, they are built in groups (default: half of RAM)
      Returns the field at pts, the wall time, and a dictionary of checks
      (unknowns, boundary residual, extinction cross-section, optical-theorem mismatch)
    """
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
    probe = Cube(J, k, classes=[], **layout)
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
        for ic, c in enumerate(batch_cls):                       # factorise and solve one class at a time (one R in RAM)
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

    # Checks, not timed: boundary residual and the optical theorem
    far = lambda dirs: sum(probe.far(coef, dirs @ Q) @ Q.T for Q, coef in sols)
    info["ext"], info["sca"] = cross_sections(far, k, d, p, RULE)
    info["residual"] = float(np.sqrt(r2 / b2))
    info["optical"] = abs(info["ext"] - info["sca"]) / info["sca"]
    info["rows"] = probe.rows
    return E, t, info
