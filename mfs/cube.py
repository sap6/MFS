"""
Scattering by the perfectly conducting unit cube [-1/2, 1/2]³ (Section 5 of
the paper). By Theorem 1 the sources lie on the skeleton: the twelve bisector
planes of the edges. The layout has
  * edge lines: on the bisector plane of each edge, lines parallel to the edge
    at the depths d_j = 0.2 q^j down to d_min = 0.2 · 2^{-J}, each carrying an
    electric and a magnetic dipole density (Legendre modes in a stretched variable);
  * vertex points: electric and magnetic dipoles on the body diagonals at the same depths;
  * base points: an 8 × 8 grid of electric dipoles on each face of an inner cube;
  * singular tails: nine more lines below d_min (and nine points below the
    deepest vertex point), tied together with weights from the edge and
    vertex exponents, so that they add one set of coefficients each.

The cube is symmetric under the three mirrors x → -x, y → -y, z → -z, so a
wave splits into eight parity classes, each solved with the sources and
boundary points of one octant, and under the rotation about the body
diagonal, so that only four class matrices are ever needed.
"""
import itertools
import time

import numpy as np

from . import lines
from .core import cross_sections, dyadic, factor, magnetic, sphere_rule

H = 0.5                                                          # half side length
REFL = [np.array(r, float) for r in itertools.product((1, -1), repeat=3)]
CLASSES = [tuple(r) for r in REFL]                               # parities (sx, sy, sz), in the same order
ROT = np.eye(3)[[1, 2, 0]]                                       # a third of a turn about the body diagonal
NU, TAU = 0.45417, 0.8400                                        # smallest electric and magnetic exponents of the corner
TAIL_E = {"e": (2 / 3, 5 / 3), "m": (4 / 3, 7 / 3)}              # edge tail exponents: (along the edge, along the bisector)
RULE = sphere_rule(40, 80)                                       # quadrature on the unit sphere, for the optical theorem


def sign(cls, R):
    """Sign of the class cls under the reflection R"""
    return float(np.prod([s for s, r in zip(cls, R) if r < 0]))


def class_part(x, k, d, p, cls):
    """Part of the wave p e^{ik d·x} with the parities cls, at the points x (M, 3)"""
    return sum(sign(cls, R) * (p * R)[None, :] * np.exp(1j * k * ((x * R) @ d))[:, None] for R in REFL) / 8


def classes_of(k, d, p):
    """The parity classes that the wave (d, p) excites"""
    x = np.random.default_rng(0).uniform(-0.5, 0.5, (40, 3))
    return [c for c in CLASSES if np.abs(class_part(x, k, d, p, c)).max() > 1e-12]


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
            cj = tuple(c[(i + j) % 3] for i in range(3))
            if cj in todo:
                todo.remove(cj)
                orbit[c].append(j)
    return orbit


def graded_nodes(levels, ng, hmax=0.1):
    """Gauss-Legendre nodes on (-H, H): ng per panel, panels halving towards both ends, none longer than hmax"""
    dist = np.concatenate(([0.0], H * 0.5 ** np.arange(levels, -1, -1)))   # distance from the edge
    brk = [0.0]
    for a, b in zip(dist[:-1], dist[1:]):
        n = int(np.ceil((b - a) / hmax))
        brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
    brk = np.array(brk)
    ends = np.concatenate((-H + brk, (H - brk[::-1])[1:]))
    g, w = np.polynomial.legendre.leggauss(ng)
    s = np.concatenate([0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(ends[:-1], ends[1:])])
    ws = np.concatenate([0.5 * (b - a) * w for a, b in zip(ends[:-1], ends[1:])])
    return s, ws


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
        if tail:
            # Local frame (along the edge, along the bisector, normal to it) and weights r^e
            self.rot = np.zeros((3, 3))
            self.rot[0, ax] = 1.0
            self.rot[1, self.b] = self.rot[1, self.c] = 1 / np.sqrt(2)
            self.rot[2, self.b], self.rot[2, self.c] = 1 / np.sqrt(2), -1 / np.sqrt(2)
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
      extra_modes: modes added to every line (4 at two wavelengths, 12 at three)
      tails      : with or without the singular tails
      J = None gives the usual layout that fails: only the base grid, n_base × n_base per face
    """
    def __init__(self, J, k=2 * np.pi, q=2 ** -0.5, extra_modes=0, tails=True, classes=CLASSES, n_base=8):
        self.k, self.classes = k, [tuple(c) for c in classes]
        if J is None:
            self.depths, tails, q = np.zeros(0), False, 0.5
        else:
            nd = int(round(J * np.log(0.5) / np.log(q)))         # number of depths below the first
            self.depths = 0.2 * q ** np.arange(nd + 1)
            n0 = max(4, 2 * J - 4) + extra_modes

        # Sources in the octant x, y, z > 0
        s = (np.arange(n_base) + 0.5) / n_base * 0.6 - 0.3       # base grid on the faces of the inner cube, side 0.6
        s = s[s > 0]
        U, V = [g.ravel() for g in np.meshgrid(s, s, indexing="ij")]
        base = np.concatenate([np.insert(np.column_stack((U, V)), ax, 0.3, axis=1) for ax in range(3)])
        vertex = H - self.depths[:, None] / np.sqrt(3) * np.ones((1, 3))
        self.pe, self.pm = np.concatenate((base, vertex)), vertex
        self.lines = [Line(ax, d, n0 + int(np.ceil(2 * np.log(2 * H / d))), kind)   # modes: n0 + ⌈2 ln(2H/d)⌉ (Appendix B)
                      for d in self.depths for ax in range(3) for kind in "em"]
        if tails:                                                # the six lines at d_min (three axes, e and m) continue as tails
            self.lines += [Line(ln.ax, ln.d, ln.n, ln.kind, tail=True) for ln in self.lines[-6:]]
            self.lines += [VertexTail(self.depths[-1], kind) for kind in "em"]

        # Columns: three per slot (point, or mode of a line); a class keeps some of them
        sizes = [len(self.pe) + len(self.pm)] + [ln.n for ln in self.lines]
        self.slot0 = np.cumsum([0] + sizes)
        self.keep = {c: np.concatenate([np.ones((sizes[0], 3), bool)] + [
            ln.keep(c) if isinstance(ln, Line) else np.ones((1, 3), bool) for ln in self.lines]).ravel()
            for c in CLASSES}
        self.unknowns = {c: int(kp.sum()) for c, kp in self.keep.items()}

        # Boundary points in the octant: graded tensor grids, ng nodes per panel,
        # rows weighted by sqrt(area) · (distance to the edge)^{1/2}
        ng = 4 if q > 0.6 else 3
        d_min = self.depths[-1] if len(self.depths) else 0.2
        levels = max(2, int(np.ceil(np.log(d_min / (2 * H)) / np.log(0.5))) + 1)
        s, w = graded_nodes(levels, ng)
        U, V = [g.ravel() for g in np.meshgrid(s[s > 0], s[s > 0], indexing="ij")]
        W = np.outer(w[s > 0], w[s > 0]).ravel()
        self.x = np.concatenate([np.insert(np.column_stack((U, V)), ax, H, axis=1) for ax in range(3)])
        self.T = np.concatenate([np.repeat(np.delete(np.eye(3), ax, axis=0)[None], len(U), axis=0) for ax in range(3)])
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
        self.fac = {c: factor(np.ascontiguousarray(A[i][:, :self.unknowns[c]])) for i, c in enumerate(self.classes)}
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
        A = np.zeros((len(classes), len(x) * nc, nmax), complex)
        sgn = np.array([[sign(c, R) for R in REFL] for c in classes])
        refl, x = np.array(REFL), np.ascontiguousarray(x)

        # Point dipoles: electric, then magnetic
        lines.points_add(x, T, refl, self.pe, self.pm, self.k, A, np.ascontiguousarray(cmap[:, 0]), sgn)

        # Lines, each electric line with its magnetic twin, and the vertex tails
        for i, ln in enumerate(self.lines):
            lo = 3 * self.slot0[1 + i]
            if isinstance(ln, VertexTail):
                for r, R in enumerate(REFL):
                    B = ln.block(self.k, x * R, np.broadcast_to(T * R, (len(x), nc, 3)))
                    for j in range(len(classes)):
                        A[j][:, cmap[j, lo]:cmap[j, lo] + 3] += sgn[j, r] * B
            elif ln.kind == "e":                                 # the magnetic twin follows in the next 3 n columns
                w = None if ln.w is None else np.stack((ln.w, self.lines[i + 1].w))
                lines.tree_add(x, T, refl, ln.ax, ln.b, ln.c, ln.pos, ln.lam, ln.d, ln.n, self.k, H, A,
                               cmap[:, lo:lo + 6 * ln.n], sgn, rot=ln.rot, w=w)
        return A

    def solve(self, d, p, only=None):
        """Coefficients of each class for the wave (d, p), the relative boundary residual, and |b|²"""
        coef, r2, b2 = {}, 0.0, 0.0
        for c in classes_of(self.k, d, p):
            if only is not None and c not in only:
                continue
            b = -np.einsum("mab,mb->ma", self.T, class_part(self.x, self.k, d, p, c)).reshape(-1) * np.repeat(self.rw, 2)
            coef[c], res = self.fac[c].solve(b, residual=True)
            r2, b2 = r2 + float(res) ** 2, b2 + np.linalg.norm(b) ** 2
        return coef, np.sqrt(r2 / b2), b2

    def field(self, coef, pts, batch=1000):
        """Scattered field at the points pts (P, 3), in batches"""
        cl, out = list(coef), []
        for i in range(0, len(pts), batch):
            A = self.matrix(pts[i:i + batch], None, cl)
            out.append(sum(A[j][:, :self.unknowns[c]] @ coef[c] for j, c in enumerate(cl)).reshape(-1, 3))
        return np.concatenate(out)

    def far(self, coef, dirs, rfar=1.0e6):
        """Far-field pattern F, with E ~ F e^{ikr}/r"""
        return self.field(coef, rfar * dirs) * rfar * np.exp(-1j * self.k * rfar)


def solve(J, d=(0, 0, -1.0), p=(0, 1.0, 0), pts=None, k=2 * np.pi, max_gb=6.0, **layout):
    """
    One plane wave at level J, using the mirror planes and the rotation about the body diagonal
      pts: points where the scattered field is wanted
      max_gb: if the class matrices would exceed this, they are built in groups
      Returns the field at pts, the wall time, and a dictionary of checks
      (unknowns, boundary residual, extinction cross-section, optical-theorem mismatch)
    """
    d, p = np.asarray(d, float), np.asarray(p, float)
    t0 = time.time()
    orbit = orbits(classes_of(k, d, p))
    todo, E, sols, info = list(orbit), 0, [], dict(unknowns=[], t_assemble=0.0, t_factor=0.0)
    while todo:
        g = len(todo)
        if g > 1:                                                # how many class matrices fit into max_gb
            probe = Cube(J, k, classes=[], **layout)
            g = max(1, min(g, int(max_gb * 2 ** 30 / (16.0 * probe.rows * probe.unknowns[todo[0]]))))
        cube = Cube(J, k, classes=todo[:g], **layout)
        todo = todo[g:]
        for j in range(3):                                       # the wave turned back by j thirds of a turn
            only = [c for c in cube.classes if j in orbit[c]]
            if only:
                Q = np.linalg.matrix_power(ROT, j)
                coef, res, b2 = cube.solve(Q.T @ d, Q.T @ p, only)
                if pts is not None:
                    E = E + cube.field(coef, pts @ Q) @ Q.T
                sols.append((cube, Q, coef, res, b2))
        info["unknowns"] += [cube.unknowns[c] for c in cube.classes]
        info["t_assemble"] += cube.t_assemble
        info["t_factor"] += cube.t_factor
        cube.fac = None                                          # the checks below need only the layout: free the factors
    t = time.time() - t0

    # Checks, not timed: boundary residual and the optical theorem
    far = lambda dirs: sum(cube.far(coef, dirs @ Q) @ Q.T for cube, Q, coef, _, _ in sols)
    info["ext"], info["sca"] = cross_sections(far, k, d, p, RULE)
    info["residual"] = float(np.sqrt(sum(r ** 2 * b2 for *_, r, b2 in sols) / sum(s[-1] for s in sols)))
    info["optical"] = abs(info["ext"] - info["sca"]) / info["sca"]
    info["rows"] = cube.rows
    return E, t, info
