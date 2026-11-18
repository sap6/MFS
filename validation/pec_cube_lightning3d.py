"""Lightning-type MFS for scattering by a perfectly conducting unit cube.

Sources
  * a coarse grid of point dipoles on an inner cube (smooth part of the field)
  * for every edge, line dipole distributions parallel to the edge at depths
    d_j = depth0 * q**j, j = 0..J, reaching very close to the edge. The dipole
    density on a line is a polynomial in the stretched variable
        u(t) = log((H + t)/(H - t)) / log((2H - d_j)/d_j),
    which resolves every scale between d_j and the edge length near both
    vertices with O(log(1/d_j)) unknowns.
  * point dipoles on the inward body diagonals at the same depths.

The line integrals are evaluated with a target-adapted quadrature in u: panels
grow dyadically away from the point of the line nearest to the target, so a
line at depth 1e-6 costs a few hundred kernel evaluations per target.

  * below the deepest free depth, the "singular tail" (TailLine): further
    lines tied together by the power laws of the edge expansion.

All line and vertex sources lie on the bisector planes of the edges, i.e. on
the skeleton of the cube, which is where the continued field has its sources
(DERIVATIONS.md, section A).

Rows are weighted by sqrt(area) * dist_to_edge**beta.

This module is the library: the source classes (Line, TailLine, VertexTail,
Sources), their kernels and the graded boundary points. The solver is
cube_mfs.py, which uses the three mirror planes of the cube and solves any
incident plane wave.
"""
import itertools

import numpy as np

import mie
from pec_sphere_convergence import dyadic

K = 2 * np.pi                    # cube side = one wavelength
D = np.array([0.0, 0.0, -1.0])   # propagation direction
P = np.array([0.0, 1.0, 0.0])    # polarisation
H = 0.5                          # half side length


def base_sources(n, depth):
    """n x n midpoint grid on each face of an inner cube at `depth` below the faces."""
    h = H - depth
    s = (np.arange(n) + 0.5) / n * 2 * h - h
    U, V = np.meshgrid(s, s, indexing="ij")
    out = []
    for ax, sign in itertools.product(range(3), (-1, 1)):
        b, c = [i for i in range(3) if i != ax]
        y = np.zeros((U.size, 3))
        y[:, ax], y[:, b], y[:, c] = sign * h, U.ravel(), V.ravel()
        out.append(y)
    return np.concatenate(out)


def incident(pts):
    return mie.plane_wave(pts, K, D, P)

LEG = np.polynomial.legendre

try:
    from line_kernel import line_block
    FAST = True
except ImportError:                       # fall back to the NumPy implementation
    FAST = False
try:                                      # shared-node assembly (line_tree.py); MFS_LINES=old selects line_kernel
    import os
    from line_tree import tree_block
    TREE = FAST and os.environ.get("MFS_LINES", "tree") != "old"
except ImportError:
    TREE = False
_W1 = np.ones((1, 1, 3))


def dyadic_pairs(dvec):
    """Dyadic for separation vectors dvec (..., 3) -> (..., 3, 3)."""
    r = np.linalg.norm(dvec, axis=-1)
    kr = K * r
    g = np.exp(1j * kr) / r
    fI = (1 + 1j / kr - 1 / kr ** 2) * g
    fR = (1 + 3j / kr - 3 / kr ** 2) * g
    rh = dvec / r[..., None]
    return fI[..., None, None] * np.eye(3) - fR[..., None, None] * rh[..., :, None] * rh[..., None, :]


def graded_nodes_capped(levels, p, q, hmax=0.1):
    """Gauss-Legendre nodes on (-H, H): panels graded by ratio q towards both
    ends, with no panel longer than hmax."""
    dist = np.concatenate(([0.0], H * q ** np.arange(levels, -1, -1)))   # distance from the edge
    brk = [0.0]
    for a, b in zip(dist[:-1], dist[1:]):
        n = int(np.ceil((b - a) / hmax))
        brk.extend(a + (b - a) * np.arange(1, n + 1) / n)
    brk = np.array(brk)
    full = np.concatenate((-H + brk, (H - brk[::-1])[1:]))
    g, w = LEG.leggauss(p)
    xs = [0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(full[:-1], full[1:])]
    ws = [0.5 * (b - a) * w for a, b in zip(full[:-1], full[1:])]
    return np.concatenate(xs), np.concatenate(ws)


def magnetic_pairs(dvec):
    """Electric field of magnetic dipoles, curl(g m)/k: (..., 3, 3) acting on m."""
    r = np.linalg.norm(dvec, axis=-1)
    f = (1j - 1 / (K * r)) * np.exp(1j * K * r) / r
    rh = dvec / r[..., None]
    out = np.zeros(dvec.shape[:-1] + (3, 3), complex)
    out[..., 0, 1], out[..., 0, 2] = -rh[..., 2], rh[..., 1]
    out[..., 1, 0], out[..., 1, 2] = rh[..., 2], -rh[..., 0]
    out[..., 2, 0], out[..., 2, 1] = -rh[..., 1], rh[..., 0]
    return f[..., None, None] * out


def cube_surface(levels, p, q, hmax=0.1):
    s, w = graded_nodes_capped(levels, p, q, hmax)
    U, V = np.meshgrid(s, s, indexing="ij")
    W = np.outer(w, w).ravel()
    pts, wts, t1, t2 = [], [], [], []
    for ax, sign in itertools.product(range(3), (-1, 1)):
        b, c = [i for i in range(3) if i != ax]
        x = np.zeros((U.size, 3))
        x[:, ax], x[:, b], x[:, c] = sign * H, U.ravel(), V.ravel()
        pts.append(x)
        wts.append(W)
        t1.append(np.tile(np.eye(3)[b], (U.size, 1)))
        t2.append(np.tile(np.eye(3)[c], (U.size, 1)))
    return (np.concatenate(pts), np.concatenate(wts),
            np.stack((np.concatenate(t1), np.concatenate(t2)), axis=1))


def cube_cells(levels, p, q, hmax=0.1):
    """As cube_surface, plus for every node the index of its quadrature cell
    and the smaller side length of that cell."""
    s, w = graded_nodes_capped(levels, p, q, hmax)
    n1 = len(s) // p
    pid = np.repeat(np.arange(n1), p)                      # 1D panel of each node
    plen = np.repeat(np.add.reduceat(w, np.arange(0, len(s), p)), p)
    I, Jj = np.meshgrid(pid, pid, indexing="ij")
    A, B = np.meshgrid(plen, plen, indexing="ij")
    cell = (I * n1 + Jj).ravel()
    hmin = np.minimum(A, B).ravel()
    x, wts, T = cube_surface(levels, p, q, hmax)
    cells = np.concatenate([cell + f * n1 * n1 for f in range(6)])
    return x, wts, T, cells, np.tile(hmin, 6)


def edge_dist(x):
    return H - np.sort(np.abs(x), axis=1)[:, 1]


class Line:
    """Line dipole distribution parallel to coordinate axis `ax`."""

    def __init__(self, ax, s1, s2, d, n_modes, angle=0.0, kind="e", ng=12, cap=0.17):
        self.ax, self.d, self.n, self.ng, self.cap = ax, d, n_modes, ng, cap
        self.kernel = dyadic_pairs if kind == "e" else magnetic_pairs
        self.b, self.c = [i for i in range(3) if i != ax]
        # position in the cross-section: distance d from the edge, at `angle`
        # (radians) from the inward bisector
        self.yb = s1 * (H - d * np.cos(np.pi / 4 + angle))
        self.yc = s2 * (H - d * np.sin(np.pi / 4 + angle))
        self.te = H - d
        self.lam = np.log((2 * H - d) / d)
        self.g, self.gw = LEG.leggauss(ng)

    def u_of_t(self, t):
        return np.log((H + t) / (H - t)) / self.lam

    def t_of_u(self, u):
        return H * np.tanh(0.5 * u * self.lam)

    def quadrature(self, x):
        """Nodes u (M, nq) and weights (M, nq), adapted to each target."""
        t0 = np.clip(x[:, self.ax], -self.te, self.te)
        r = np.sqrt((x[:, self.b] - self.yb) ** 2 + (x[:, self.c] - self.yc) ** 2
                    + (x[:, self.ax] - t0) ** 2)
        u0 = self.u_of_t(t0)
        w0 = np.minimum(r * (2 * H / (H ** 2 - t0 ** 2)) / self.lam, self.cap)   # width in u
        n_dy = int(np.ceil(np.log2(self.cap / w0.min()))) + 1
        n_cap = int(np.ceil(2 / self.cap))
        # panel offsets from u0, in units that grow dyadically then linearly
        off = np.concatenate(([0.0], 2.0 ** np.arange(n_dy)))[None, :] * w0[:, None]
        off = np.minimum(off, self.cap)
        off = np.concatenate((off, self.cap * (1 + np.arange(1, n_cap + 1))[None, :]
                              * np.ones((len(x), 1))), axis=1)
        off = np.maximum.accumulate(off, axis=1)
        us, ws = [], []
        for sgn in (-1, 1):
            a = np.clip(u0[:, None] + sgn * off[:, :-1], -1, 1)
            b = np.clip(u0[:, None] + sgn * off[:, 1:], -1, 1)
            mid, half = 0.5 * (a + b), 0.5 * (b - a) * sgn        # half >= 0
            us.append((mid[:, :, None] + (half * sgn)[:, :, None] * self.g).reshape(len(x), -1))
            ws.append((half[:, :, None] * self.gw).reshape(len(x), -1))
        return np.concatenate(us, axis=1), np.concatenate(ws, axis=1)

    def block(self, x, T=None):
        """(M, nc, n_modes, 3) field block; nc = 2 with tangents T, else 3.
        Targets are grouped by how close they are to the line, so that only
        nearby targets pay for the deep dyadic refinement."""
        if TREE:
            mag = self.kernel is magnetic_pairs
            Be, Bm = tree_block(x, T, self.ax, self.b, self.c, (self.yb, self.yc), self.lam, self.d, self.n, K, H,
                                we=None if mag else _W1, wm=_W1 if mag else None)
            return Bm if mag else Be
        if FAST:
            out = np.empty((len(x), 3 if T is None else 2, self.n, 3), complex)
            line_block(np.ascontiguousarray(x), np.ascontiguousarray(T) if T is not None else np.zeros((1, 2, 3)),
                       T is not None, self.ax, self.b, self.c, self.yb, self.yc, self.te, self.lam,
                       self.d, self.n, self.g, self.gw, self.cap, self.kernel is magnetic_pairs, K, H, out)
            return out
        t0 = np.clip(x[:, self.ax], -self.te, self.te)
        r = np.sqrt((x[:, self.b] - self.yb) ** 2 + (x[:, self.c] - self.yc) ** 2
                    + (x[:, self.ax] - t0) ** 2)
        w0 = np.minimum(r * (2 * H / (H ** 2 - t0 ** 2)) / self.lam, self.cap)
        key = np.ceil(np.log2(self.cap / w0) / 3).astype(int)          # buckets of 3 dyadic levels
        out = np.empty((len(x), 3 if T is None else 2, self.n, 3), complex)
        for kval in np.unique(key):
            idx = np.nonzero(key == kval)[0]
            for i in range(0, len(idx), 400):
                sub = idx[i:i + 400]
                out[sub] = self._block(x[sub], None if T is None else T[sub])
        return out

    def _block(self, x, T=None):
        u, w = self.quadrature(x)
        t = self.t_of_u(u)
        w = w * H * 0.5 * self.lam / np.cosh(0.5 * u * self.lam) ** 2   # dt = (dt/du) du
        dvec = np.empty(u.shape + (3,))
        dvec[..., self.ax] = x[:, None, self.ax] - t
        dvec[..., self.b] = (x[:, self.b] - self.yb)[:, None]
        dvec[..., self.c] = (x[:, self.c] - self.yc)[:, None]
        G = self.kernel(dvec)                                           # (M, nq, 3, 3)
        if T is not None:
            G = np.einsum("mab,mqbc->mqac", T, G)
        V = LEG.legvander(u, self.n - 1) * (w / self.d)[..., None]      # (M, nq, n)
        return np.einsum("mqk,mqac->makc", V, G)


class TailLine:
    """Singular tail of the line sources of one edge, below the deepest free
    line (depth d). Lines at depths d * qt**i, i = 1..n_tail, all carrying the
    modes of the deepest free line, are tied together with fixed weights
    (d_i / d)**e: one set of modes per exponent instead of one per depth. The
    geometric remainder below the last depth is lumped into the last line.

    Dipole components are taken in the frame of the edge: along the edge
    (axis), along the bisector (radial) and normal to the bisector plane. The
    tail stands for the surface current that the continued field leaves on
    the bisector half-plane (DERIVATIONS.md, section B). A sheet density
    rho^(e - 1) gives line weights (d_i / d)**e, and the edge expansion gives
      electric axis    J_s   = [H_rho]  ~ rho^(-1/3)   e = 2/3   (next 8/3)
      electric radial  J_rho = [H_s]    ~ rho^(2/3)    e = 5/3   (next 11/3)
      magnetic axis    K_s   = [E_rho]  ~ rho^(1/3)    e = 4/3   (next 8/3)
      magnetic radial  K_rho = [E_s]    ~ rho^(4/3)    e = 7/3   (next 11/3)
    A surface current has no component normal to the bisector plane, so the
    normal dipoles are left out unless `normal` is set. One exponent per
    component is enough: the second terms change nothing at J = 5, and the
    integer exponents used in earlier runs (1, 2, 1, 2) were not physical."""

    def __init__(self, base, n_tail, qt, exps_axis, exps_radial, normal=False, taper=False):
        self.base, self.taper = base, bool(taper)
        self.ax, self.yb, self.yc, self.kernel = base.ax, base.yb, base.yc, base.kernel
        ns = max(len(exps_axis), len(exps_radial))
        self.nm, self.n = base.n, base.n * ns
        self.deg = np.tile(np.arange(base.n), ns)
        self.axis_index = 0
        ecomp = np.full((ns, 3), np.nan)
        ecomp[:len(exps_axis), 0] = exps_axis
        ecomp[:len(exps_radial), 1] = exps_radial
        if normal:
            ecomp[:len(exps_radial), 2] = exps_radial
        self.colmask = np.repeat(~np.isnan(ecomp), base.n, axis=0)        # (n, 3)
        ecomp = np.where(np.isnan(ecomp), 2.0, ecomp)
        ratio = qt ** np.arange(1, n_tail + 1)
        self.w = ratio[None, :, None] ** ecomp[:, None, :]              # (slots, n_tail, 3)
        self.w[:, -1, :] /= 1 - qt ** ecomp
        s1, s2 = np.sign(base.yb), np.sign(base.yc)
        self.pos = [(s1 * (H - base.d * r * np.cos(np.pi / 4)), s2 * (H - base.d * r * np.sin(np.pi / 4)))
                    for r in ratio]
        self.rot = np.zeros((3, 3))                                      # local (axis, radial, normal) <- Cartesian
        self.rot[0, base.ax] = 1.0
        self.rot[1, base.b], self.rot[1, base.c] = s1 / np.sqrt(2), s2 / np.sqrt(2)
        self.rot[2, base.b], self.rot[2, base.c] = s1 / np.sqrt(2), -s2 / np.sqrt(2)

    def block(self, x, T=None):
        b = self.base
        if TREE and not self.taper:
            mag = b.kernel is magnetic_pairs
            Be, Bm = tree_block(x, T, b.ax, b.b, b.c, np.array(self.pos), b.lam, b.d, b.n, K, H, rot=self.rot,
                                we=None if mag else self.w, wm=self.w if mag else None)
            return Bm if mag else Be
        ns = self.n // self.nm
        out = np.zeros((len(x), 3 if T is None else 2, self.n, 3), complex)
        tmp = np.empty((len(x), 3 if T is None else 2, self.nm, 3), complex)
        xc = np.ascontiguousarray(x)
        Tc = np.ascontiguousarray(T) if T is not None else np.zeros((1, 2, 3))
        for i, (yb, yc) in enumerate(self.pos):
            line_block(xc, Tc, T is not None, b.ax, b.b, b.c, yb, yc, b.te, b.lam, b.d, b.n, b.g, b.gw,
                       b.cap, b.kernel is magnetic_pairs, K, H, tmp, self.taper)
            loc = tmp @ self.rot.T                                        # source components in the edge frame
            for ie in range(ns):
                out[:, :, ie * self.nm:(ie + 1) * self.nm] += loc * self.w[ie, i]
        return out


class VertexTail:
    """Singular tail of the point dipoles of one vertex: points on the inward
    body diagonal at distances d * qt**i below the deepest free one, tied
    together with weights (d_i / d)**e. Near the vertex the field behaves
    like r^(nu - 1) (electric type, nu = 0.45417) and H like r^(tau - 1)
    (tau = 0.840), see cone_exponents.py and DERIVATIONS.md D; dipoles of
    strength r^e on a geometric row give a field r^(e - 3), so e = nu + 2
    for electric and tau + 2 for magnetic dipoles. Optional (vtail_n > 0)."""

    NU, TAU = (0.45417, 1.23087), (0.8400, 1.20635)

    def __init__(self, sg, d, n_tail, qt, n_exp, magnetic):
        self.ax, self.yb, self.yc = 2, float(sg[0]), float(sg[1])       # stored like a line along z through the vertex
        self.kernel = magnetic_pairs if magnetic else dyadic_pairs
        exps = np.array((self.TAU if magnetic else self.NU)[:n_exp]) + 2.0
        self.n = len(exps)
        ratio = qt ** np.arange(1, n_tail + 1)
        self.w = ratio[:, None] ** exps[None, :]
        self.w[-1] /= 1 - qt ** exps
        self.pts = np.asarray(sg, float)[None, :] * (H - d * ratio[:, None] / np.sqrt(3.0))

    def block(self, x, T=None):
        G = self.kernel(x[:, None, :] - self.pts[None, :, :])            # (m, n_tail, 3, 3)
        B = np.einsum("mtbc,te->mebc", G, self.w)
        return np.einsum("mab,mebc->maec", T, B) if T is not None else B.transpose(0, 2, 1, 3)


def same_place(l1, l2):
    """True if l1 is an electric and l2 a magnetic line (or tail) with the
    same geometry and modes, so that one pass over the kernel serves both."""
    if type(l1) is not type(l2) or type(l1) not in (Line, TailLine):
        return False
    if isinstance(l1, TailLine):
        if l1.taper or l2.taper or l1.pos != l2.pos:
            return False
        l1, l2 = l1.base, l2.base
    return (l1.kernel is dyadic_pairs and l2.kernel is magnetic_pairs
            and (l1.ax, l1.yb, l1.yc, l1.d, l1.n) == (l2.ax, l2.yb, l2.yc, l2.d, l2.n))


def pair_block(l1, l2, x, T=None):
    """Blocks of an electric line (or tail) and of its magnetic twin."""
    if isinstance(l1, TailLine):
        b = l1.base
        return tree_block(x, T, b.ax, b.b, b.c, np.array(l1.pos), b.lam, b.d, b.n, K, H, rot=l1.rot, we=l1.w, wm=l2.w)
    return tree_block(x, T, l1.ax, l1.b, l1.c, (l1.yb, l1.yc), l1.lam, l1.d, l1.n, K, H, we=_W1, wm=_W1)


class Sources:
    def __init__(self, J, n_base=6, depth0=0.2, q=0.25, gamma=2.0, n0=4, angles=(0.0,), kinds="e", vrays=1,
                 tail_n=0, tail_q=0.5, tail_e=(2 / 3,), tail_m=(4 / 3,), tail_normal=0,
                 tail_er=(5 / 3,), tail_mr=(7 / 3,), vtail_n=0, vtail_exp=1, depths=None):
        self.depths = depth0 * q ** np.arange(J + 1) if J >= 0 else np.array([])
        if depths is not None and J >= 0:                    # any decreasing list of depths (cube_mfs.tapered)
            self.depths = np.asarray(depths, float)
        pts = [base_sources(n_base, depth0)]
        # inward rays from every vertex: the body diagonal, optionally three more
        dirs = [np.array([1.0, 1, 1])]
        if vrays > 1:
            dirs += [np.array(v, float) for v in ((2, 1, 1), (1, 2, 1), (1, 1, 2))]
        vpts = []
        for d in self.depths:
            for dr in dirs:
                dr = dr / np.linalg.norm(dr)
                vpts.append(np.array([np.array(sg) * (H - d * dr)
                                      for sg in itertools.product((-1, 1), repeat=3)]))
        vpts = np.concatenate(vpts) if vpts else np.zeros((0, 3))
        self.points = np.concatenate(pts + [vpts])
        self.mpoints = vpts if "m" in kinds else np.zeros((0, 3))
        self.lines = []
        for d in self.depths:
            n_modes = n0 + int(np.ceil(gamma * np.log(2 * H / d)))
            for ax in range(3):
                for s1, s2 in itertools.product((-1, 1), repeat=2):
                    for ang in angles:
                        for kind in kinds:
                            self.lines.append(Line(ax, s1, s2, d, n_modes, np.deg2rad(ang), kind))
        if tail_n and J >= 0:
            deepest = [ln for ln in self.lines if ln.d == self.depths[-1]]
            for ln in deepest:
                ea = tail_m if ln.kernel is magnetic_pairs else tail_e
                er = tuple(e + 1 for e in ea) if tail_normal else (
                    tail_mr if ln.kernel is magnetic_pairs else tail_er)
                self.lines.append(TailLine(ln, tail_n, tail_q, ea, er, tail_normal))
        if vtail_n and J >= 0:
            for sg in itertools.product((-1, 1), repeat=3):
                for kind in kinds:
                    self.lines.append(VertexTail(sg, self.depths[-1], vtail_n, tail_q, vtail_exp, kind == "m"))
        self.n_unknowns = 3 * (len(self.points) + len(self.mpoints) + sum(l.n for l in self.lines))

    def matrix(self, x, T=None, batch=1200):
        nc = 3 if T is None else 2
        A = np.empty((len(x), nc, self.n_unknowns // 3, 3), complex)
        for i in range(0, len(x), batch):
            sl = slice(i, i + batch)
            G = dyadic(x[sl], self.points, K)
            A[sl, :, :len(self.points)] = (np.einsum("mab,mnbc->manc", T[sl], G) if T is not None
                                           else G.transpose(0, 2, 1, 3))
            lo = len(self.points)
            if len(self.mpoints):
                G = magnetic_pairs(x[sl][:, None, :] - self.mpoints[None, :, :])
                A[sl, :, lo:lo + len(self.mpoints)] = (np.einsum("mab,mnbc->manc", T[sl], G) if T is not None
                                                       else G.transpose(0, 2, 1, 3))
                lo += len(self.mpoints)
            i, Tsl = 0, None if T is None else T[sl]
            while i < len(self.lines):
                line = self.lines[i]
                mate = self.lines[i + 1] if i + 1 < len(self.lines) else None
                if TREE and mate is not None and same_place(line, mate):
                    A[sl, :, lo:lo + line.n], A[sl, :, lo + line.n:lo + line.n + mate.n] = pair_block(line, mate, x[sl], Tsl)
                    lo += line.n + mate.n
                    i += 2
                else:
                    A[sl, :, lo:lo + line.n] = line.block(x[sl], Tsl)
                    lo += line.n
                    i += 1
        return A.reshape(len(x) * nc, self.n_unknowns)
