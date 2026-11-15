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

Rows are weighted by sqrt(area) * dist_to_edge**beta.
"""
import glob
import itertools
import json
import sys
import time

import numpy as np
from scipy.linalg import lstsq

from pec_cube_lightning import H, K, base_sources, graded_nodes, incident
from pec_sphere_convergence import dyadic, fibonacci_sphere

LEG = np.polynomial.legendre

try:
    from line_kernel import line_block
    FAST = True
except ImportError:                       # fall back to the NumPy implementation
    FAST = False


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


def cube_surface(levels, p, q):
    s, w = graded_nodes_capped(levels, p, q)
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


def cube_cells(levels, p, q):
    """As cube_surface, plus for every node the index of its quadrature cell
    and the smaller side length of that cell."""
    s, w = graded_nodes_capped(levels, p, q)
    n1 = len(s) // p
    pid = np.repeat(np.arange(n1), p)                      # 1D panel of each node
    plen = np.repeat(np.add.reduceat(w, np.arange(0, len(s), p)), p)
    I, Jj = np.meshgrid(pid, pid, indexing="ij")
    A, B = np.meshgrid(plen, plen, indexing="ij")
    cell = (I * n1 + Jj).ravel()
    hmin = np.minimum(A, B).ravel()
    x, wts, T = cube_surface(levels, p, q)
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
    (d_i / d)**e for each exponent e: one set of modes per exponent instead of
    one per depth. The exponents follow the edge expansion of the field
    (rho^(2n/3) and integer powers): e = 2/3, 1, 4/3, ... for dipoles along
    the edge and e + 1 for transverse dipoles. The geometric remainder below the last depth is lumped into
    the last line."""

    def __init__(self, base, n_tail, qt, exps):
        self.base, self.exps = base, tuple(exps)
        self.ax, self.yb, self.yc, self.kernel = base.ax, base.yb, base.yc, base.kernel
        self.nm = base.n
        self.n = base.n * len(self.exps)
        self.deg = np.tile(np.arange(base.n), len(self.exps))
        ratio = qt ** np.arange(1, n_tail + 1)
        # exponent e for the dipole component along the edge, e + 1 for the
        # two transverse components (their near field is one power more singular)
        ecomp = np.array(self.exps)[:, None] + np.where(np.arange(3) == base.ax, 0.0, 1.0)[None, :]
        self.w = ratio[None, :, None] ** ecomp[:, None, :]              # (n_exp, n_tail, 3)
        self.w[:, -1, :] /= 1 - qt ** ecomp
        s1, s2 = np.sign(base.yb), np.sign(base.yc)
        self.pos = [(s1 * (H - base.d * r * np.cos(np.pi / 4)), s2 * (H - base.d * r * np.sin(np.pi / 4)))
                    for r in ratio]

    def block(self, x, T=None):
        b = self.base
        out = np.zeros((len(x), 3 if T is None else 2, self.n, 3), complex)
        tmp = np.empty((len(x), 3 if T is None else 2, self.nm, 3), complex)
        xc = np.ascontiguousarray(x)
        Tc = np.ascontiguousarray(T) if T is not None else np.zeros((1, 2, 3))
        for i, (yb, yc) in enumerate(self.pos):
            line_block(xc, Tc, T is not None, b.ax, b.b, b.c, yb, yc, b.te, b.lam, b.d, b.n, b.g, b.gw,
                       b.cap, b.kernel is magnetic_pairs, K, H, tmp)
            for ie in range(len(self.exps)):
                out[:, :, ie * self.nm:(ie + 1) * self.nm] += tmp * self.w[ie, i]
        return out


class Sources:
    def __init__(self, J, n_base=6, depth0=0.2, q=0.25, gamma=2.0, n0=4, angles=(0.0,), kinds="e", vrays=1,
                 tail_n=0, tail_q=0.5, tail_exps=(2 / 3, 1.0)):
        self.depths = depth0 * q ** np.arange(J + 1) if J >= 0 else np.array([])
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
            self.lines += [TailLine(ln, tail_n, tail_q, tail_exps) for ln in deepest]
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
            for line in self.lines:
                A[sl, :, lo:lo + line.n] = line.block(x[sl], None if T is None else T[sl])
                lo += line.n
        return A.reshape(len(x) * nc, self.n_unknowns)


class SymSources(Sources):
    """Sources reduced by the two mirror symmetries of the problem
    (E_inc = y_hat exp(-ikz)): even under x -> -x, odd under y -> -y.
    Only a quarter of the sources is kept; each basis function is the
    symmetrised sum over the four images, evaluated by reflecting targets."""

    GROUP = [(np.array([1.0, 1, 1]), 1.0), (np.array([-1.0, 1, 1]), 1.0),
             (np.array([1.0, -1, 1]), -1.0), (np.array([-1.0, -1, 1]), -1.0)]

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.points = self.points[(self.points[:, 0] > 0) & (self.points[:, 1] > 0)]
        self.mpoints = self.mpoints[(self.mpoints[:, 0] > 0) & (self.mpoints[:, 1] > 0)] if len(self.mpoints) else self.mpoints
        lines, keep = [], [np.ones((len(self.points) + len(self.mpoints), 3), bool)]
        for ln in self.lines:
            mag = ln.kernel is magnetic_pairs
            if ln.ax == 2:                       # z-lines: orbit of four, keep x>0, y>0
                if ln.yb > 0 and ln.yc > 0:
                    lines.append(ln)
                    keep.append(np.ones((ln.n, 3), bool))
                continue
            if ln.yb < 0:                        # x-lines: b is y; y-lines: b is x
                continue
            sigma = 1.0 if ln.ax == 0 else -1.0  # parity of the reflection that maps the line to itself
            par = (-1.0) ** getattr(ln, 'deg', np.arange(ln.n))[:, None]
            sc = np.where(np.arange(3) == ln.ax, -1.0, 1.0)[None, :] * (-1.0 if mag else 1.0)
            lines.append(ln)
            keep.append(sigma * par * sc > 0)
        self.lines = lines
        self.keep = np.concatenate(keep).reshape(-1)
        self._n_all = 3 * (len(self.points) + len(self.mpoints) + sum(l.n for l in self.lines))
        self.n_unknowns = int(self.keep.sum())

    def matrix(self, x, T=None, batch=1200):
        out = 0
        n_keep, self.n_unknowns = self.n_unknowns, self._n_all
        try:
            for R, sig in self.GROUP:
                if T is not None:
                    A = Sources.matrix(self, x * R, T * R, batch)
                else:
                    A = Sources.matrix(self, x * R, None, batch).reshape(len(x), 3, -1)
                    A = (A * R[None, :, None]).reshape(len(x) * 3, -1)
                out = out + sig * A[:, self.keep]
        finally:
            self.n_unknowns = n_keep
        return out


def quarter(x):
    return (x[:, 0] > 0) & (x[:, 1] > 0)


def solve(src, col_levels, col_p, col_q, beta, cond=1e-12, moments=0):
    x, w, T, cells, hmin = cube_cells(col_levels, col_p, col_q)
    if isinstance(src, SymSources):
        sel = quarter(x)
        x, w, T, cells, hmin = x[sel], w[sel], T[sel], cells[sel], hmin[sel]
        cells = np.unique(cells, return_inverse=True)[1]
    t0 = time.time()
    A = src.matrix(x, T)
    b = -np.einsum("mab,mb->ma", T, incident(x)).reshape(-1)
    if moments:
        # test against cell averages, scaled to mimic an H^{-1/2} norm:
        # row = sqrt(h_c / A_c) * integral over the cell
        import scipy.sparse as sp
        nc = cells.max() + 1
        area = np.bincount(cells, weights=w, minlength=nc)
        hc = np.zeros(nc)
        hc[cells] = hmin
        S = sp.csr_matrix((w * np.sqrt(hc / area)[cells], (cells, np.arange(len(x)))), shape=(nc, len(x)))
        S2 = sp.kron(S, sp.identity(2), format="csr")
        A, b = S2 @ A, S2 @ b
    else:
        rw = np.repeat(np.sqrt(w) * edge_dist(x) ** beta, 2)
        A *= rw[:, None]
        b = b * rw
    scale = np.linalg.norm(A, axis=0)
    A /= scale
    t1 = time.time()
    c, _, rank, _ = lstsq(A, b, lapack_driver="gelsy", cond=cond, overwrite_a=True, overwrite_b=True)
    return c / scale, {"rows": len(b), "unknowns": src.n_unknowns, "rank": int(rank),
                       "t_assemble": t1 - t0, "t_solve": time.time() - t1,
                       "coef_norm": float(np.linalg.norm(c))}


def load_bem():
    files = sorted(glob.glob("bem_cube_h*.npz"), key=lambda f: float(f.split("_h")[1][:-4]))
    return np.load(files[0]) if files else None


def main(Js, beta=1 / 3, q=0.25, gamma=2.0, col_q=1 / 3, col_p=2, n_base=6, tag="run", angles="0", kinds="e", moments=0,
         col_hmax=0.1, sym=0, vrays=1, cond=1e-12, tail=0, tail_q=0.5, tail_exps="0.6667:1",
         col_extra=1):
    angles = tuple(float(a) for a in str(angles).split(":"))
    texps = tuple(round(3 * float(a)) / 3 for a in str(tail_exps).split(":"))
    bem = load_bem()
    pts = bem["pts"] if bem is not None else 2.0 * fibonacci_sphere(400)
    results, prev = [], None
    for J in Js:
        src = (SymSources if sym else Sources)(J, n_base=n_base, q=q, gamma=gamma, angles=angles, kinds=kinds, vrays=vrays,
                                               tail_n=tail, tail_q=tail_q, tail_exps=texps)
        d_min = src.depths[-1] if J >= 0 else 0.2
        col_levels = max(2, int(np.ceil(np.log(d_min / (2 * H)) / np.log(col_q))) + col_extra)
        c, info = solve(src, col_levels, col_p, col_q, beta, moments=moments, cond=cond)
        E = (src.matrix(pts) @ c).reshape(-1, 3)
        info.update(J=J, d_min=float(d_min), beta=beta, q=q, gamma=gamma, col_q=col_q, col_p=col_p,
                    n_base=n_base)
        if bem is not None:
            info["vs_bem"] = float(np.linalg.norm(E - bem["E"]) / np.linalg.norm(bem["E"]))
        if prev is not None:
            info["change"] = float(np.linalg.norm(E - prev) / np.linalg.norm(E))
        prev = E
        results.append(info)
        print(f"J={J:2d} d_min={d_min:.1e} unk={info['unknowns']:5d} rows={info['rows']:6d} "
              f"rank={info['rank']:5d} vs_bem={info.get('vs_bem', float('nan')):.2e} "
              f"change={info.get('change', float('nan')):.2e} |c|={info['coef_norm']:.1e} "
              f"asm={info['t_assemble']:.0f}s solve={info['t_solve']:.0f}s", flush=True)
        json.dump(results, open(f"results_pec_cube_lightning3d_{tag}.json", "w"), indent=1)
        np.save(f"cube_lightning3d_field_{tag}.npy", E)
    return results


if __name__ == "__main__":
    if sys.argv[1] == "check":
        # line quadrature against brute force
        line = Line(0, 1, 1, 0.01, 12)
        rng = np.random.default_rng(0)
        x = np.column_stack((rng.uniform(-H, H, 40), np.full(40, H), H - 10 ** rng.uniform(-3, -0.4, 40)))
        fast = line.block(x)
        g, w = LEG.leggauss(40)
        edges = np.linspace(-1, 1, 4001)
        u = (0.5 * (edges[:-1] + edges[1:])[:, None] + 0.5 * np.diff(edges)[:, None] * g).ravel()
        wu = (0.5 * np.diff(edges)[:, None] * w).ravel()
        t = line.t_of_u(u)
        wt = wu * H * 0.5 * line.lam / np.cosh(0.5 * u * line.lam) ** 2
        y = np.zeros((len(t), 3))
        y[:, 0], y[:, 1], y[:, 2] = t, line.yb, line.yc
        V = LEG.legvander(u, line.n - 1) * (wt / line.d)[:, None]
        slow = np.einsum("qk,mqac->makc", V, dyadic(x, y, K))
        print("line quadrature max rel error:",
              np.abs(fast - slow).max() / np.abs(slow).max())
    else:
        Js = [int(a) for a in sys.argv[1].split(",")]
        kw = dict(a.split("=") for a in sys.argv[2:])
        main(Js, **{k: (v if k in ("tag", "angles", "kinds", "tail_exps") else float(v) if "." in v or "e" in v else int(v))
                    for k, v in kw.items()})
