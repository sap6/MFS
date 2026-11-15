"""MFS for a perfectly conducting unit cube with edge-adapted sources.

Compared with pec_cube_lightning.py (point dipoles clustered at the edges),
the sources near each edge are *line* dipole distributions parallel to the
edge, at depths d_j = depth0 * q**j below it. Each line carries a smooth
(piecewise polynomial) dipole density, on panels graded towards the vertices.
The number of unknowns per layer then grows like log(1/d_j) instead of 1/d_j.

Least-squares rows are weighted by sqrt(area) * dist_to_edge**beta, so that the
(integrable) edge singularity does not dominate the fit.

Accuracy is measured against a Bempp reference on r = 2 (bem_cube_h*.npz) and
by the change between successive refinement levels.
"""
import glob
import itertools
import json
import sys
import time

import numpy as np
from scipy.linalg import lstsq

from pec_cube_lightning import H, K, base_sources, cube_surface, incident
from pec_sphere_convergence import dyadic, fibonacci_sphere

GL = np.polynomial.legendre


def edge_dist(x):
    a = np.sort(np.abs(x), axis=1)
    return H - a[:, 1]


def line_panels(d, pb, eta):
    """Quadrature nodes t (along the edge) and density basis matrix B (nq, nb)
    for a line at depth d. Panels are graded dyadically towards both ends."""
    lo = d                                     # stay a distance d from the vertex
    brk = [lo]
    while brk[-1] * 2 < H * 0.75:
        brk.append(brk[-1] * 2)
    brk.append(H)                              # distances from the vertex, half edge
    brk = np.array(brk)
    half = -H + brk                            # coordinates on [-H, 0]
    edges = np.concatenate((half, -half[-2::-1]))
    ts, cols = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        nq = max(pb + 2, int(np.ceil((b - a) / (eta * d))))
        g, w = GL.leggauss(nq)
        ts.append(0.5 * (a + b) + 0.5 * (b - a) * g)
        cols.append(np.stack([w * GL.Legendre.basis(m)(g) for m in range(pb)], axis=1)
                    * 0.5 * (b - a) / d)       # scaled so that entries stay O(1)
    nq_tot = sum(len(t) for t in ts)
    B = np.zeros((nq_tot, pb * len(ts)))
    r = 0
    for i, c in enumerate(cols):
        B[r:r + len(c), i * pb:(i + 1) * pb] = c
        r += len(c)
    return np.concatenate(ts), B


def edge_lines(levels, depth0, q, pb, eta):
    """List of (quadrature points (nq,3), basis matrix B) for all edges/layers."""
    out = []
    for j in range(levels + 1):
        d = depth0 * q ** j
        t, B = line_panels(d, pb, eta)
        off = H - d / np.sqrt(2)
        for ax in range(3):
            b, c = [i for i in range(3) if i != ax]
            for s1, s2 in itertools.product((-1, 1), repeat=2):
                y = np.zeros((len(t), 3))
                y[:, ax], y[:, b], y[:, c] = t, s1 * off, s2 * off
                out.append((y, B))
    return out


def vertex_points(levels, depth0, q):
    return np.array([np.array(sg) * (H - depth0 * q ** j / np.sqrt(3))
                     for j in range(levels + 1)
                     for sg in itertools.product((-1, 1), repeat=3)])


class Sources:
    def __init__(self, n_base, levels, depth0=0.2, q=0.5, pb=3, eta=0.7):
        pts = [base_sources(n_base, depth0)]
        if levels >= 0:
            pts.append(vertex_points(levels, depth0, q))
        self.points = np.concatenate(pts)
        self.lines = edge_lines(levels, depth0, q, pb, eta) if levels >= 0 else []
        self.n_unknowns = 3 * (len(self.points) + sum(B.shape[1] for _, B in self.lines))

    def matrix(self, x, T=None, budget=2.5e6):
        """Field matrix at points x. With tangents T (M,2,3): rows are the two
        tangential components (2M rows); otherwise all three (3M rows)."""
        nc = 3 if T is None else 2
        A = np.empty((len(x), nc, self.n_unknowns // 3, 3), complex)

        def fill(lo, y, B):
            step = max(1, int(budget / len(y)))
            for i in range(0, len(x), step):
                G = dyadic(x[i:i + step], y, K)
                if T is not None:
                    G = np.einsum("mab,mnbc->manc", T[i:i + step], G)
                else:
                    G = G.transpose(0, 2, 1, 3)
                if B is not None:
                    G = np.einsum("manc,nk->makc", G, B)
                A[i:i + step, :, lo:lo + G.shape[2]] = G
            return lo + (len(y) if B is None else B.shape[1])

        lo = fill(0, self.points, None)
        for y, B in self.lines:
            lo = fill(lo, y, B)
        return A.reshape(len(x) * nc, self.n_unknowns)


def solve(src, col_levels, col_p, beta, cond=1e-12):
    x, w, T = cube_surface(col_levels, col_p)
    rw = np.sqrt(w) * edge_dist(x) ** beta
    A = src.matrix(x, T) * np.repeat(rw, 2)[:, None]
    b = (-np.einsum("mab,mb->ma", T, incident(x)) * rw[:, None]).reshape(-1)
    scale = np.linalg.norm(A, axis=0)
    A /= scale
    t0 = time.time()
    c, _, rank, _ = lstsq(A, b, lapack_driver="gelsy", cond=cond, overwrite_a=True, overwrite_b=True)
    return c / scale, {"rows": len(b), "unknowns": src.n_unknowns, "rank": int(rank),
                       "t_solve": time.time() - t0, "coef_norm": float(np.linalg.norm(c))}


class Test:
    def __init__(self, levels=10, p=3, beta=1 / 3):
        self.x, self.w, self.T = cube_surface(levels, p)
        self.rhs = -np.einsum("mab,mb->ma", self.T, incident(self.x))
        self.dist = edge_dist(self.x)
        self.beta = beta
        files = sorted(glob.glob("bem_cube_h*.npz"), key=lambda f: float(f.split("_h")[1][:-4]))
        self.bem = np.load(files[0]) if files else None
        self.pts = self.bem["pts"] if files else 2.0 * fibonacci_sphere(400)

    def evaluate(self, src, c):
        res = (src.matrix(self.x, self.T) @ c).reshape(-1, 2) - self.rhs
        m2 = np.sum(np.abs(res) ** 2, axis=1)
        r2 = np.sum(np.abs(self.rhs) ** 2, axis=1)
        wb = self.w * self.dist ** (2 * self.beta)
        E = (src.matrix(self.pts) @ c).reshape(-1, 3)
        out = {"bc_l2": float(np.sqrt(np.sum(self.w * m2) / np.sum(self.w * r2))),
               "bc_weighted": float(np.sqrt(np.sum(wb * m2) / np.sum(wb * r2)))}
        if self.bem is not None:
            out["vs_bem"] = float(np.linalg.norm(E - self.bem["E"]) / np.linalg.norm(self.bem["E"]))
        return out, E


def run(configs, beta, tag):
    test = Test(beta=beta)
    results, prev = [], None
    for cfg in configs:
        src = Sources(cfg["n_base"], cfg["levels"], pb=cfg.get("pb", 3), eta=cfg.get("eta", 0.7),
                      q=cfg.get("q", 0.5), depth0=cfg.get("depth0", 0.2))
        c, info = solve(src, cfg["levels"] + cfg.get("col_extra", 2), cfg.get("col_p", 3), beta)
        err, E = test.evaluate(src, c)
        info.update(err, **cfg, beta=beta)
        if prev is not None:
            info["change"] = float(np.linalg.norm(E - prev) / np.linalg.norm(E))
        prev = E
        results.append(info)
        print(f"L={cfg['levels']:2d} unk={info['unknowns']:5d} rows={info['rows']:6d} rank={info['rank']:5d} "
              f"bc_l2={info['bc_l2']:.2e} bc_w={info['bc_weighted']:.2e} "
              f"vs_bem={info.get('vs_bem', float('nan')):.2e} change={info.get('change', float('nan')):.2e} "
              f"t={info['t_solve']:.0f}s", flush=True)
    json.dump(results, open(f"results_pec_cube_edge_{tag}.json", "w"), indent=1)
    return results


if __name__ == "__main__":
    beta = float(sys.argv[1]) if len(sys.argv) > 1 else 1 / 3
    lmax = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    run([{"n_base": 6, "levels": L} for L in range(-1, lmax + 1)], beta, f"beta{beta:.2f}")
