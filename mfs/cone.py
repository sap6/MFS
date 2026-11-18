"""
Exponents of the field at a conical point (Appendix A of the paper): the tip
of the almond (an elliptic cone) and the vertex of the cube.

Near the tip of a perfectly conducting cone the field is a sum of terms
  electric type: E = ∇(r^nu Y),   E ~ r^(nu - 1)
  magnetic type: H = ∇(r^tau Z),  H ~ r^(tau - 1)
where Y (Z) is a Dirichlet (Neumann) eigenfunction of the Laplace-Beltrami
operator on the part of the unit sphere outside the cone, with eigenvalue
nu(nu + 1) (tau(tau + 1)). We find nu and tau by the method of particular
solutions (Betcke and Trefethen 2005): the functions P_nu(-x·s), with poles s
inside the cone, solve the eigenvalue equation on the sphere, and nu is an
exponent when a combination of them vanishes on the boundary of the cone
without vanishing inside. This is again a method of fundamental solutions.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import lpmv

T1 = 0.583333                                     # tip of the almond (units of L), and the half-angles of its cone (Woo et al. 1993)
TAN_Y = 4.83345 * T1 / (2.08335 ** 2 * 0.96)
TAN_Z = 1.61115 * T1 / (2.08335 ** 2 * 0.96)


def P(nu, x):
    """Legendre function P_nu(x)"""
    return lpmv(0, nu, x)


def dP(nu, x):
    """P_nu'(x)"""
    return nu * (P(nu - 1, x) - x * P(nu, x)) / (1 - x * x)


def d2P(nu, x):
    """P_nu''(x), from Legendre's equation"""
    return (nu * (nu + 1) * P(nu, x) - 2 * x * dP(nu, x)) / (x * x - 1)


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def smallest_sigma(B, C):
    """Smallest singular value of the boundary part B of an orthonormal basis of the span of [B; C]: zero at an exponent"""
    A = np.vstack((B, C))
    U, s, _ = np.linalg.svd(A / np.linalg.norm(A, axis=0), full_matrices=False)
    return np.linalg.svd(U[:len(B), s > 1e-12 * s[0]], compute_uv=False)[-1]


def minima(sigma, lo, hi, step, tol=None):
    """
    Local minima of sigma(nu) on a grid in (lo, hi), refined (integers are skipped);
    with tol, only those that reach below tol times the median of sigma
    """
    grid = np.arange(lo, hi, step)
    grid = grid[np.abs(grid - np.round(grid)) > 1e-6]
    sig = np.array([sigma(v) for v in grid])
    out = []
    for i in range(1, len(grid) - 1):
        if sig[i] < sig[i - 1] and sig[i] <= sig[i + 1]:
            r = minimize_scalar(sigma, bounds=(grid[i - 1], grid[i + 1]), method="bounded", options={"xatol": 1e-10})
            if tol is None or r.fun < tol * np.median(sig):
                out.append(r.x)
    return out


class EllipticCone:
    """
    Cone along -x with the tip at the origin and cross-section (y/tan_y)² + (z/tan_z)² = x².
    Poles on the focal segment of the cross-section, clustered towards its ends. The class
    (sy, sz) gives the parities of the potential under y → -y and z → -z.
    """
    def __init__(self, tan_y=TAN_Y, tan_z=TAN_Z, n_src=24, n_bnd=400, n_int=60, seed=0):
        psi = (np.arange(n_bnd) + 0.5) / n_bnd * (np.pi / 2)            # one quadrant of the boundary curve
        self.xb = unit(np.column_stack((-np.ones(n_bnd), tan_y * np.cos(psi), tan_z * np.sin(psi))))
        tang = unit(np.column_stack((0 * psi, -tan_y * np.sin(psi), tan_z * np.cos(psi))))
        self.nb = np.cross(self.xb, unit(tang - np.sum(tang * self.xb, axis=1)[:, None] * self.xb))   # conormal
        f = 0.98 * np.sqrt(tan_y ** 2 - tan_z ** 2)
        u = np.cos((np.arange(n_src) + 0.5) / n_src * (np.pi / 2))
        self.src = unit(np.column_stack((-np.ones(n_src), f * u, 0 * u)))
        rng = np.random.default_rng(seed)
        pts = unit(rng.normal(size=(40 * n_int, 3)))
        pts = np.abs(pts)
        pts[:, 0] *= np.sign(rng.normal(size=len(pts)))
        inside = (pts[:, 0] < 0) & ((pts[:, 1] / tan_y) ** 2 + (pts[:, 2] / tan_z) ** 2 < pts[:, 0] ** 2)
        self.xi = pts[~inside][:n_int]

    def columns(self, nu, x, sy, sz, normal=None):
        """Symmetrised particular solutions at x (values, or derivatives along normal); odd in z: the derivative across the segment"""
        ez, cols = np.array([0.0, 0, 1]), []
        for s in self.src:
            col = 0
            for ry, sign in ((1, 1.0), (-1, sy)):
                sj = s * [1.0, ry, 1.0]
                c = np.clip(-(x @ sj), -1 + 1e-15, 1 - 1e-15)
                if sz > 0:
                    col = col + sign * (P(nu, c) if normal is None else -dP(nu, c) * (normal @ sj))
                elif normal is None:
                    col = col + sign * (-dP(nu, c) * (x @ ez))
                else:
                    col = col + sign * (d2P(nu, c) * (normal @ sj) * (x @ ez) - dP(nu, c) * (normal @ ez))
            cols.append(col)
        return np.column_stack(cols)

    def exponents(self, sy, sz, neumann=False, hi=3.6):
        """Exponents nu (Dirichlet) or tau (neumann=True) of the class (sy, sz), below hi"""
        def sigma(nu):
            B = self.columns(nu, self.xb, sy, sz, self.nb if neumann else None)
            return smallest_sigma(B / np.sqrt(len(B)), self.columns(nu, self.xi, sy, sz) / np.sqrt(len(self.xi)))
        return minima(sigma, 0.02, hi + 1e-9, 0.02, tol=1e-3)


class CubeVertex:
    """
    Vertex of the solid octant x, y, z < 0; the field lives in the other seven. Poles inside the
    spherical triangle, clustered exponentially towards its corners (a lightning layout on the
    sphere), each with the fundamental solution and its two tangential derivatives.
    """
    def __init__(self, n=20, sigma=3.0, dmax=0.7, ni=200, seed=0):
        V, C = -np.eye(3), unit(-np.ones(3))
        self.poles, self.tans = [], []
        for v in V:
            t = unit(C - (C @ v) * v)
            for j in range(1, n + 1):
                dl = dmax * np.exp(-sigma * (np.sqrt(n) - np.sqrt(j)))
                p = np.cos(dl) * v + np.sin(dl) * t
                tp = unit(t - (t @ p) * p)
                self.poles.append(p)
                self.tans.append((tp, np.cross(p, tp)))

        # Boundary: the three sides, Gauss panels graded towards the corners
        g, w = np.polynomial.legendre.leggauss(4)
        edges = np.concatenate(([0.0], (np.pi / 4) * 0.6 ** np.arange(18, -1, -1)))
        th = np.concatenate([0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(edges[:-1], edges[1:])])
        wt = np.concatenate([0.5 * (b - a) * w for a, b in zip(edges[:-1], edges[1:])])
        th, wt = np.concatenate((th, np.pi / 2 - th[::-1])), np.concatenate((wt, wt[::-1]))
        self.xb = np.concatenate([np.outer(np.cos(th), V[i]) + np.outer(np.sin(th), V[(i + 1) % 3]) for i in range(3)])
        self.nb = np.concatenate([np.tile(V[(i + 2) % 3], (len(th), 1)) for i in range(3)])
        self.wb = np.tile(wt, 3)
        self.rho = np.arccos(np.clip(np.max(-self.xb, axis=1), -1, 1))   # distance to the nearest corner
        pts = unit(np.random.default_rng(seed).normal(size=(8 * ni, 3)))
        self.xi = pts[~np.all(pts < 0, axis=1)][:ni]

    def columns(self, nu, x, normal=None):
        out = []
        for s, tt in zip(self.poles, self.tans):
            c = np.clip(-(x @ s), -1 + 1e-15, 1 - 1e-15)
            if normal is None:
                out += [P(nu, c)] + [-dP(nu, c) * (x @ t) for t in tt]
            else:
                ns = np.sum(normal * s, axis=1)
                out += [-dP(nu, c) * ns] + [d2P(nu, c) * ns * (x @ t) - dP(nu, c) * (normal @ t) for t in tt]
        return np.column_stack(out)

    def exponents(self, neumann=False):
        """Exponents nu (Dirichlet) or tau (neumann=True) of the vertex, below 2.6"""
        def sigma(nu):
            # Neumann rows are weighted by the distance to the corner, where the normal derivative grows like 1/distance
            wb = np.sqrt(self.wb) * (self.rho if neumann else 1.0)
            B = self.columns(nu, self.xb, self.nb if neumann else None) * wb[:, None]
            return smallest_sigma(B, self.columns(nu, self.xi) / np.sqrt(len(self.xi)))
        return minima(sigma, 0.06, 2.6, 0.04)
