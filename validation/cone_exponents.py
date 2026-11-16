"""Singularity exponents of the electromagnetic field at the tip of a perfectly
conducting elliptic cone (the tangent cone of the NASA almond's tip).

Near the tip the field is a sum of two families (e.g. Van Bladel, Singular
Electromagnetic Fields and Sources, ch. 5; De Smedt and Van Bladel 1986):

    electric type:  E = grad( r^nu  Y(theta, phi) ),   E ~ r^(nu - 1),
    magnetic type:  H = grad( r^tau Z(theta, phi) ),   H ~ r^(tau - 1),  E ~ r^tau,

where Y (Z) is a Dirichlet (Neumann) eigenfunction of the Laplace-Beltrami
operator on the part of the unit sphere outside the cone, with eigenvalue
nu (nu + 1) (tau (tau + 1)).

The eigenvalues are computed with the method of particular solutions in its
stable form (Betcke and Trefethen 2005): for a trial nu the functions

    P_nu( -cos(angle between x and s_j) ),   s_j inside the cone,

satisfy the eigenvalue equation exactly on the sphere away from s_j (they are
the fundamental solutions of the sphere, with a logarithmic pole at s_j); nu is
an exponent when a combination vanishes on the boundary curve without
vanishing in the domain. The poles s_j lie on the focal segment of the
elliptic cross-section, clustered towards its ends, as the sources of the
almond itself do; classes odd in z use the derivative of the pole across the
segment. Each mirror-parity class (y -> -y, z -> -z) is computed separately.

    python3 cone_exponents.py            # almond tip, all classes, with checks
    python3 cone_exponents.py cube       # vertex of a cube (field outside a solid octant)

The cube vertex uses the same method with poles inside the spherical triangle
cut out by the cube, clustered exponentially towards its three corners along
the medians (a lightning-type layout on the sphere), each pole carrying the
fundamental solution and its two tangential derivatives. Check: the first
Dirichlet exponent, 0.45417, is the known value for the exterior of a cube
corner (the "Fichera corner").
"""
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import lpmv

# half-angles of the almond's tip cone, from the ogive formulas of Woo et al. (1993):
# slopes dy/dx = 4.83345 * T1 / (2.08335^2 * 0.96), dz/dx = 1.61115 * T1 / (2.08335^2 * 0.96), T1 = 0.58333
T1 = 0.583333
TAN_Y = 4.83345 * T1 / (2.08335 ** 2 * 0.96)
TAN_Z = 1.61115 * T1 / (2.08335 ** 2 * 0.96)


def leg(nu, x):
    return lpmv(0, nu, x)


def dleg(nu, x):
    """d/dx P_nu(x)."""
    return nu * (leg(nu - 1, x) - x * leg(nu, x)) / (1 - x * x)


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


class Cone:
    """Cone along -x with the tip at the origin and cross-section
    (y / tan_y)^2 + (z / tan_z)^2 = x^2; the field lives outside it."""

    def __init__(self, tan_y, tan_z, n_src=24, n_bnd=400, n_int=60, focal=0.98, seed=0):
        self.ty, self.tz = tan_y, tan_z
        psi = (np.arange(n_bnd) + 0.5) / n_bnd * (np.pi / 2)              # one quadrant; classes handle the rest
        self.xb = unit(np.column_stack((-np.ones(n_bnd), tan_y * np.cos(psi), tan_z * np.sin(psi))))
        tang = unit(np.column_stack((np.zeros(n_bnd), -tan_y * np.sin(psi), tan_z * np.cos(psi))))
        tang = unit(tang - np.sum(tang * self.xb, axis=1)[:, None] * self.xb)
        self.nb = np.cross(self.xb, tang)                                 # conormal: tangent to the sphere, normal to the curve
        # poles: on the focal segment (in the plane z = 0), Chebyshev clustering towards the foci
        f = focal * np.sqrt(max(tan_y ** 2 - tan_z ** 2, 0.0))
        u = np.cos((np.arange(n_src) + 0.5) / n_src * (np.pi / 2))        # y >= 0 half; mirror images come from the class
        self.src = unit(np.column_stack((-np.ones(n_src), f * u, np.zeros(n_src))))
        self.circular = f == 0.0
        if self.circular:
            self.src = self.src[:1]
        rng = np.random.default_rng(seed)
        pts = unit(rng.normal(size=(40 * n_int, 3)))
        pts = np.abs(pts) * np.array([1.0, 1.0, 1.0])                     # first quadrant in (y, z); x of either sign
        pts[:, 0] *= np.sign(rng.normal(size=len(pts)))
        inside = (pts[:, 0] < 0) & ((pts[:, 1] / tan_y) ** 2 + (pts[:, 2] / tan_z) ** 2 < pts[:, 0] ** 2)
        self.xi = pts[~inside][:n_int]

    def columns(self, nu, x, sy, sz, normal=None):
        """Symmetrised fundamental solutions at points x (value, or derivative along `normal`)."""
        cols = []
        for s in self.src:
            col = 0
            for ry, rz, sign in ((1, 1, 1.0), (-1, 1, sy), (1, -1, sz), (-1, -1, sy * sz)):
                if sz > 0:
                    sj = s * np.array([1.0, ry, 1.0])
                    if rz < 0:
                        continue                                           # the pole lies in the plane z = 0
                    c = np.clip(-(x @ sj), -1 + 1e-15, 1 - 1e-15)
                    col = col + sign * (leg(nu, c) if normal is None else -dleg(nu, c) * (normal @ sj))
                else:                                                      # odd in z: derivative of the pole across the segment
                    sj = s * np.array([1.0, ry, 1.0])
                    if rz < 0:
                        continue
                    ez = np.array([0.0, 0.0, 1.0])
                    c = np.clip(-(x @ sj), -1 + 1e-15, 1 - 1e-15)
                    if normal is None:
                        col = col + sign * (-dleg(nu, c) * (x @ ez))
                    else:
                        d2 = (nu * (nu + 1) * leg(nu, c) - 2 * c * dleg(nu, c)) / (c * c - 1)   # P'' from the Legendre equation
                        col = col + sign * (d2 * (normal @ sj) * (x @ ez) - dleg(nu, c) * (normal @ ez))
            if self.circular and (sy < 0 or sz < 0):
                raise ValueError("circular cone: only the axisymmetric class is implemented")
            cols.append(col)
        return np.column_stack(cols)

    def sigma(self, nu, sy, sz, neumann=False):
        """Smallest singular value of the boundary part of an orthonormal basis of
        the trial space sampled on boundary and interior: zero at an exponent."""
        B = self.columns(nu, self.xb, sy, sz, self.nb if neumann else None)
        I = self.columns(nu, self.xi, sy, sz)
        A = np.vstack((B / np.sqrt(len(B)), I / np.sqrt(len(I))))
        A = A / np.linalg.norm(A, axis=0)
        Q, R = np.linalg.qr(A)
        keep = np.abs(np.diag(R)) > 1e-13 * np.abs(np.diag(R)).max()
        if not keep.all():                                                # drop numerically dependent poles
            U, s, _ = np.linalg.svd(A, full_matrices=False)
            Q = U[:, s > 1e-13 * s[0]]
        return np.linalg.svd(Q[:len(B)], compute_uv=False)[-1]

    def exponents(self, sy, sz, neumann=False, lo=0.02, hi=3.6, step=0.02):
        """All exponents of the class in (lo, hi): local minima of sigma that reach (nearly) zero."""
        grid = np.arange(lo, hi + 1e-9, step)
        grid = grid[np.abs(grid - np.round(grid)) > 1e-6]                 # P_nu(-cos) loses its pole at integer nu
        sig = np.array([self.sigma(v, sy, sz, neumann) for v in grid])
        out = []
        for i in range(1, len(grid) - 1):
            if sig[i] < sig[i - 1] and sig[i] <= sig[i + 1]:
                r = minimize_scalar(lambda v: self.sigma(v, sy, sz, neumann), bounds=(grid[i - 1], grid[i + 1]),
                                    method="bounded", options={"xatol": 1e-10})
                if r.fun < 1e-3 * np.median(sig):
                    out.append((r.x, r.fun))
        return out


def d2leg(nu, x):
    """P_nu'' from Legendre's equation."""
    return (nu * (nu + 1) * leg(nu, x) - 2 * x * dleg(nu, x)) / (x * x - 1)


class CubeVertex:
    """Vertex at the origin of the solid octant x, y, z < 0; the field lives on the
    other seven octants. On the unit sphere the body is the triangle with corners
    -e_x, -e_y, -e_z and sides on the coordinate planes."""

    def __init__(self, n=20, sigma=3.0, dmax=0.7, ni=200, bq=0.6, blev=18, seed=0):
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
        g, w = np.polynomial.legendre.leggauss(4)                     # sides: panels graded towards both corners
        edges = np.concatenate(([0.0], (np.pi / 4) * bq ** np.arange(blev, -1, -1)))
        th = np.concatenate([0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(edges[:-1], edges[1:])])
        wt = np.concatenate([0.5 * (b - a) * w for a, b in zip(edges[:-1], edges[1:])])
        th, wt = np.concatenate((th, np.pi / 2 - th[::-1])), np.concatenate((wt, wt[::-1]))
        self.xb = np.concatenate([np.outer(np.cos(th), V[i]) + np.outer(np.sin(th), V[(i + 1) % 3]) for i in range(3)])
        self.nb = np.concatenate([np.tile(V[(i + 2) % 3], (len(th), 1)) for i in range(3)])   # conormal of each side
        self.wb = np.tile(wt, 3)
        self.rho = np.arccos(np.clip(np.max(-self.xb, axis=1), -1, 1))                       # distance to the nearest corner
        pts = unit(np.random.default_rng(seed).normal(size=(8 * ni, 3)))
        self.xi = pts[~np.all(pts < 0, axis=1)][:ni]

    def columns(self, nu, x, normal=None):
        out = []
        for s, tt in zip(self.poles, self.tans):
            c = np.clip(-(x @ s), -1 + 1e-15, 1 - 1e-15)
            if normal is None:
                out += [leg(nu, c)] + [-dleg(nu, c) * (x @ t) for t in tt]
            else:
                ns = np.sum(normal * s, axis=1)
                out += [-dleg(nu, c) * ns] + [d2leg(nu, c) * ns * (x @ t) - dleg(nu, c) * (normal @ t) for t in tt]
        return np.column_stack(out)

    def sigma(self, nu, neumann=False):
        # Neumann rows are weighted by the distance to the corner (the normal derivative grows like 1/distance there)
        B = self.columns(nu, self.xb, self.nb if neumann else None) * (np.sqrt(self.wb) * (self.rho if neumann else 1.0))[:, None]
        I = self.columns(nu, self.xi) / np.sqrt(len(self.xi))
        A = np.vstack((B, I))
        A = A / np.linalg.norm(A, axis=0)
        U, s, _ = np.linalg.svd(A, full_matrices=False)
        return np.linalg.svd(U[:len(B), s > 1e-12 * s[0]], compute_uv=False)[-1]

    def exponents(self, neumann=False, lo=0.06, hi=2.6, step=0.04):
        grid = np.arange(lo, hi, step)
        grid = grid[np.abs(grid - np.round(grid)) > 1e-6]
        sig = np.array([self.sigma(v, neumann) for v in grid])
        out = []
        for i in range(1, len(grid) - 1):
            if sig[i] < sig[i - 1] and sig[i] <= sig[i + 1]:
                r = minimize_scalar(lambda v: self.sigma(v, neumann), bounds=(grid[i - 1], grid[i + 1]),
                                    method="bounded", options={"xatol": 1e-9})
                out.append((r.x, r.fun))
        return out


def circular_exponent(alpha, neumann=False):
    """First axisymmetric exponent of a circular cone of half-angle alpha (exact condition)."""
    x = np.cos(np.pi - alpha)
    f = (lambda v: dleg(v, x)) if neumann else (lambda v: leg(v, x))
    return brentq(f, 0.05, 0.9999) if not neumann else brentq(f, 1.01, 2.5)


CLASSES = (("even y, even z", 1, 1), ("odd y, even z", -1, 1), ("even y, odd z", 1, -1), ("odd y, odd z", -1, -1))


def table(cone, hi=3.6):
    res = {}
    for name, sy, sz in CLASSES:
        for neumann in (False, True):
            ex = cone.exponents(sy, sz, neumann, hi=hi)
            res[(sy, sz, neumann)] = [e for e, _ in ex]
            print(f"  {name:15s} {'Neumann  (magnetic, tau)' if neumann else 'Dirichlet (electric, nu)'}: "
                  + "  ".join(f"{e:.5f}" for e, _ in ex) + "    (sigma_min: " + " ".join(f"{s:.0e}" for _, s in ex) + ")")
    return res


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "cube":
        print("cube vertex (exterior of a solid octant); known first Dirichlet exponent: 0.45417")
        for n in (8, 14, 20):
            v = CubeVertex(n=n)
            for neumann in (False, True):
                print(f"  poles per corner {n:2d}, {'Neumann  (magnetic, tau)' if neumann else 'Dirichlet (electric, nu)'}: "
                      + "  ".join(f"{e:.5f} ({s:.0e})" for e, s in v.exponents(neumann)))
        sys.exit()
    print("check 1: circular cone, half-angle 20 degrees (exact: zero of P_nu(cos 160 deg))")
    c = Cone(np.tan(np.deg2rad(20.0)), np.tan(np.deg2rad(20.0)))
    print(f"  particular solutions: {c.exponents(1, 1)[0][0]:.8f}   exact: {circular_exponent(np.deg2rad(20.0)):.8f}")
    print("check 2: nearly circular elliptic cone (20 and 19.9 degrees) against the circular one")
    c = Cone(np.tan(np.deg2rad(20.0)), np.tan(np.deg2rad(19.9)))
    print(f"  particular solutions: {c.exponents(1, 1, hi=1.0)[0][0]:.6f}   circular 20: {circular_exponent(np.deg2rad(20.0)):.6f}"
          f"   circular 19.9: {circular_exponent(np.deg2rad(19.9)):.6f}")
    print(f"almond tip: half-angles {np.rad2deg(np.arctan(TAN_Y)):.2f} and {np.rad2deg(np.arctan(TAN_Z)):.2f} degrees")
    print(f"  circular cones with these half-angles: nu = {circular_exponent(np.arctan(TAN_Y)):.4f} and "
          f"{circular_exponent(np.arctan(TAN_Z)):.4f}")
    for n_src in (12, 24, 36):
        print(f" poles per half focal segment: {n_src}")
        table(Cone(TAN_Y, TAN_Z, n_src=n_src), hi=3.6 if n_src == 24 else 1.6)
