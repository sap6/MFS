"""MFS for exterior scattering by a perfectly conducting unit sphere,
validated against the Mie series.

Sources: N electric dipoles (3 unknowns each) on a Fibonacci lattice of radius
R < 1. Boundary condition: tangential E_sc = -tangential E_inc at M
collocation points (2M equations), solved in the least-squares sense.

Errors reported (all relative, discrete L2):
  bc    tangential boundary residual at independent test points on r = 1
  near  E_sc vs Mie on r = 1.05
  mid   E_sc vs Mie on r = 1.5
  far   far-field pattern vs Mie (evaluated at r = 1e4)

Usage: python pec_sphere_convergence.py [--quick]
"""
import json
import sys
import time

import numpy as np
from scipy.linalg import lstsq

import mie

D = np.array([0.0, 0.0, -1.0])  # propagation direction
P = np.array([0.0, 1.0, 0.0])   # polarisation


def fibonacci_sphere(n):
    j = np.arange(n)
    z = 1 - (2 * j + 1) / n
    th = 2 * np.pi * j / ((1 + 5 ** 0.5) / 2)
    s = np.sqrt(1 - z ** 2)
    return np.column_stack((s * np.cos(th), s * np.sin(th), z))


def dyadic(x, y, k):
    """(I + grad grad / k^2) e^{ikr}/r, shape (len(x), len(y), 3, 3)."""
    dvec = x[:, None] - y[None]
    r = np.linalg.norm(dvec, axis=2)
    kr = k * r
    g = np.exp(1j * kr) / r
    fI = (1 + 1j / kr - 1 / kr ** 2) * g
    fR = (1 + 3j / kr - 3 / kr ** 2) * g
    rh = dvec / r[..., None]
    return (fI[..., None, None] * np.eye(3)
            - fR[..., None, None] * rh[..., :, None] * rh[..., None, :])


def tangent_basis(n):
    ref = np.where(np.abs(n[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    t1 = np.cross(n, ref)
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    return np.stack((t1, np.cross(n, t1)), axis=1)  # (M, 2, 3)


def evaluate(pts, y, alpha, k, batch=400):
    out = np.empty((len(pts), 3), complex)
    for i in range(0, len(pts), batch):
        out[i:i + batch] = np.einsum("pnab,nb->pa", dyadic(pts[i:i + batch], y, k), alpha)
    return out


def rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


class Reference:
    """Mie reference data on fixed test surfaces, computed once per k."""

    def __init__(self, k, ntest=1500):
        self.k = k
        self.u = fibonacci_sphere(ntest) @ _random_rotation()
        self.radii = {"near": 1.05, "mid": 1.5, "far": 1e4}
        self.exact = {name: mie.scattered_field(r * self.u, k, D, P)
                      for name, r in self.radii.items()}
        self.T = tangent_basis(self.u)
        self.rhs = -np.einsum("mab,mb->ma", self.T, mie.plane_wave(self.u, k, D, P))


def _random_rotation():
    q, _ = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))
    return q


def solve(k, N, R, ref, M=None):
    M = M or 2 * N
    x = fibonacci_sphere(M)
    y = R * fibonacci_sphere(N)
    T = tangent_basis(x)
    A = np.einsum("mab,mnbc->manc", T, dyadic(x, y, k)).reshape(2 * M, 3 * N)
    b = -np.einsum("mab,mb->ma", T, mie.plane_wave(x, k, D, P)).reshape(-1)
    t0 = time.time()
    alpha, _, rank, sv = lstsq(A, b)
    t_solve = time.time() - t0
    alpha = alpha.reshape(N, 3)

    res = {"k": k, "N": N, "M": M, "R": R, "rank": int(rank),
           "cond": float(sv[0] / sv[-1]), "alpha_norm": float(np.linalg.norm(alpha)),
           "t_solve": t_solve}
    bc = np.einsum("mab,mb->ma", ref.T, evaluate(ref.u, y, alpha, k))
    res["bc"] = rel(bc, ref.rhs)
    for name, r in ref.radii.items():
        res[name] = rel(evaluate(r * ref.u, y, alpha, k), ref.exact[name])
    return res


def main(quick=False):
    results = []
    ks = [2 * np.pi, 4 * np.pi]
    Ns = [50, 100, 200, 300] if quick else [50, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200]
    Rs_fixed = [0.5, 0.8] if quick else [0.3, 0.5, 0.7, 0.9]
    Rs_sweep = np.round(np.linspace(0.1, 0.95, 6 if quick else 18), 3)
    N_sweep = 200 if quick else 400

    for k in ks:
        ref = Reference(k)
        for R in Rs_fixed:
            for N in Ns:
                r = solve(k, N, R, ref)
                r["sweep"] = "N"
                results.append(r)
                print(f"k={k:6.3f} R={R:.2f} N={N:5d}  bc={r['bc']:.1e} near={r['near']:.1e} "
                      f"mid={r['mid']:.1e} far={r['far']:.1e} |a|={r['alpha_norm']:.1e} "
                      f"cond={r['cond']:.1e} t={r['t_solve']:.1f}s", flush=True)
        for R in Rs_sweep:
            r = solve(k, N_sweep, float(R), ref)
            r["sweep"] = "R"
            results.append(r)
            print(f"k={k:6.3f} R={R:.3f} N={N_sweep:5d}  bc={r['bc']:.1e} near={r['near']:.1e} "
                  f"mid={r['mid']:.1e} far={r['far']:.1e} |a|={r['alpha_norm']:.1e} "
                  f"cond={r['cond']:.1e}", flush=True)

    out = "results_pec_sphere_quick.json" if quick else "results_pec_sphere.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
