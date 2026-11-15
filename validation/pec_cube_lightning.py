"""Exploratory: MFS for exterior scattering by a perfectly conducting unit cube,
comparing uniform source placement with "lightning"-style placement in which
dipole sources cluster exponentially towards the edges and vertices
(after Gopal & Trefethen's lightning solvers, carried over to 3D edges).

There is no closed-form solution for the cube, so accuracy is measured by
  bc_l2 / bc_max  relative tangential boundary residual at independent test
                  points on the faces (area-weighted L2, and max)
  far             far-field change relative to the most refined run

Usage: python pec_cube_lightning.py [--quick]
"""
import itertools
import json
import sys
import time

import numpy as np
from scipy.linalg import lstsq

import mie
from pec_sphere_convergence import dyadic, evaluate, fibonacci_sphere

K = 2 * np.pi                    # cube side = one wavelength
D = np.array([0.0, 0.0, -1.0])
P = np.array([0.0, 1.0, 0.0])
H = 0.5                          # half side length


# ---------- boundary sampling -------------------------------------------------

def graded_nodes(levels, p, q=0.5):
    """Gauss-Legendre nodes/weights on (-H, H), panels graded towards both ends.
    levels = 0 gives uniform panels."""
    if levels == 0:
        brk = np.linspace(0, H, 5)
    else:
        brk = np.concatenate(([0.0], H * q ** np.arange(levels, -1, -1)))  # distance from edge
    half_brk = H - brk[::-1]                     # 0 .. H measured from the centre
    brk_full = np.concatenate((-half_brk[::-1], half_brk[1:]))
    g, w = np.polynomial.legendre.leggauss(p)
    xs, ws = [], []
    for a, b in zip(brk_full[:-1], brk_full[1:]):
        xs.append(0.5 * (a + b) + 0.5 * (b - a) * g)
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(xs), np.concatenate(ws)


def cube_surface(levels, p):
    """Points, quadrature weights and two unit tangents on the cube faces."""
    s, w = graded_nodes(levels, p)
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


# ---------- source placement --------------------------------------------------

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


def lightning_sources(levels, depth0, q=0.5, spacing=1.0):
    """Sources clustered towards the 12 edges and 8 vertices.

    Layer j sits at distance d_j = depth0 * q**j from the edge (along the inward
    bisector), with spacing along the edge proportional to d_j, plus one source
    per layer on each inward body diagonal."""
    out = []
    for j in range(levels + 1):
        d = depth0 * q ** j
        off = H - d / np.sqrt(2)
        nt = max(1, int(round((2 * H - 2 * d) / (spacing * d))))
        t = (np.arange(nt) + 0.5) / nt * (2 * H - 2 * d) - (H - d)
        for ax in range(3):
            b, c = [i for i in range(3) if i != ax]
            for s1, s2 in itertools.product((-1, 1), repeat=2):
                y = np.zeros((nt, 3))
                y[:, ax], y[:, b], y[:, c] = t, s1 * off, s2 * off
                out.append(y)
        for sg in itertools.product((-1, 1), repeat=3):
            out.append((np.array(sg) * (H - d / np.sqrt(3)))[None, :])
    return np.concatenate(out)


# ---------- solve -------------------------------------------------------------

def incident(pts):
    return mie.plane_wave(pts, K, D, P)


def solve(y, col_levels, col_p, test, inc=incident):
    x, w, T = cube_surface(col_levels, col_p)
    sw = np.sqrt(w)
    A = np.einsum("mab,mnbc->manc", T, dyadic(x, y, K))
    A = (A * sw[:, None, None, None]).reshape(2 * len(x), 3 * len(y))
    b = (-np.einsum("mab,mb->ma", T, inc(x)) * sw[:, None]).reshape(-1)
    scale = np.linalg.norm(A, axis=0)            # column normalisation
    t0 = time.time()
    c, _, rank, _ = lstsq(A / scale, b, lapack_driver="gelsy", cond=1e-13)
    t_solve = time.time() - t0
    alpha = (c / scale).reshape(len(y), 3)

    xt, wt, Tt = test["x"], test["w"], test["T"]
    res_t = np.einsum("mab,mb->ma", Tt, evaluate(xt, y, alpha, K)) - test["rhs"]
    num = np.sum(np.abs(res_t) ** 2, axis=1)
    out = {"N": len(y), "M": len(x), "rank": int(rank), "t_solve": t_solve,
           "bc_l2": float(np.sqrt(np.sum(wt * num) / test["rhs_norm2"])),
           "bc_max": float(np.sqrt(num.max())),
           "coef_norm": float(np.linalg.norm(c))}
    return out, evaluate(test["far_pts"], y, alpha, K)


def make_test(inc=incident):
    x, w, T = cube_surface(8, 3)                 # independent, finer than any collocation set
    rhs = -np.einsum("mab,mb->ma", T, inc(x))
    return {"x": x, "w": w, "T": T, "rhs": rhs,
            "rhs_norm2": float(np.sum(w * np.sum(np.abs(rhs) ** 2, axis=1))),
            "far_pts": 1e4 * fibonacci_sphere(400)}


def main(quick=False):
    test = make_test()
    depth0 = 0.2
    runs = []
    # uniform placement: inner-cube grids of increasing size
    for n in ([4, 6, 8] if quick else [4, 6, 8, 11]):
        runs.append(("uniform", f"n={n}", base_sources(n, depth0), 0, max(3, int(np.ceil(n / 2.4)))))
    # lightning placement: fixed coarse base grid + edge/vertex clusters
    for L in ([0, 1, 2] if quick else [0, 1, 2, 3]):
        y = np.concatenate((base_sources(6, depth0), lightning_sources(L, depth0)))
        runs.append(("lightning", f"levels={L}", y, L + 2, 3))

    results, fars = [], []
    for kind, label, y, col_levels, col_p in runs:
        r, far = solve(y, col_levels, col_p, test)
        r.update(kind=kind, label=label)
        results.append(r)
        fars.append(far)
        print(f"{kind:9s} {label:9s} N={r['N']:5d} M={r['M']:5d} rank={r['rank']:5d} "
              f"bc_l2={r['bc_l2']:.2e} bc_max={r['bc_max']:.2e} t={r['t_solve']:.1f}s", flush=True)

    ref = fars[-1]                               # most refined lightning run
    for r, far in zip(results, fars):
        r["far"] = float(np.linalg.norm(far - ref) / np.linalg.norm(ref))
    print("\nfar-field difference from the most refined lightning run:")
    for r in results:
        print(f"{r['kind']:9s} {r['label']:9s} N={r['N']:5d}  far={r['far']:.2e}")

    out = "results_pec_cube_quick.json" if quick else "results_pec_cube.json"
    json.dump(results, open(out, "w"), indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
