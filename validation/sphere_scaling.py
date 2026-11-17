"""Sphere: how the number of sources grows with the wavenumber.

For k = 2 pi, 4 pi, 6 pi, 8 pi (diameter 2 to 8 wavelengths) and sources on
the sphere R = 0.5, prints the error against Mie for increasing N and the
smallest N that reaches 1e-3 and 1e-6.

    python3 sphere_scaling.py
"""
import time

import numpy as np

import mie
from lsq import qr_solve
from pec_sphere_convergence import D, P, dyadic, evaluate, fibonacci_sphere, rel, tangent_basis

NS = (100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400, 3200)


def mfs(N, k, pts, R=0.5):
    t0 = time.time()
    x, y = fibonacci_sphere(2 * N), R * fibonacci_sphere(N)
    T = tangent_basis(x)
    A = np.empty((4 * N, 3 * N), complex, order="F")
    for i in range(0, 2 * N, 400):                      # in slabs, to bound the memory of the kernel
        sl = slice(i, i + 400)
        A[2 * i:2 * i + 2 * len(x[sl])] = np.einsum("mab,mnbc->manc", T[sl], dyadic(x[sl], y, k)).reshape(-1, 3 * N)
    b = -np.einsum("mab,mb->ma", T, mie.plane_wave(x, k, D, P)).reshape(-1)
    alpha = qr_solve(A, b, cond=1e-14)[0].reshape(N, 3)
    return evaluate(pts, y, alpha, k), time.time() - t0


if __name__ == "__main__":
    pts = 2.0 * fibonacci_sphere(400)
    mfs(50, 2 * np.pi, pts)
    for mult in (2, 4, 6, 8):
        k = mult * np.pi
        exact = mie.scattered_field(pts, k, D, P)
        first = {}
        for N in NS:
            E, t = mfs(N, k, pts)
            e = rel(E, exact)
            print(f"k = {mult} pi  N = {N:5d}  unknowns = {3 * N:5d}  time = {t:6.1f}s  error = {e:.2e}", flush=True)
            for tol in (1e-3, 1e-6):
                if e < tol and tol not in first:
                    first[tol] = (N, t)
            if e < 1e-9:
                break
        print(f"k = {mult} pi (ka = {k:.1f}): first N below 1e-3: {first.get(1e-3)}, below 1e-6: {first.get(1e-6)}", flush=True)
