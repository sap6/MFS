"""
Scattering by the perfectly conducting unit sphere (Section 3.1 of the paper).
Electric dipoles on a sphere of radius R < 1, boundary points on the unit
sphere, both on Fibonacci lattices, and the Mie series as the exact solution.
"""
import numpy as np

from .core import QR, electric_block, fibonacci, plane_wave, tangent_basis


def solve(N, k=2 * np.pi, R=0.5, d=(0, 0, -1.0), p=(0, 1.0, 0), M=None):
    """
    Least-squares MFS for the sphere
      N: number of source points (3N unknowns)
      R: radius of the source sphere
      M: number of boundary points (default 2N, so 4N equations)
      Returns the source points y (N, 3) and the coefficients alpha (N, 3)
    """
    M = M or 2 * N

    # Boundary and source points
    x = fibonacci(M)                             # (M, 3), also the unit normals
    y = R * fibonacci(N)                         # (N, 3)

    # Form F-ordered matrix A directly: tangential components of the dipole fields, (2M, 3N)
    T = np.ascontiguousarray(tangent_basis(x))   # (M, 2, 3)
    A = np.empty((2 * M, 3 * N), complex, order="F")
    electric_block(x, T, y, float(k), A)

    # Form vector b: minus the tangential incident field
    b = -np.einsum("mab,mb->ma", T, plane_wave(x, k, d, p)).reshape(-1)

    # Solve the least-squares problem for alpha
    alpha = QR(A, overwrite=True).solve(b).reshape(N, 3)
    return y, alpha
