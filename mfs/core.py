"""
Building blocks shared by every shape: the incident wave, the fields of
electric and magnetic dipoles, points on a sphere, tangent bases, the
least-squares solver and the far field.
"""
import numpy as np
from scipy.linalg import blas, lapack, lstsq, solve_triangular


def plane_wave(x, k, d, p):
    """
    E(x) = p e^{ik d·x}, d orthogonal to p
      x: points, shape (M, 3)
      Returns an array of shape (M, 3)
    """
    return np.asarray(p)[None, :] * np.exp(1j * k * (x @ np.asarray(d)))[:, None]


def wave(theta=0.0, phi=0.0, psi=90.0):
    """
    Plane wave arriving from the direction (theta, phi), in degrees. psi turns
    the polarisation from the theta-unit vector towards the phi-unit vector.
    The default is the standard wave of the paper: along -z, polarised along y.
      Returns the propagation direction d and the polarisation p
    """
    t, f, s = np.radians([theta, phi, psi])
    rhat = np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
    that = np.array([np.cos(t) * np.cos(f), np.cos(t) * np.sin(f), -np.sin(t)])
    fhat = np.array([-np.sin(f), np.cos(f), 0.0])
    return -rhat, np.cos(s) * that + np.sin(s) * fhat


def dyadic(x, y, k):
    """
    G(x, y) = (I + ∇∇/k²) e^{ikr}/r, the field of an electric dipole at y
      x: field points, shape (M, 3)
      y: source points, shape (N, 3)
      Returns an array of shape (M, N, 3, 3)
    """
    d = x[:, None] - y[None]                       # (M, N, 3)
    r = np.linalg.norm(d, axis=2)                  # (M, N)
    kr = k * r
    g = np.exp(1j * kr) / r
    fI = (1 + 1j / kr - 1 / kr ** 2) * g
    fR = (1 + 3j / kr - 3 / kr ** 2) * g
    rhat = d / r[..., None]
    return fI[..., None, None] * np.eye(3) - fR[..., None, None] * rhat[..., :, None] * rhat[..., None, :]


def magnetic(x, y, k):
    """
    (1/k) ∇ × (e^{ikr}/r β), the field of a magnetic dipole at y, as a matrix acting on β
      x: field points, shape (M, 3)
      y: source points, shape (N, 3)
      Returns an array of shape (M, N, 3, 3)
    """
    d = x[:, None] - y[None]
    r = np.linalg.norm(d, axis=2)
    f = (1j - 1 / (k * r)) * np.exp(1j * k * r) / r
    e = d / r[..., None]
    out = np.zeros(d.shape + (3,), complex)
    out[..., 0, 1], out[..., 0, 2] = -e[..., 2], e[..., 1]
    out[..., 1, 0], out[..., 1, 2] = e[..., 2], -e[..., 0]
    out[..., 2, 0], out[..., 2, 1] = -e[..., 1], e[..., 0]
    return f[..., None, None] * out


def fibonacci(n):
    """
    Evenly distribute n points on the unit sphere (Fibonacci lattice)
    https://arxiv.org/pdf/0912.4540.pdf
    """
    i = np.arange(n)
    z = 1 - (2 * i + 1) / n
    theta = 2 * np.pi * i / ((1 + 5 ** 0.5) / 2)
    s = np.sqrt(1 - z ** 2)
    return np.column_stack((s * np.cos(theta), s * np.sin(theta), z))


def tangent_basis(n):
    """
    n: unit normals, shape (M, 3)
    Returns two unit tangents at each point, shape (M, 2, 3)
    """
    ref = np.where(np.abs(n[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    t1 = np.cross(n, ref)
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    return np.stack((t1, np.cross(n, t1)), axis=1)


def field(pts, y, alpha, k, batch=400):
    """
    E_N(x) = Σ_j G(x, y_j) α_j, evaluated in batches
      pts  : (P, 3) evaluation points
      y    : (N, 3) source points
      alpha: (N, 3) coefficients
    """
    out = np.empty((len(pts), 3), complex)
    for i in range(0, len(pts), batch):
        out[i:i + batch] = np.einsum("pnab,nb->pa", dyadic(pts[i:i + batch], y, k), alpha)
    return out


def rel(a, b):
    """Relative error |a - b| / |b| in the discrete 2-norm"""
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


# Least squares

class QR:
    """
    Householder QR of A with implicit column scaling, kept in the factored form
    of LAPACK so that each right-hand side costs one application of Q^H
    """
    def __init__(self, A):
        A = np.array(A, dtype=complex, order="F")    # a copy that LAPACK may overwrite
        self.scale = np.linalg.norm(A, axis=0)
        A /= self.scale
        self.n = A.shape[1]
        lwork = int(lapack.zgeqrf_lwork(*A.shape)[0].real)
        self.qr, self.tau, _, _ = lapack.zgeqrf(A, lwork=lwork, overwrite_a=True)
        self.R = np.triu(self.qr[:self.n])
        self.dR = np.abs(np.diag(self.R))

    def solve(self, b, residual=False, cond=1e-12):
        """Minimiser of |A c - b| (and the norm of the residual, if asked for)"""
        B = np.array(b.reshape(len(b), -1), dtype=complex, order="F")
        lwork = int(lapack.zunmqr("L", "C", self.qr, self.tau, B, -1)[1][0].real)
        y = lapack.zunmqr("L", "C", self.qr, self.tau, B, lwork, overwrite_c=True)[0].reshape(b.shape)
        if self.dR.min() > cond * self.dR.max():
            c = solve_triangular(self.R, y[:self.n])
        else:                                          # numerically rank deficient: rank-revealing solve
            c = lstsq(self.R, y[:self.n], lapack_driver="gelsy", cond=cond)[0]
        c = c / (self.scale if c.ndim == 1 else self.scale[:, None])
        return (c, np.linalg.norm(y[self.n:], axis=0)) if residual else c


class CSNE:
    """
    Corrected semi-normal equations (Björck 1987): Cholesky factor R of the
    Gram matrix A^H A, the solution of R^H R c = A^H b, and two steps of
    refinement with the residual b - A c. About half the work of QR, and as
    accurate while cond(A)² is well below 1/eps. Raises LinAlgError otherwise.
    """
    def __init__(self, A, max_cond=1e7):
        self.A, self.n = A, A.shape[1]
        G = np.conj(blas.zherk(1.0, A.T, trans=0, lower=0))    # upper triangle of A^H A
        self.scale = np.sqrt(np.real(np.diag(G)))              # column norms
        G /= self.scale[:, None] * self.scale[None, :]
        self.R, info = lapack.zpotrf(G, lower=0, overwrite_a=True)
        self.dR = np.abs(np.diag(self.R))
        if info != 0 or self.dR.max() > max_cond * 1e-3 * self.dR.min():
            raise np.linalg.LinAlgError("too ill-conditioned for the semi-normal equations")

    def _normal(self, r):
        # (A^H A)^{-1} A^H r, with the column scaling undone
        rhs = np.conj(self.A.T @ np.conj(r))
        s = self.scale if rhs.ndim == 1 else self.scale[:, None]
        return lapack.zpotrs(self.R, rhs / s, lower=0)[0] / s

    def solve(self, b, residual=False, cond=None):
        c = self._normal(b)
        for _ in range(2):
            c = c + self._normal(b - self.A @ c)
        return (c, np.linalg.norm(b - self.A @ c, axis=0)) if residual else c


def factor(A):
    """Semi-normal equations if A is well conditioned, QR otherwise"""
    try:
        return CSNE(A)
    except np.linalg.LinAlgError:
        return QR(A)


# Far field and cross-sections

def sphere_rule(nt=40, nphi=80):
    """
    Gauss-Legendre rule in cos(theta), trapezoidal rule in phi
      Returns directions (nt * nphi, 3) and weights (nt * nphi,)
    """
    g, w = np.polynomial.legendre.leggauss(nt)
    phi = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    ct, ph = np.meshgrid(g, phi, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    dirs = np.column_stack(((st * np.cos(ph)).ravel(), (st * np.sin(ph)).ravel(), ct.ravel()))
    return dirs, np.repeat(w, nphi) * 2 * np.pi / nphi


def cross_sections(far, k, d, p, rule=None):
    """
    Extinction and scattering cross-sections from the far-field pattern F, where
    E ~ F(x̂) e^{ikr}/r. For a perfect conductor the two are equal (optical
    theorem), so their mismatch measures the error without a reference.
      far: function of directions (D, 3) returning F, shape (D, 3)
    """
    dirs, w = rule or sphere_rule()
    F = far(np.vstack((dirs, d[None])))
    ext = 4 * np.pi / k * float(np.imag(np.conj(p) @ F[-1]))
    sca = float(np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1)))
    return ext, sca
