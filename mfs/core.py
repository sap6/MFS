"""
Building blocks shared by every shape: the incident wave, the fields of
electric and magnetic dipoles, points on a sphere, tangent bases, the
least-squares solver and the far field.
"""
import numpy as np
import scipy.linalg as sla
from numba import njit, prange
from scipy.linalg import blas, lapack


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
    G(x, y) = (I + ∇∇/k²) e^{ikr}/r, the field of an electric dipole at y (with 1/(4π) absorbed into the moment)
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
    (1/k) ∇ × (e^{ikr}/r β), the field of a magnetic dipole at y (1/(4π) absorbed), as a matrix acting on β
      x: field points, shape (M, 3)
      y: source points, shape (N, 3)
      Returns an array of shape (M, N, 3, 3)
    """
    d = x[:, None] - y[None]                       # (M, N, 3)
    r = np.linalg.norm(d, axis=2)                  # (M, N)
    f = (1j - 1 / (k * r)) * np.exp(1j * k * r) / r
    e = d / r[..., None]
    return f[..., None, None] * np.cross(np.eye(3), e[..., None, :])   # the matrix of β → e × β


@njit(fastmath=True, inline="always")
def _electric_factors(dx0, dx1, dx2, k, inv_k):
    r2 = dx0 * dx0 + dx1 * dx1 + dx2 * dx2
    inv_r = 1.0 / np.sqrt(r2)
    kr = k * (r2 * inv_r)
    ik = inv_r * inv_k
    ik2 = ik * ik
    ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
    return (1.0 + 1j * ik - ik2) * ph, (1.0 + 3j * ik - 3.0 * ik2) * ph * (inv_r * inv_r)


@njit(parallel=True, fastmath=True, cache=True)
def electric_block(x, T, y, k, out):
    """out[m * nc + a, 3 * n + c] = T[m, a] · G(x[m], y[n]) e_c, assembled columnwise into F-ordered out"""
    M, nc, N = x.shape[0], T.shape[1], y.shape[0]
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    for n in prange(N):
        yn0, yn1, yn2 = y[n, 0], y[n, 1], y[n, 2]
        c0 = 3 * n
        for m in range(M):
            dx0, dx1, dx2 = x[m, 0] - yn0, x[m, 1] - yn1, x[m, 2] - yn2
            fI, fR = _electric_factors(dx0, dx1, dx2, k, inv_k)
            mt, row = m * tmul, m * nc
            for a in range(nc):
                t0, t1, t2 = T[mt, a, 0], T[mt, a, 1], T[mt, a, 2]
                td = fR * (t0 * dx0 + t1 * dx1 + t2 * dx2)
                out[row + a, c0] = fI * t0 - td * dx0
                out[row + a, c0 + 1] = fI * t1 - td * dx1
                out[row + a, c0 + 2] = fI * t2 - td * dx2


@njit(parallel=True, fastmath=True, cache=True)
def _field_kernel(pts, y, alpha, k, out):
    P, N = pts.shape[0], y.shape[0]
    inv_k = 1.0 / k
    for i in prange(P):
        px0, px1, px2 = pts[i, 0], pts[i, 1], pts[i, 2]
        s0, s1, s2 = 0.0j, 0.0j, 0.0j
        for j in range(N):
            dx0, dx1, dx2 = px0 - y[j, 0], px1 - y[j, 1], px2 - y[j, 2]
            fI, fR = _electric_factors(dx0, dx1, dx2, k, inv_k)
            a0, a1, a2 = alpha[j, 0], alpha[j, 1], alpha[j, 2]
            da = fR * (dx0 * a0 + dx1 * a1 + dx2 * a2)
            s0 += fI * a0 - da * dx0
            s1 += fI * a1 - da * dx1
            s2 += fI * a2 - da * dx2
        out[i, 0], out[i, 1], out[i, 2] = s0, s1, s2


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
    Two unit tangents t1 = n × e / |n × e|, t2 = n × t1 at each point (e = z, or x near the poles)
      n: unit normals, shape (M, 3)
      Returns an array of shape (M, 2, 3)
    """
    ref = np.where(np.abs(n[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    t1 = np.cross(n, ref)
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    return np.stack((t1, np.cross(n, t1)), axis=1)


def field(pts, y, alpha, k):
    """
    E_N(x) = Σ_j G(x, y_j) α_j, evaluated without forming the 4-tensor G
      pts  : evaluation points, shape (P, 3)
      y    : source points, shape (N, 3)
      alpha: coefficients, shape (N, 3)
      Returns an array of shape (P, 3)
    """
    pts = np.ascontiguousarray(pts, dtype=float)
    out = np.empty((len(pts), 3), complex)
    if len(pts) == 0 or len(y) == 0:
        out[:] = 0.0
        return out
    _field_kernel(pts, np.ascontiguousarray(y, dtype=float), np.ascontiguousarray(alpha, dtype=complex), float(k), out)
    return out


def rel(a, b):
    """Relative error |a - b| / |b| in the discrete 2-norm"""
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


# Symmetry: a wave splits into parts of definite parity under a group of reflections

def sign(cls, R):
    """Sign of the parity class cls (one sign per axis) under the reflection R (a vector of signs)"""
    return float(np.prod([s for s, r in zip(cls, R) if r < 0]))


def class_part(x, k, d, p, refl, cls):
    """
    (1/|G|) Σ_R sign(cls, R) R E(R x): the part of the wave E = p e^{ik d·x} with the parities cls
      refl: the reflections of the group G, shape (|G|, 3)
      d, p: shape (3,) for one wave, or (W, 3) for W waves
      Returns an array of shape (M, 3), or (M, W, 3)
    """
    D, P = np.atleast_2d(d), np.atleast_2d(p)                  # (W, 3)
    E = 0
    for R in refl:
        phase = np.exp(1j * k * (x * R) @ D.T)                   # e^{ik d·Rx}, shape (M, W)
        E = E + sign(cls, R) * R * P[None] * phase[:, :, None]
    return (E[:, 0] if np.ndim(d) == 1 else E) / len(refl)


def classes_of(k, d, p, refl):
    """The parity classes that the wave (d, p) excites; a class is labelled by its signs, like the reflections"""
    x = np.random.default_rng(0).uniform(-0.5, 0.5, (40, 3))
    tol = 1e-12 * max(float(np.linalg.norm(p)), 1e-300)
    return [tuple(int(v) for v in R) for R in refl if np.abs(class_part(x, k, d, p, refl, tuple(R))).max() > tol]


# Least squares

class QR:
    """Householder QR of A with its columns scaled to unit norm; Q is kept as its reflectors, not formed"""
    def __init__(self, A, overwrite=False):
        self.scale = np.linalg.norm(A, axis=0)
        A = A if (overwrite and A.flags.f_contiguous) else np.array(A, dtype=complex, order="F")
        A /= self.scale
        lwork = int(lapack.zgeqrf(A, lwork=-1)[2][0].real)
        self.qr, self.tau, _, _ = lapack.zgeqrf(A, lwork=lwork, overwrite_a=1)

    def solve(self, b, residual=False):
        """Minimiser of |A c - b| (and the norm of the residual, if asked for)"""
        B = np.array(b.reshape(len(b), -1), dtype=complex, order="F")
        lwork = int(lapack.zunmqr("L", "C", self.qr, self.tau, B, -1)[1][0].real)
        y = lapack.zunmqr("L", "C", self.qr, self.tau, B, lwork, overwrite_c=True)[0].reshape(b.shape)   # Q^H b
        n, s = self.qr.shape[1], (self.scale if y.ndim == 1 else self.scale[:, None])
        c = sla.solve_triangular(self.qr[:n, :n], y[:n], check_finite=False) / s                         # R c = Q^H b
        return (c, np.linalg.norm(y[n:], axis=0)) if residual else c


class CSNE:
    """
    Corrected semi-normal equations (Björck 1987): Cholesky factor R of the
    Gram matrix A^H A, the solution of R^H R c = A^H b, and two steps of
    refinement with the residual b - A c. About half the work of QR, and as
    accurate while cond(A)² is well below 1/eps. Raises LinAlgError otherwise.
    """
    def __init__(self, A, max_cond=3e7):
        self.A = A
        if A.flags.f_contiguous:
            G = blas.zherk(1.0, A, trans=2, lower=0)           # upper triangle of A^H A, without copying A
        else:                                                  # A is C-ordered, so A.T is F-ordered: (A^H A)^T = conj(A^H A)
            G = blas.zherk(1.0, A.T, trans=0, lower=0)
            np.conj(G, out=G)
        self.scale = np.sqrt(np.real(np.diag(G)))              # column norms
        G /= self.scale[:, None]
        G /= self.scale[None, :]
        self.R = sla.cho_factor(G, overwrite_a=True, check_finite=False)   # raises LinAlgError if not positive definite
        dR = np.abs(np.diag(self.R[0]))
        if dR.max() > max_cond * 1e-3 * dR.min():             # the diagonal of R underestimates cond(A) by up to 1e3
            raise np.linalg.LinAlgError("too ill-conditioned for the semi-normal equations")

    def _fwd(self, c):
        if self.A.flags.f_contiguous:
            return blas.zgemv(1.0, self.A, c) if c.ndim == 1 else blas.zgemm(1.0, self.A, c)
        return self.A @ c

    def _normal(self, r):
        if self.A.flags.f_contiguous:
            rhs = blas.zgemv(1.0, self.A, r, trans=2) if r.ndim == 1 else blas.zgemm(1.0, self.A, r, trans_a=2)
        else:
            rhs = np.conj(self.A.T @ np.conj(r))               # A^H r, without copying A
        s = self.scale if rhs.ndim == 1 else self.scale[:, None]
        return sla.cho_solve(self.R, rhs / s, check_finite=False) / s

    def solve(self, b, residual=False):
        c = self._normal(b)
        for _ in range(2):
            r = b - self._fwd(c)
            c = c + self._normal(r)
        if not residual:
            return c
        r = b - self._fwd(c)
        return c, np.linalg.norm(r, axis=0)


def factor(A, overwrite=False):
    """Semi-normal equations if A is well conditioned, QR otherwise"""
    try:
        return CSNE(A)
    except np.linalg.LinAlgError:
        pass
    return QR(A, overwrite=overwrite)


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
    F = far(np.vstack((dirs, d)))
    ext = 4 * np.pi / k * float(np.imag(np.conj(p) @ F[-1]))
    sca = float(np.sum(w * np.sum(np.abs(F[:-1]) ** 2, axis=1)))
    return ext, sca
