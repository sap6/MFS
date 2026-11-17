"""Least squares by Householder QR without forming Q.

scipy.linalg.qr(mode="economic") builds the explicit factor Q, which costs
about as much as the factorisation itself. Here the reflectors stay in the
factored form returned by LAPACK (geqrf) and are applied to the right-hand
side with unmqr, so a solve costs one factorisation; further right-hand
sides (other incident waves) cost one application each.
"""
import numpy as np
from scipy.linalg import blas, lapack, lstsq, solve_triangular


class QR:
    def __init__(self, A):
        """Factorises A (complex, rows >= columns). A is overwritten if it is
        Fortran-ordered; otherwise LAPACK works on a copy."""
        self.n = A.shape[1]
        lwork = int(lapack.zgeqrf_lwork(*A.shape)[0].real)          # blocked algorithm needs the optimal workspace
        self.qr, self.tau, _, info = lapack.zgeqrf(A, lwork=lwork, overwrite_a=True)
        if info != 0:
            raise np.linalg.LinAlgError(f"zgeqrf: info = {info}")
        self.R = np.triu(self.qr[:self.n])
        self.dR = np.abs(np.diag(self.R))

    def qtb(self, b):
        """Q^H b in full (the first n entries are the right-hand side of the
        triangular system, the rest is the residual in the basis of Q's
        complement); b is (rows,) or (rows, nrhs)."""
        B = np.asfortranarray(b.reshape(len(b), -1), dtype=complex)
        if B is b or np.shares_memory(B, b):
            B = B.copy(order="F")
        lwork = int(lapack.zunmqr("L", "C", self.qr, self.tau, B, -1)[1][0].real)
        y, _, info = lapack.zunmqr("L", "C", self.qr, self.tau, B, lwork, overwrite_c=True)
        if info != 0:
            raise np.linalg.LinAlgError(f"zunmqr: info = {info}")
        return y.reshape(b.shape)

    def solve(self, b, cond=1e-12, residual=False):
        """Minimiser of |A c - b| and the numerical rank (and, on request,
        the norm of the residual per right-hand side). Plain back
        substitution if R is well enough conditioned; otherwise the
        rank-revealing solve of the triangular problem (same minimiser)."""
        y = self.qtb(b)
        if self.dR.min() > cond * self.dR.max():
            c, rank = solve_triangular(self.R, y[:self.n]), self.n
        else:
            c, _, rank, _ = lstsq(self.R, y[:self.n], lapack_driver="gelsy", cond=cond)
        if residual:
            r = np.sqrt(np.linalg.norm(y[self.n:], axis=0) ** 2
                        + np.linalg.norm(self.R @ c - y[:self.n], axis=0) ** 2)
            return c, int(rank), r
        return c, int(rank)


class CSNE:
    """Least squares by the corrected semi-normal equations (Bjorck 1987):
    Cholesky factor of the Gram matrix A^H A, solution of R^H R x = A^H b,
    then `refine` steps of refinement with the residual b - A x. About half
    the work of Householder QR for tall matrices, and as accurate as long as
    cond(A)^2 * eps is well below one. The columns are scaled to unit norm
    implicitly (their norms are the diagonal of the Gram matrix), so A, a
    C-ordered array, is neither copied nor changed. Raises LinAlgError if
    the Gram matrix is not numerically positive definite or the diagonal of
    its factor suggests a condition number above about `max_cond`."""

    def __init__(self, A, max_cond=1e7, refine=2):
        self.A, self.n, self.refine = A, A.shape[1], refine
        G = np.conj(blas.zherk(1.0, A.T, trans=0, lower=0))        # upper triangle of A^H A, no copy of A
        self.scale = np.sqrt(np.real(np.diag(G)))
        G /= self.scale[:, None] * self.scale[None, :]
        c, info = lapack.zpotrf(G, lower=0, overwrite_a=True)
        if info != 0:
            raise np.linalg.LinAlgError("A^H A is not numerically positive definite")
        self.c = c
        self.dR = np.abs(np.diag(c))
        if self.dR.max() > max_cond * 1e-3 * self.dR.min():      # the diagonal ratio underestimates cond by up to ~1e3
            raise np.linalg.LinAlgError("A is too ill conditioned for the semi-normal equations")

    def _ne(self, r):
        """(A^H A)^{-1} A^H r, with the scaling undone."""
        rhs = np.conj(self.A.T @ np.conj(r))
        sc = self.scale if rhs.ndim == 1 else self.scale[:, None]
        return lapack.zpotrs(self.c, rhs / sc, lower=0)[0] / sc

    def solve(self, b, cond=None, residual=False):
        """Minimiser of |A x - b| (and, on request, the residual norms)."""
        b = np.asarray(b, dtype=complex)
        x = self._ne(b)
        for _ in range(self.refine):
            x = x + self._ne(b - self.A @ x)
        if residual:
            return x, np.linalg.norm(b - self.A @ x, axis=0)
        return x


class ScaledQR:
    """Householder QR of the column-scaled matrix, with the interface of CSNE."""

    def __init__(self, A):
        Ac = np.asfortranarray(A)
        if Ac is A or np.shares_memory(Ac, A):
            Ac = Ac.copy(order="F")
        self.scale = np.linalg.norm(Ac, axis=0)
        Ac /= self.scale
        self.qr = QR(Ac)
        self.n, self.dR = self.qr.n, self.qr.dR

    def solve(self, b, cond=1e-12, residual=False):
        out = self.qr.solve(np.asarray(b, dtype=complex), cond, residual=True)
        sc = self.scale if out[0].ndim == 1 else self.scale[:, None]
        return (out[0] / sc, out[2]) if residual else out[0] / sc


def factor(A, method="auto"):
    """Factorisation of a least-squares matrix with implicit column scaling:
    the corrected semi-normal equations if A is well conditioned, Householder
    QR otherwise (or as asked: method = "qr", "csne"). The result has
    .solve(b, residual=...) returning the unscaled minimiser, .scale (the
    column norms), .n and .dR."""
    if method != "qr":
        try:
            return CSNE(A)
        except np.linalg.LinAlgError:
            if method == "csne":
                raise
    return ScaledQR(A)


def qr_solve(A, b, cond=1e-12):
    return QR(A).solve(b, cond)
