"""Least squares by Householder QR without forming Q.

scipy.linalg.qr(mode="economic") builds the explicit factor Q, which costs
about as much as the factorisation itself. Here the reflectors stay in the
factored form returned by LAPACK (geqrf) and are applied to the right-hand
side with unmqr, so a solve costs one factorisation; further right-hand
sides (other incident waves) cost one application each.
"""
import numpy as np
from scipy.linalg import lapack, lstsq, solve_triangular


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


def qr_solve(A, b, cond=1e-12):
    return QR(A).solve(b, cond)
