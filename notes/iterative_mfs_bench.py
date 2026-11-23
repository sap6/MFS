"""
Matrix-free block-QR preconditioned LSMR verification on polyhedral skeletons
(Section 8 of the paper): evaluates A x and A^H v directly via the panel/moment
trees in mfs.lines without forming the full collocation matrix A, and compares
diagonal column scaling, disjoint skeleton-block QR, and coarse-bordered
skeleton-block QR on the unit cube and zero-symmetry scalene tetrahedron.

    python notes/iterative_mfs_bench.py --case cube6
    python notes/iterative_mfs_bench.py --case scalene4
    python notes/iterative_mfs_bench.py --scalene-dir 5 && python notes/iterative_mfs_bench.py --scalene-iter 5
    python notes/iterative_mfs_bench.py --summary
"""
import gc
import os
import sys
import time

import numpy as np
import scipy.linalg as sla
from scipy.linalg import blas

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO)

from mfs import core, cube, lines, tetrahedron

NPZ = os.path.join(REPO, "data", "iterative_check.npz")
D_W, P_W = core.wave()
PTS = 2.0 * core.fibonacci(400)


def lsmr_ne_step(matvec, rmatvec, b, tol_ne=1e-4, maxiter=250, state=None, step_limit=None):
    """Resumable Golub-Kahan LSMR returning (done, y, itn, rel_ne, rel_r, state_dict)."""
    bnorm = float(np.linalg.norm(b))
    if bnorm == 0.0:
        return True, np.zeros(0, np.complex128), 0, 0.0, 0.0, None
    if state is None:
        u = b.astype(np.complex128, copy=True)
        beta = bnorm
        u /= beta
        v = rmatvec(u)
        alpha = float(np.linalg.norm(v))
        v /= alpha
        ar0 = alpha * beta
        alphabar, zetabar = alpha, ar0
        rho, rhobar, cbar, sbar = 1.0, 1.0, 1.0, 0.0
        h, hbar, y = v.copy(), np.zeros_like(v), np.zeros_like(v)
        betadd, betad, rhodold, tautildeold, thetatilde, zeta, d = beta, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
        start_itn = 1
    else:
        u, v, h, hbar, y = state["u"], state["v"], state["h"], state["hbar"], state["y"]
        (alpha, beta, ar0, alphabar, zetabar, rho, rhobar, cbar, sbar,
         betadd, betad, rhodold, tautildeold, thetatilde, zeta, d) = [float(x) for x in state["sc"]]
        start_itn = int(state["next_itn"])

    end_itn = maxiter if step_limit is None else min(maxiter, start_itn + step_limit - 1)
    rnorm = float(np.sqrt(d))
    for itn in range(start_itn, end_itn + 1):
        u = matvec(v) - alpha * u
        beta = float(np.linalg.norm(u))
        if beta > 0.0:
            u /= beta
            v = rmatvec(u) - beta * v
            alpha = float(np.linalg.norm(v))
            if alpha > 0.0:
                v /= alpha

        chat, shat, alphahat = _sym_ortho(alphabar, 0.0)
        rhoold = rho
        c, s, rho = _sym_ortho(alphahat, beta)
        thetanew = s * alpha
        alphabar = c * alpha

        rhobarold = rhobar
        zetaold = zeta
        thetabar = sbar * rho
        rhotemp = cbar * rho
        cbar, sbar, rhobar = _sym_ortho(rhotemp, thetanew)
        zeta = cbar * zetabar
        zetabar = -sbar * zetabar

        hbar = h - (thetabar * rho / (rhoold * rhobarold)) * hbar
        y += (zeta / (rho * rhobar)) * hbar
        h = v - (thetanew / rho) * h

        betaacute = chat * betadd
        betacheck = -shat * betadd
        betahat = c * betaacute
        betadd = -s * betaacute

        thetatildeold = thetatilde
        ctildeold, stildeold, rhotildeold = _sym_ortho(rhodold, thetabar)
        thetatilde = stildeold * rhobar
        rhodold = ctildeold * rhobar
        betad = -stildeold * betad + ctildeold * betahat

        tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold
        taud = (zeta - thetatilde * tautildeold) / rhodold
        d += betacheck * betacheck
        rnorm = float(np.sqrt(d + (betad - taud) ** 2 + betadd * betadd))
        arnorm = abs(zetabar)
        if arnorm / ar0 <= tol_ne:
            return True, y, itn, arnorm / ar0, rnorm / bnorm, None
    if end_itn >= maxiter:
        return True, y, maxiter, abs(zetabar) / ar0, rnorm / bnorm, None
    st_out = {
        "u": u, "v": v, "h": h, "hbar": hbar, "y": y, "next_itn": np.array(end_itn + 1, np.int64),
        "sc": np.array([alpha, beta, ar0, alphabar, zetabar, rho, rhobar, cbar, sbar,
                        betadd, betad, rhodold, tautildeold, thetatilde, zeta, d], float),
    }
    return False, y, end_itn, abs(zetabar) / ar0, rnorm / bnorm, st_out


def lsmr_ne(matvec, rmatvec, b, tol_ne=1e-4, maxiter=250):
    """
    Fong & Saunders (2011) Golub-Kahan LSMR for min_y ||b - A y||_2, stopping cleanly when
    the relative normal-equation residual ||A^H r_k||_2 / ||A^H b||_2 <= tol_ne.
    Identical iterates to scipy.sparse.linalg.lsmr (verified in tests/test_mfs.py).
    """
    _, y, itn, rel_ne, rel_r, _ = lsmr_ne_step(matvec, rmatvec, b, tol_ne=tol_ne, maxiter=maxiter)
    return y, itn, rel_ne, rel_r


def _sym_ortho(a, b):
    if b == 0.0:
        return (1.0 if a >= 0 else -1.0), 0.0, abs(a)
    if a == 0.0:
        return 0.0, (1.0 if b >= 0 else -1.0), abs(b)
    if abs(b) > abs(a):
        tau = a / b
        s = (1.0 if b >= 0 else -1.0) / float(np.sqrt(1.0 + tau * tau))
        return s * tau, s, b / s
    tau = b / a
    c = (1.0 if a >= 0 else -1.0) / float(np.sqrt(1.0 + tau * tau))
    return c, c * tau, a / c


# --- Matrix-free forward & adjoint operators and block-by-block preconditioner builders ---

def cube_mf_ops(cb, c):
    """Matrix-free A x and A^H v for parity class c of Cube cb (never forms A)"""
    Tw = np.ascontiguousarray(cb.T * cb.rw[:, None, None])
    sgn = np.array([[core.sign(c, R) for R in cube.REFL]])
    kp = cb.keep[c]
    line_entries = [(i, ln) for i, ln in enumerate(cb.lines) if isinstance(ln, cube.Line) and ln.kind == "e"]

    def matvec(x_vec):
        E = cb.field({c: x_vec}, cb.x)
        return np.einsum("mac,mc->ma", Tw, E).ravel()

    def rmatvec(v_vec):
        u = np.ascontiguousarray(np.einsum("mac,ma->mc", Tw, np.conj(v_vec.reshape(-1, 2)))[None, :, :])
        ne, nm = len(cb.pe), len(cb.pm)
        n_vt = 9 if cb.tails else 0
        y_pt = np.ascontiguousarray(np.vstack((cb.pe, cb.lines[-2].pts)) if n_vt else cb.pe)
        out_a = np.zeros((1, len(y_pt), 3), np.complex128)
        out_b = np.zeros((1, len(y_pt), 3), np.complex128)
        lines.eval_points_adj(cb.x, cube.REFL, sgn, u, y_pt, float(cb.k), out_a, out_b)

        c_all = np.zeros(3 * cb.slot0[-1], np.complex128)
        c_all[:3 * ne] = out_a[0, :ne].ravel()
        if nm:
            c_all[3 * ne:3 * (ne + nm)] = out_b[0, ne - nm:ne].ravel()
        if n_vt:
            vt_e, vt_m = cb.lines[-2], cb.lines[-1]
            lo = 3 * cb.slot0[-3]
            c_all[lo:lo + 3] = vt_e.w @ out_a[0, ne:]
            c_all[lo + 3:lo + 6] = vt_m.w @ out_b[0, ne:]

        for i, ln in line_entries:
            lo = 3 * cb.slot0[1 + i]
            out_pair = np.zeros((1, 2, 1, ln.n, 3), np.complex128)
            w = None if ln.w is None else np.stack((ln.w, cb.lines[i + 1].w))
            lines.eval_line_adj(cb.x, cube.REFL, ln.r_act, sgn, ln.t_hat, ln.mids, ln.lam, ln.d, ln.n,
                                float(cb.k), cube.H, u, out_pair, rot=ln.rot, w=w)
            c_all[lo:lo + 6 * ln.n] = out_pair.reshape(6 * ln.n)
        return c_all[kp]

    return matvec, rmatvec


def cube_assemble_cols(cb, c, glob_cols):
    """Assembles ONLY the selected global columns of class c (never allocates the full matrix A)"""
    Tw = np.ascontiguousarray(cb.T * cb.rw[:, None, None])
    M = cb.rows
    nb = len(glob_cols)
    cmap = np.full((1, 3 * cb.slot0[-1]), -1, np.int64)
    cmap[0, glob_cols] = np.arange(nb)
    A = np.moveaxis(np.zeros((M, nb, 1), np.complex128, order="F"), 2, 0)
    sgn = np.array([[core.sign(c, R) for R in cube.REFL]])
    if np.any(cmap[0, :3 * cb.slot0[1]] >= 0):
        # Assemble points only if present in glob_cols
        A_pt = np.moveaxis(np.zeros((M, 3 * cb.slot0[1], 1), np.complex128, order="F"), 2, 0)
        lines.points_add(cb.x, Tw, cube.REFL, cb.pe, cb.pm, cb.k, A_pt, np.zeros(1, np.int64), sgn)
        act = np.nonzero(cmap[0, :3 * cb.slot0[1]] >= 0)[0]
        A[0][:, cmap[0, act]] = A_pt[0][:, act]
        del A_pt
    for i, ln in enumerate(cb.lines):
        if isinstance(ln, cube.Line) and ln.kind == "e":
            lo = 3 * cb.slot0[1 + i]
            sub = cmap[:, lo:lo + 6 * ln.n]
            if np.any(sub >= 0):
                w = None if ln.w is None else np.stack((ln.w, cb.lines[i + 1].w))
                lines.line_add(cb.x, Tw, cube.REFL, ln.r_act, ln.t_hat, ln.mids, ln.lam, ln.d, ln.n, cb.k, cube.H,
                               A, sub, sgn, rot=ln.rot, w=w)
    if cb.tails:
        lo = 3 * cb.slot0[-3]
        if np.any(cmap[0, lo:lo + 6] >= 0):
            vt_e, vt_m = cb.lines[-2], cb.lines[-1]
            A_vt = np.moveaxis(np.zeros((M, 6, 1), np.complex128, order="F"), 2, 0)
            lines.radial_tails_add(cb.x, Tw, cube.REFL, np.ascontiguousarray(vt_e.pts[None]),
                                   np.ascontiguousarray(vt_e.w[None]), np.ascontiguousarray(vt_m.w[None]),
                                   cb.k, A_vt, np.zeros(1, np.int64), np.array([3], np.int64), sgn)
            act = np.nonzero(cmap[0, lo:lo + 6] >= 0)[0]
            A[0][:, cmap[0, lo + act]] = A_vt[0][:, act]
    return A[0]


def scalene_mf_ops(pr):
    """Matrix-free A x and A^H v for zero-symmetry ScaleneTetrahedron pr (never forms A)"""
    Tw = np.ascontiguousarray(pr.T * pr.rw[:, None, None])

    def matvec(x_vec):
        return np.einsum("mac,mc->ma", Tw, pr.field(x_vec, pr.x)).ravel()

    def rmatvec(v_vec):
        u = np.ascontiguousarray(np.einsum("mac,ma->mc", Tw, np.conj(v_vec.reshape(-1, 2)))[None, :, :])
        ne, nm = len(pr.pe), len(pr.pm)
        y_pt = np.ascontiguousarray(np.vstack([pr.pe] + [rpts.reshape(-1, 3) for rpts, *_ in pr.v_tails]))
        out_a = np.zeros((1, len(y_pt), 3), np.complex128)
        out_b = np.zeros((1, len(y_pt), 3), np.complex128)
        lines.eval_points_adj(pr.x, tetrahedron._REFL1, tetrahedron._SGN1, u, y_pt, pr.k, out_a, out_b)

        z = np.zeros(pr.n_unknowns, np.complex128)
        z[:3 * ne] = out_a[0, :ne].ravel()
        if nm:
            z[3 * ne:3 * (ne + nm)] = out_b[0, ne - nm:ne].ravel()
        off_y, c_vt = ne, pr.vt_col0
        for rpts, we, wm, nc_e, nc_m in pr.v_tails:
            R_cnt, nstep = rpts.shape[0], rpts.shape[1]
            n_r = R_cnt * nstep
            z[c_vt:c_vt + nc_e] = np.einsum("et,rtc->rec", we, out_a[0, off_y:off_y + n_r].reshape(R_cnt, nstep, 3)).ravel()
            z[c_vt + nc_e:c_vt + nc_e + nc_m] = np.einsum("mt,rtc->rmc", wm, out_b[0, off_y:off_y + n_r].reshape(R_cnt, nstep, 3)).ravel()
            off_y += n_r
            c_vt += nc_e + nc_m

        for (H_e, t_hat, mids, lam, dp, n, F, w_pair, is_tail), cm in zip(pr.lines, pr.line_cmaps):
            out_pair = np.zeros((1, 2, 1, n, 3), np.complex128)
            lines.eval_line_adj(pr.x, tetrahedron._REFL1, tetrahedron._RACT1, tetrahedron._SGN1,
                                t_hat, mids, lam, dp, n, pr.k, H_e, u, out_pair,
                                rot=F, w=w_pair if is_tail else None)
            act = cm[0] >= 0
            z[cm[0, act]] += out_pair.reshape(6 * n)[act]
        return z

    return matvec, rmatvec


def scalene_assemble_cols(pr, cols):
    """Assembles ONLY the selected columns of ScaleneTetrahedron pr (never allocates full A)"""
    Tw = np.ascontiguousarray(pr.T * pr.rw[:, None, None])
    M = pr.rows
    nb = len(cols)
    g2b = np.full(pr.n_unknowns, -1, np.int64)
    g2b[cols] = np.arange(nb)
    A = np.moveaxis(np.zeros((M, nb, 1), np.complex128, order="F"), 2, 0)
    n_pt_cols = 3 * (len(pr.pe) + len(pr.pm))
    if np.any(g2b[:n_pt_cols] >= 0):
        A_pt = np.moveaxis(np.zeros((M, n_pt_cols, 1), np.complex128, order="F"), 2, 0)
        lines.points_add(pr.x, Tw, tetrahedron._REFL1, pr.pe, pr.pm, pr.k, A_pt, np.zeros(1, np.int64), tetrahedron._SGN1)
        act = np.nonzero(g2b[:n_pt_cols] >= 0)[0]
        A[0][:, g2b[act]] = A_pt[0][:, act]
        del A_pt
    for (H_e, t_hat, mids, lam, dp, n, F, w_pair, is_tail), cm in zip(pr.lines, pr.line_cmaps):
        cm_sub = np.where(cm >= 0, g2b[np.maximum(0, cm)], -1)
        if np.any(cm_sub >= 0):
            lines.line_add(pr.x, Tw, tetrahedron._REFL1, tetrahedron._RACT1, t_hat, mids, lam, dp, n, pr.k, H_e,
                           A, cm_sub, tetrahedron._SGN1, rot=F, w=w_pair if is_tail else None)
    c_vt = pr.vt_col0
    for rpts, we, wm, nc_e, nc_m in pr.v_tails:
        tot_v = nc_e + nc_m
        if np.any(g2b[c_vt:c_vt + tot_v] >= 0):
            A_vt = np.moveaxis(np.zeros((M, tot_v, 1), np.complex128, order="F"), 2, 0)
            lines.radial_tails_add(pr.x, Tw, tetrahedron._REFL1, rpts, we, wm, pr.k, A_vt,
                                   np.zeros(1, np.int64), np.array([nc_e], np.int64), tetrahedron._SGN1)
            act = np.nonzero(g2b[c_vt:c_vt + tot_v] >= 0)[0]
            A[0][:, g2b[c_vt + act]] = A_vt[0][:, act]
            del A_vt
        c_vt += tot_v
    return A[0]


def factor_bordered_qr(assemble_fn, idx0, idxe):
    N = len(idx0) + sum(len(ix) for ix in idxe)
    diag_norm = np.empty(N)
    A0 = assemble_fn(idx0)
    sc0 = np.linalg.norm(A0, axis=0)
    diag_norm[idx0] = sc0
    A0 /= sc0
    Q0, R0 = np.linalg.qr(A0, mode="reduced")
    del A0
    n0 = len(idx0)
    p_bytes = (n0 * (n0 + 1) // 2) * 16 + sc0.nbytes

    Bs, Rs, scs = [], [], []
    for idx in idxe:
        Ae = assemble_fn(idx)
        sc = np.linalg.norm(Ae, axis=0)
        diag_norm[idx] = sc
        Ae /= sc
        Be = Q0.conj().T @ Ae
        Ae -= Q0 @ Be
        Re = np.linalg.qr(Ae, mode="r")
        del Ae
        Bs.append(Be)
        Rs.append(Re)
        scs.append(sc)
        nb = len(idx)
        p_bytes += Be.nbytes + (nb * (nb + 1) // 2) * 16 + sc.nbytes
    del Q0
    gc.collect()
    return R0, sc0, Bs, Rs, scs, p_bytes, diag_norm


def make_bordered_solvers(idx0, idxe, R0, sc0, Bs, Rs, scs):
    def p_inv(y):
        x = np.empty_like(y)
        r0 = y[idx0].copy()
        for idx, sc, B, R in zip(idxe, scs, Bs, Rs):
            ue = sla.solve_triangular(R, y[idx], lower=False, check_finite=False)
            r0 -= B @ ue
            x[idx] = ue / sc
        x[idx0] = sla.solve_triangular(R0, r0, lower=False, check_finite=False) / sc0
        return x

    def p_invH(z):
        out = np.empty_like(z)
        u0 = sla.solve_triangular(R0, z[idx0] / sc0, trans="C", lower=False, check_finite=False)
        out[idx0] = u0
        for idx, sc, B, R in zip(idxe, scs, Bs, Rs):
            out[idx] = sla.solve_triangular(R, z[idx] / sc - B.conj().T @ u0, trans="C", lower=False, check_finite=False)
        return out

    return p_inv, p_invH


def build_bordered_qr(assemble_fn, idx0, idxe):
    """
    Coarse-bordered skeleton block-QR preconditioner assembled one block at a time (never forming full A):
      Block 0 (coarse core) -> Q0 R0; each edge wall e -> (Ae - Q0 Be) = Qe Re, then Q0 is freed!
      Returns (p_inv, p_invH, p_bytes, diag_norm)
    """
    R0, sc0, Bs, Rs, scs, p_bytes, diag_norm = factor_bordered_qr(assemble_fn, idx0, idxe)
    p_inv, p_invH = make_bordered_solvers(idx0, idxe, R0, sc0, Bs, Rs, scs)
    return p_inv, p_invH, p_bytes, diag_norm


def load_state():
    return dict(np.load(NPZ)) if os.path.exists(NPZ) else {}


SPLIT_TMP = "/tmp/iterative_mfs_split.npz"


def run_cube_one(J, k_mult=1, part="all"):
    st = load_state()
    cb_w = cube.Cube(2, classes=[])
    mv_w, rmv_w = cube_mf_ops(cb_w, (-1, 1, 1))
    _ = mv_w(np.ones(cb_w.unknowns[(-1, 1, 1)], np.complex128))
    _ = rmv_w(np.ones(cb_w.rows, np.complex128))
    del cb_w

    k_val = k_mult * 2 * np.pi
    c = (-1, 1, 1)
    cb = cube.Cube(J, k=k_val, classes=[])
    N, M = cb.unknowns[c], cb.rows
    dense_mb = 16.0 * M * N / 2 ** 20
    kp = np.nonzero(cb.keep[c])[0]
    g2l = np.full(3 * cb.slot0[-1], -1, np.int64)
    g2l[kp] = np.arange(N)

    b = -np.einsum("mab,mb->ma", cb.T, core.class_part(cb.x, cb.k, D_W, P_W, cube.REFL, c)).reshape(-1) * np.repeat(cb.rw, 2)
    mv, rmv = cube_mf_ops(cb, c)

    if part in ("all", "blk"):
        cb_dir = cube.Cube(J, k=k_val, classes=[c])
        x_dir, _ = cb_dir.fac[c].solve(b, residual=True)
        E_dir = cb.field({c: x_dir}, PTS)
        del cb_dir
        gc.collect()

        coarse_glob = [np.arange(3 * cb.slot0[1]), np.arange(3 * cb.slot0[-3], 3 * cb.slot0[-1])]
        edge_globs = [[] for _ in range(3)]
        for i, ln in enumerate(cb.lines):
            if isinstance(ln, cube.Line) and ln.kind == "e":
                gl = np.arange(3 * cb.slot0[1 + i], 3 * cb.slot0[1 + i] + 6 * ln.n)
                if not ln.tail and ln.d >= cb.depths[1] - 1e-12:
                    coarse_glob.append(gl)
                else:
                    edge_globs[ln.ax].append(gl)
        cg = np.concatenate(coarse_glob)
        idx0 = g2l[cg[g2l[cg] >= 0]]
        idxe = [g2l[np.concatenate(eg)[g2l[np.concatenate(eg)] >= 0]] for eg in edge_globs]

        t0 = time.time()
        p_inv, p_invH, p_bytes, diag_sc = build_bordered_qr(lambda loc: cube_assemble_cols(cb, c, kp[loc]), idx0, idxe)
        t_prec = time.time() - t0

        t0 = time.time()
        y_blk, itn_blk, ne_blk, r_blk = lsmr_ne(lambda v: mv(p_inv(v)), lambda w: p_invH(rmv(w)), b, tol_ne=1e-4, maxiter=120)
        t_iter = time.time() - t0
        E_blk = cb.field({c: p_inv(y_blk)}, PTS)
        rel_blk = core.rel(E_blk, E_dir)
        p_mb = p_bytes / 2 ** 20
        if part == "blk":
            np.savez(SPLIT_TMP, E_dir=E_dir, diag_sc=diag_sc,
                     blk=np.array([p_mb, itn_blk, t_prec, t_iter, rel_blk], float))
            print(f"Cube k={k_mult}*2pi J={J} [blk]: P_blk={p_mb:.2f} MiB, itn={itn_blk:2d} "
                  f"(prec={t_prec:.2f}s, iter={t_iter:.2f}s), rel_E={rel_blk:.2e}", flush=True)
            return
    else:
        sp = np.load(SPLIT_TMP)
        E_dir, diag_sc = sp["E_dir"], sp["diag_sc"]
        p_mb, itn_blk, t_prec, t_iter, rel_blk = sp["blk"]
        itn_blk = int(itn_blk)

    y_dg, itn_dg, ne_dg, r_dg = lsmr_ne(lambda v: mv(v / diag_sc), lambda w: rmv(w) / diag_sc, b,
                                        tol_ne=1e-4, maxiter=60 if J >= 8 else 100)
    rel_dg = core.rel(cb.field({c: y_dg / diag_sc}, PTS), E_dir)

    print(f"Cube k={k_mult}*2pi J={J}: M={M}, N={N}, A_dense={dense_mb:.1f} MiB vs P_blk={p_mb:.2f} MiB ({dense_mb/p_mb:.1f}x) | "
          f"Blk-LSMR: itn={itn_blk:2d} (prec={t_prec:.2f}s, iter={t_iter:.2f}s), rel_E={rel_blk:.2e} | "
          f"Diag-LSMR: itn={itn_dg:3d}, rel_E={rel_dg:.2e}", flush=True)

    prev = {(int(r[1]), int(r[0])): list(r) for r in st.get("cube", np.zeros((0, 12)))}
    prev[k_mult, J] = [J, k_mult, M, N, dense_mb, p_mb, itn_blk, t_prec, t_iter, rel_blk, itn_dg, rel_dg]
    st["cube"] = np.array([prev[k] for k in sorted(prev)])
    np.savez(NPZ, **st)
    if part == "diag" and os.path.exists(SPLIT_TMP):
        os.remove(SPLIT_TMP)


def _scalene_blocks(pr):
    n_pt = 3 * (len(pr.pe) + len(pr.pm))
    vt_cols = np.arange(pr.vt_col0, pr.n_unknowns)
    coarse_list = [np.arange(n_pt), vt_cols]
    edge_lists = [[] for _ in range(6)]
    n_depths = len(pr.depths)
    for d_i in range(n_depths):
        for e_i in range(6):
            cm = pr.line_cmaps[d_i * 6 + e_i]
            cols = cm[0, cm[0] >= 0]
            if d_i < 2:
                coarse_list.append(cols)
            else:
                edge_lists[e_i].append(cols)
    for e_i in range(6):
        cm = pr.line_cmaps[n_depths * 6 + e_i]
        edge_lists[e_i].append(cm[0, cm[0] >= 0])
    return np.concatenate(coarse_list), [np.concatenate(el) for el in edge_lists]


def run_scalene(J, part="all", step_limit=30):
    st = load_state()
    tmp_path = f"/tmp/iterative_mfs_scalene_J{J}.npz"
    if part == "dir":
        E_dir, t_dir, _ = tetrahedron.solve_scalene(J, D_W, P_W, PTS, q=2 ** -0.5)
        np.savez(tmp_path, E_dir=E_dir, t_dir=float(t_dir))
        print(f"Scalene J={J} [dir]: t_dir={t_dir:.2f}s", flush=True)
        return

    pr_w = tetrahedron.ScaleneTetrahedron(2, q=2 ** -0.5)
    mv_w, rmv_w = scalene_mf_ops(pr_w)
    _ = mv_w(np.ones(pr_w.n_unknowns, np.complex128))
    _ = rmv_w(np.ones(pr_w.rows, np.complex128))
    del pr_w

    pr = tetrahedron.ScaleneTetrahedron(J, q=2 ** -0.5)
    M, N = pr.rows, pr.n_unknowns
    dense_mb = 16.0 * M * N / 2 ** 20
    b = -np.einsum("mab,mb->ma", pr.T, core.plane_wave(pr.x, pr.k, D_W, P_W)).reshape(-1) * np.repeat(pr.rw, 2)
    idx0, idxe = _scalene_blocks(pr)

    if part == "prec":
        sp = dict(np.load(tmp_path))
        t0 = time.time()
        R0, sc0, Bs, Rs, scs, p_bytes, _ = factor_bordered_qr(lambda cols: scalene_assemble_cols(pr, cols), idx0, idxe)
        t_prec = time.time() - t0
        p_mb = p_bytes / 2 ** 20
        sp.update({"R0": R0, "sc0": sc0, "p_mb": float(p_mb), "t_prec": float(t_prec), "t_iter": 0.0, "has_lsmr": np.array(0)})
        for i, (B, R, sc) in enumerate(zip(Bs, Rs, scs)):
            sp[f"eB{i}"], sp[f"eR{i}"], sp[f"esc{i}"] = B, R, sc
        np.savez(tmp_path, **sp)
        print(f"Scalene J={J} [prec]: P_blk={p_mb:.1f} MiB, t_prec={t_prec:.2f}s", flush=True)
        return

    if part == "all":
        E_dir, t_dir, _ = tetrahedron.solve_scalene(J, D_W, P_W, PTS, q=2 ** -0.5)
        gc.collect()
        t0 = time.time()
        p_inv, p_invH, p_bytes, _ = build_bordered_qr(lambda cols: scalene_assemble_cols(pr, cols), idx0, idxe)
        t_prec = time.time() - t0
        p_mb = p_bytes / 2 ** 20
        mv, rmv = scalene_mf_ops(pr)
        t0 = time.time()
        y_blk, itn_blk, ne_blk, r_blk = lsmr_ne(lambda v: mv(p_inv(v)), lambda w: p_invH(rmv(w)), b, tol_ne=1e-4, maxiter=80)
        t_iter = time.time() - t0
    else:
        sp = dict(np.load(tmp_path))
        E_dir, t_dir = sp["E_dir"], float(sp["t_dir"])
        if "eR0" not in sp:
            t0 = time.time()
            R0, sc0, Bs, Rs, scs, p_bytes, _ = factor_bordered_qr(lambda cols: scalene_assemble_cols(pr, cols), idx0, idxe)
            t_prec = time.time() - t0
            p_mb = p_bytes / 2 ** 20
            t_iter_prev, lsmr_st = 0.0, None
        else:
            R0, sc0, p_mb, t_prec = sp["R0"], sp["sc0"], float(sp["p_mb"]), float(sp["t_prec"])
            Bs = [sp[f"eB{i}"] for i in range(6)]
            Rs = [sp[f"eR{i}"] for i in range(6)]
            scs = [sp[f"esc{i}"] for i in range(6)]
            t_iter_prev = float(sp.get("t_iter", 0.0))
            lsmr_st = ({k[5:]: sp[k] for k in sp if k.startswith("lsmr_")}
                       if int(sp.get("has_lsmr", 0)) else None)
        p_inv, p_invH = make_bordered_solvers(idx0, idxe, R0, sc0, Bs, Rs, scs)
        mv, rmv = scalene_mf_ops(pr)
        t0 = time.time()
        done, y_blk, itn_blk, ne_blk, r_blk, st_next = lsmr_ne_step(
            lambda v: mv(p_inv(v)), lambda w: p_invH(rmv(w)), b,
            tol_ne=1e-4, maxiter=80, state=lsmr_st, step_limit=step_limit)
        t_iter = t_iter_prev + (time.time() - t0)
        if not done:
            sp.update({"R0": R0, "sc0": sc0, "p_mb": float(p_mb), "t_prec": float(t_prec),
                       "t_iter": float(t_iter), "has_lsmr": np.array(1)})
            for i, (B, R, sc) in enumerate(zip(Bs, Rs, scs)):
                sp[f"eB{i}"], sp[f"eR{i}"], sp[f"esc{i}"] = B, R, sc
            for k, v in st_next.items():
                sp[f"lsmr_{k}"] = v
            np.savez(tmp_path, **sp)
            print(f"Scalene J={J} [iter {itn_blk}]: NE_res={ne_blk:.2e}, t_iter={t_iter:.2f}s (checkpoint saved)", flush=True)
            return

    x_blk = p_inv(y_blk)
    E_blk = pr.field(x_blk, PTS)
    rel_blk = core.rel(E_blk, E_dir)

    print(f"Scalene J={J}: M={M}, N={N}, A_dense={dense_mb:.1f} MiB vs P_blk={p_mb:.1f} MiB ({dense_mb/p_mb:.1f}x) | "
          f"Blk-LSMR: itn={itn_blk:2d} (prec={t_prec:.2f}s, iter={t_iter:.2f}s vs dir={t_dir:.2f}s), "
          f"NE_res={ne_blk:.2e}, rel_E_vs_dir={rel_blk:.2e}", flush=True)

    prev = {int(r[0]): list(r) for r in st.get("scalene", np.zeros((0, 9)))}
    prev[J] = [J, M, N, dense_mb, p_mb, itn_blk, t_prec, t_iter, rel_blk]
    st["scalene"] = np.array([prev[k] for k in sorted(prev)])
    np.savez(NPZ, **st)
    if os.path.exists(tmp_path):
        os.remove(tmp_path)


def show_summary():
    st = load_state()
    if "cube" in st:
        print("=== Cube Matrix-Free Block-QR LSMR (class (-1,1,1), tol_NE = 1e-4) ===")
        for r in st["cube"]:
            print(f"  J={int(r[0])} (k={int(r[1])}*2pi): M={int(r[2]):5d}, N={int(r[3]):4d}, A={r[4]:6.1f} MiB vs P={r[5]:5.2f} MiB ({r[4]/r[5]:4.1f}x) | "
                  f"Blk itn={int(r[6]):2d} (prec={r[7]:.2f}s, iter={r[8]:.2f}s, rel_E={r[9]:.2e}) vs Diag itn={int(r[10]):3d} (rel_E={r[11]:.2e})")
    if "scalene" in st:
        print("\n=== Zero-Symmetry Scalene Tetrahedron Matrix-Free Block-QR LSMR (tol_NE = 1e-4) ===")
        for r in st["scalene"]:
            print(f"  J={int(r[0])}: M={int(r[1]):5d}, N={int(r[2]):4d}, A={r[3]:6.1f} MiB vs P={r[4]:5.1f} MiB ({r[3]/r[4]:4.1f}x) | "
                  f"Blk itn={int(r[5]):2d} (prec={r[6]:.2f}s, iter={r[7]:.2f}s, rel_E_vs_dir={r[8]:.2e})")


if __name__ == "__main__":
    part = sys.argv[sys.argv.index("--part") + 1] if "--part" in sys.argv else "all"
    step = int(sys.argv[sys.argv.index("--step") + 1]) if "--step" in sys.argv else 30
    if "--scalene-dir" in sys.argv:
        run_scalene(int(sys.argv[sys.argv.index("--scalene-dir") + 1]), part="dir")
    elif "--scalene-prec" in sys.argv:
        run_scalene(int(sys.argv[sys.argv.index("--scalene-prec") + 1]), part="prec")
    elif "--scalene-iter" in sys.argv:
        run_scalene(int(sys.argv[sys.argv.index("--scalene-iter") + 1]), part="iter", step_limit=step)
    elif "--case" in sys.argv:
        cs = sys.argv[sys.argv.index("--case") + 1]
        cases = {"cube4": (4, 1), "cube6": (6, 1), "cube8": (8, 1), "cube6_k2": (6, 2)}
        if cs in cases:
            run_cube_one(*cases[cs], part=part)
        elif cs.startswith("scalene"):
            run_scalene(int(cs.replace("scalene", "")), part=part, step_limit=step)
    else:
        show_summary()
