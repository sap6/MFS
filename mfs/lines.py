"""
Fields of the line and point sources of the polyhedral solvers (Section 5 and
Appendix B of the paper), assembled straight into the matrices of the parity
classes, together with direct near- and far-field evaluation from solved modes.

A line at depth d carries a dipole density that is a polynomial in the stretched
variable u in [-1, 1], with t(u) = H tanh(λu/2) along the edge and λ = ln((2H - d)/d).
Each matrix entry is the integral of a dipole field against one Legendre mode, and
most of the work is shared between targets:
  * [-1, 1] is cut into a binary tree of panels whose nodes and weights against all
    modes are computed once per line (exact moments in t on the four coarsest levels,
    so a distant target needs only PM kernel values for the whole line);
  * a target descends the tree until a panel is no longer than its distance from it;
  * the electric and magnetic line share the kernel values, the nine lines of a tail
    are summed with their weights, and the symmetry images of a target are summed
    with the signs of each parity class before anything is written;
  * a line is mapped to itself by the midpoint reflection or half-turn, so half of
    the symmetry images repeat the other half (mode n times (-1)^n) and are not computed.
"""
import numpy as np
from numba import njit, prange
from scipy.interpolate import BarycentricInterpolator

LMOM = 4          # the four coarsest levels use moments (panel widths 2, 1, 1/2, 1/4 in u)
PM, NG = 12, 8    # nodes of a moment panel and of a Gauss panel, for a line of one wavelength
LHARD = 40        # deepest level of the descent

_TREES = {}       # the tables of each line, built once by build_tree
_INTERP = {}      # cached 1D barycentric interpolators on the pm reference Gauss nodes


def nodes(k, H):
    """
    Nodes per panel (moment panels, Gauss panels) for wavenumber k and lines of length up to 2H.
    The kernel phase along the coarsest panels grows like kH: 12 / 8 nodes at one wavelength
    (the cube fields change by < 1e-6 with more), 20 / 10 at two and 28 / 12 at three (Appendix B)
    """
    s = max(0.0, k * 2 * H / (2 * np.pi) - 1.0)          # wavelengths per edge, minus one
    return int(np.ceil(PM + 8 * s)), int(np.ceil(NG + 2 * s))


def build_tree(lam, d, n, H, k):
    """
    Node positions t, weights W (nodes, n) and the table of panels of one line, level by level;
    tabulated down to panels of length about d in the middle of the line
    """
    ltab = int(min(max(LMOM + 2, np.ceil(np.log2(lam / (2.0 * d)))), 12))
    pm, ng = nodes(k, H)
    key = (round(lam, 12), round(d, 15), n, ltab, pm, ng)
    if key in _TREES:
        return _TREES[key]
    leg = np.polynomial.legendre
    if pm not in _INTERP:
        gm, _ = leg.leggauss(pm)
        _INTERP[pm] = (gm, BarycentricInterpolator(gm, np.eye(pm)))
    gm, interp = _INTERP[pm]
    gg, gwg = leg.leggauss(ng)
    gf, gwf = leg.leggauss(16)                           # fine rule for the exact moments

    def u_rule(lo, hi, g, gw):                           # Gauss rule on [lo, hi] in u: nodes t and mode weights (scaled by 1/d)
        u = 0.5 * (lo + hi)[:, None] + 0.5 * (hi - lo)[:, None] * g
        w = 0.5 * (hi - lo)[:, None] * gw * H * 0.5 * lam / np.cosh(0.5 * lam * u) ** 2 / d
        return (H * np.tanh(0.5 * lam * u)).ravel(), leg.legvander(u.ravel(), n - 1) * w.ravel()[:, None]

    npan = 2 ** (ltab + 1) - 1
    off = np.zeros(npan + 1, np.int64)
    tlo, thi = np.empty(npan), np.empty(npan)
    ts, Ws = [], []
    p = 0                                                # index of the first panel of the level
    for lev in range(ltab + 1):
        m = 2 ** lev
        ulo = -1.0 + 2.0 * np.arange(m) / m
        uhi = ulo + 2.0 / m
        a, b = H * np.tanh(0.5 * lam * ulo), H * np.tanh(0.5 * lam * uhi)
        tlo[p:p + m], thi[p:p + m] = a, b
        if lev < LMOM:
            # Moment panels: pm Gauss nodes in t, weights = moments of the Lagrange polynomials against the modes
            nsub = max(8, int(np.ceil(2.0 / m / 0.03125)))   # sub-intervals of width at most 1/32 in u
            nf = nsub * 16
            s_grid = np.linspace(0.0, 1.0, nsub + 1)
            e_lo = ulo[:, None] + (uhi - ulo)[:, None] * s_grid[:-1]
            e_hi = ulo[:, None] + (uhi - ulo)[:, None] * s_grid[1:]
            tf_all, VW_all = u_rule(e_lo.ravel(), e_hi.ravel(), gf, gwf)
            tau = (2.0 * tf_all.reshape(m, nf) - (a + b)[:, None]) / (b - a)[:, None]
            L = interp(tau.ravel()).reshape(m, nf, pm)
            ts.append((0.5 * (a + b)[:, None] + 0.5 * (b - a)[:, None] * gm[None, :]).ravel())
            Ws.append((L.swapaxes(1, 2) @ VW_all.reshape(m, nf, n)).reshape(m * pm, n))
        else:
            # Gauss panels: ng nodes in u, weights = the modes times (dt/du) / d (column scaling cancels 1/d)
            tg, VW = u_rule(ulo, uhi, gg, gwg)
            ts.append(tg)
            Ws.append(VW)
        off[p + 1:p + m + 1] = off[p] + (pm if lev < LMOM else ng) * np.arange(1, m + 1)
        p += m
    out = (np.ascontiguousarray(np.concatenate(ts)), np.ascontiguousarray(np.vstack(Ws)), off, tlo, thi, ltab, gg, gwg)
    if len(_TREES) >= 32:
        _TREES.clear()
    _TREES[key] = out
    return out


def _mode_strides(cmap, ns, n):
    """Active mode start k0 and stride kstep (1, or 2 for single parity; k0 = n skips unused components)"""
    k0 = np.full((2, ns, 3), n, np.int64)
    kstep = np.ones((2, ns, 3), np.int64)
    ev = np.arange(0, n, 2) * 3
    od = np.arange(1, n, 2) * 3
    for kind in range(2):
        for s in range(ns):
            base = (kind * ns + s) * n * 3
            for cc in range(3):
                has_ev = bool(np.any(cmap[:, base + ev + cc] >= 0))
                has_od = bool(len(od) and np.any(cmap[:, base + od + cc] >= 0))
                if has_ev and has_od:
                    k0[kind, s, cc], kstep[kind, s, cc] = 0, 1
                elif has_ev:
                    k0[kind, s, cc], kstep[kind, s, cc] = 0, 2
                elif has_od:
                    k0[kind, s, cc], kstep[kind, s, cc] = 1, 2
    ncomp = 2 if np.all(k0[:, :, 2] == n) else 3
    return k0, kstep, ncomp


def _pack_cmap(cmap):
    """Compact active local/global column indices per class (at most half the columns are active in any mirror class)"""
    ncls = cmap.shape[0]
    nact = np.array([int(np.sum(cmap[cl] >= 0)) for cl in range(ncls)], np.int64)
    max_act = max(1, int(nact.max()))
    lcols = np.zeros((ncls, max_act), np.int64)
    gcols = np.zeros((ncls, max_act), np.int64)
    for cl in range(ncls):
        idx = np.nonzero(cmap[cl] >= 0)[0]
        lcols[cl, :len(idx)] = idx
        gcols[cl, :len(idx)] = cmap[cl, idx]
    return lcols, gcols, nact


@njit(parallel=True, fastmath=True, cache=True)
def _line_kernel_single(x, T, refl, r_act, t_hat, mid0, rot, tq, W, off, tlo, thi, ltab, g, gw,
                        lam, d, n, k, H, A, lcols, gcols, nact, m_lim, sgn, ncomp):
    """
    Exact 6-scalar polynomial decomposition for an electric + magnetic line pair y(t) = mid0 + t t_hat:
      * along d(t) = p_perp - (t - pt) t_hat, (Tm · d)(rot · d) is quadratic and (rot × Tm) · d is linear,
        so panel quadrature integrates six scalar radial functions in a branch-free SIMD loop;
      * targets are grouped into L2-fitting blocks of 16 and flushed over active columns into F-ordered A.
    """
    M = int(m_lim.max())
    na, nc, ncls = r_act.shape[0], T.shape[1], sgn.shape[0]
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    th0, th1, th2 = t_hat[0], t_hat[1], t_hat[2]
    m0_0, m0_1, m0_2 = mid0[0], mid0[1], mid0[2]
    max_nq = max(g.shape[0], int(off[1] - off[0]))
    max_act = lcols.shape[1]

    ft = np.empty(3)
    for cc in range(3):
        ft[cc] = rot[cc, 0] * th0 + rot[cc, 1] * th1 + rot[cc, 2] * th2

    BM = 16
    nblk = (M + BM - 1) // BM
    for blk in prange(nblk):
        m_start = blk * BM
        m_end = min(M, m_start + BM)
        b_rows = (m_end - m_start) * nc
        B_tile = np.zeros((ncls, max_act, b_rows), np.complex128)
        V_loc = np.zeros((6 * n, nc), np.complex128)

        I_acc = np.empty((12, n))
        S_pan = np.empty((12, max_nq))
        t_pan = np.empty(max_nq)
        W_loc = np.empty((max_nq, n))
        D = np.empty((nc, 3))
        q0_c = np.empty((nc, 3))
        q1_c = np.empty((nc, 3))
        q2_c = np.empty((nc, 3))
        cp_c = np.empty((nc, 3))
        ct_c = np.empty((nc, 3))
        Tm = np.empty((nc, 3))
        up = np.empty(nc)
        ut = np.empty(nc)
        fp = np.empty(3)
        xm = np.empty(3)
        stack_l = np.empty(2 * LHARD + 8, np.int64)
        stack_i = np.empty(2 * LHARD + 8, np.int64)

        for m in range(m_start, m_end):
            loc_r0 = (m - m_start) * nc
            for ia in range(na):
                r = r_act[ia]
                for j in range(3):
                    xm[j] = x[m, j] * refl[r, j]
                    for a in range(nc):
                        Tm[a, j] = T[m * tmul, a, j] * refl[r, j]
                px0, px1, px2 = xm[0] - m0_0, xm[1] - m0_1, xm[2] - m0_2
                pt = px0 * th0 + px1 * th1 + px2 * th2
                p0, p1, p2 = px0 - pt * th0, px1 - pt * th1, px2 - pt * th2
                rho2 = p0 * p0 + p1 * p1 + p2 * p2

                for a in range(nc):
                    up[a] = Tm[a, 0] * p0 + Tm[a, 1] * p1 + Tm[a, 2] * p2
                    ut[a] = Tm[a, 0] * th0 + Tm[a, 1] * th1 + Tm[a, 2] * th2
                for cc in range(ncomp):
                    fp[cc] = rot[cc, 0] * p0 + rot[cc, 1] * p1 + rot[cc, 2] * p2
                for a in range(nc):
                    u0, u1, u2 = Tm[a, 0], Tm[a, 1], Tm[a, 2]
                    upa, uta = up[a], ut[a]
                    for cc in range(ncomp):
                        f0, f1, f2 = rot[cc, 0], rot[cc, 1], rot[cc, 2]
                        D[a, cc] = u0 * f0 + u1 * f1 + u2 * f2
                        q0_c[a, cc] = upa * fp[cc]
                        q1_c[a, cc] = upa * ft[cc] + uta * fp[cc]
                        q2_c[a, cc] = uta * ft[cc]
                        c0, c1, c2 = f1 * u2 - f2 * u1, f2 * u0 - f0 * u2, f0 * u1 - f1 * u0
                        cp_c[a, cc] = c0 * p0 + c1 * p1 + c2 * p2
                        ct_c[a, cc] = c0 * th0 + c1 * th1 + c2 * th2

                I_acc[:] = 0.0
                sp, stack_l[0], stack_i[0] = 1, 0, 0
                while sp > 0:
                    sp -= 1
                    lev, idx = stack_l[sp], stack_i[sp]
                    wid = 2.0 / (1 << lev)
                    p = (1 << lev) - 1 + idx if lev <= ltab else -1
                    a0 = tlo[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * idx))
                    b0 = thi[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * (idx + 1)))
                    gap = max(a0 - pt, 0.0, pt - b0)
                    if b0 - a0 > np.sqrt(rho2 + gap * gap) and lev < LHARD:
                        stack_l[sp], stack_i[sp], stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx, lev + 1, 2 * idx + 1
                        sp += 2
                        continue
                    if p >= 0:
                        q0 = off[p]
                        nq = off[p + 1] - q0
                        for iq in range(nq):
                            t_pan[iq] = tq[q0 + iq]
                    else:
                        q0 = -1
                        nq = g.shape[0]
                        u_mid = -1.0 + wid * (idx + 0.5)
                        u_half = 0.5 * wid
                        w_fac = u_half * H * 2.0 * lam / d
                        for iq in range(nq):
                            u = u_mid + u_half * g[iq]
                            ex = np.exp(lam * u)
                            inv_ex1 = 1.0 / (ex + 1.0)
                            t_pan[iq] = H * (ex - 1.0) * inv_ex1
                            w0 = w_fac * gw[iq] * ex * (inv_ex1 * inv_ex1)
                            W_loc[iq, 0] = w0
                            if n > 1:
                                W_loc[iq, 1] = w0 * u
                            for kk in range(2, n):
                                W_loc[iq, kk] = ((2 * kk - 1) * u * W_loc[iq, kk - 1] - (kk - 1) * W_loc[iq, kk - 2]) / kk

                    for iq in range(nq):
                        tau = t_pan[iq] - pt
                        tau2 = tau * tau
                        r2 = rho2 + tau2
                        inv_r = 1.0 / np.sqrt(r2)
                        inv_r2 = inv_r * inv_r
                        kr = k * (r2 * inv_r)
                        ik = inv_r * inv_k
                        ik2 = ik * ik
                        phr = np.cos(kr) * inv_r
                        phi = np.sin(kr) * inv_r
                        c1 = 1.0 - ik2
                        c3 = 1.0 - 3.0 * ik2
                        ik3 = 3.0 * ik
                        fIr = c1 * phr - ik * phi
                        fIi = c1 * phi + ik * phr
                        fRr = (c3 * phr - ik3 * phi) * inv_r2
                        fRi = (c3 * phi + ik3 * phr) * inv_r2
                        fMr = (-ik * phr - phi) * inv_r
                        fMi = (phr - ik * phi) * inv_r
                        S_pan[0, iq], S_pan[1, iq] = fIr, fIi
                        S_pan[2, iq], S_pan[3, iq] = fRr, fRi
                        S_pan[4, iq], S_pan[5, iq] = tau * fRr, tau * fRi
                        S_pan[6, iq], S_pan[7, iq] = tau2 * fRr, tau2 * fRi
                        S_pan[8, iq], S_pan[9, iq] = fMr, fMi
                        S_pan[10, iq], S_pan[11, iq] = tau * fMr, tau * fMi

                    if q0 >= 0:
                        for iq in range(nq):
                            wq_row = q0 + iq
                            for ch_s in range(12):
                                sv = S_pan[ch_s, iq]
                                for kk in range(n):
                                    I_acc[ch_s, kk] += sv * W[wq_row, kk]
                    else:
                        for iq in range(nq):
                            for ch_s in range(12):
                                sv = S_pan[ch_s, iq]
                                for kk in range(n):
                                    I_acc[ch_s, kk] += sv * W_loc[iq, kk]

                for a in range(nc):
                    for cc in range(ncomp):
                        Dac = D[a, cc]
                        q0v, q1v, q2v = q0_c[a, cc], q1_c[a, cc], q2_c[a, cc]
                        cpv, ctv = cp_c[a, cc], ct_c[a, cc]
                        for kk in range(n):
                            er = Dac * I_acc[0, kk] - q0v * I_acc[2, kk] + q1v * I_acc[4, kk] - q2v * I_acc[6, kk]
                            ei = Dac * I_acc[1, kk] - q0v * I_acc[3, kk] + q1v * I_acc[5, kk] - q2v * I_acc[7, kk]
                            mr = cpv * I_acc[8, kk] - ctv * I_acc[10, kk]
                            mi = cpv * I_acc[9, kk] - ctv * I_acc[11, kk]
                            V_loc[kk * 3 + cc, a] = 2.0 * (er + 1j * ei)
                            V_loc[(n + kk) * 3 + cc, a] = 2.0 * (mr + 1j * mi)

                for cl in range(ncls):
                    if m < m_lim[cl]:
                        sg = sgn[cl, r]
                        na_cl = nact[cl]
                        for j_act in range(na_cl):
                            lc = lcols[cl, j_act]
                            for a in range(nc):
                                B_tile[cl, j_act, loc_r0 + a] += sg * V_loc[lc, a]

        row0 = m_start * nc
        for cl in range(ncls):
            if m_start < m_lim[cl]:
                r_cnt = (min(m_end, m_lim[cl]) - m_start) * nc
                na_cl = nact[cl]
                for j_act in range(na_cl):
                    gc = gcols[cl, j_act]
                    for rr in range(r_cnt):
                        A[cl, row0 + rr, gc] += B_tile[cl, j_act, rr]


@njit(parallel=True, fastmath=True, cache=True)
def _line_kernel_multi(x, T, refl, r_act, t_hat, mids, w, rot, tq, W, off, tlo, thi, ltab, g, gw,
                       lam, d, n, k, H, A, lcols, gcols, nact, m_lim, sgn, ncomp):
    """
    Single-tree evaluation for a bundle of parallel lines y_i(t) = mids[i] + t t_hat
    weighted by w[kind, s, i, cc] (such as the nine-line singular edge tails): sums the
    weighted kernel across all npos lines at each panel node before multiplying by W.
    """
    M = int(m_lim.max())
    na, npos, ns, nc, ncls = r_act.shape[0], mids.shape[0], w.shape[1], T.shape[1], sgn.shape[0]
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    th0, th1, th2 = t_hat[0], t_hat[1], t_hat[2]
    max_nq = max(g.shape[0], int(off[1] - off[0]))
    max_act = lcols.shape[1]
    nch = 4 * ns * nc * ncomp

    ft = np.empty(3)
    for cc in range(3):
        ft[cc] = rot[cc, 0] * th0 + rot[cc, 1] * th1 + rot[cc, 2] * th2

    BM = 16
    nblk = (M + BM - 1) // BM
    for blk in prange(nblk):
        m_start = blk * BM
        m_end = min(M, m_start + BM)
        b_rows = (m_end - m_start) * nc
        B_tile = np.zeros((ncls, max_act, b_rows), np.complex128)
        V_loc = np.zeros((2 * ns * n * 3, nc), np.complex128)

        I_tot = np.empty((nch, n))
        K_pan = np.empty((nch, max_nq))
        t_pan = np.empty(max_nq)
        W_loc = np.empty((max_nq, n))
        Tm = np.empty((nc, 3))
        xm = np.empty(3)
        p_arr = np.empty((npos, 3))
        pt_arr = np.empty(npos)
        rho2_arr = np.empty(npos)
        D = np.empty((nc, 3))
        q0_c = np.empty((npos, nc, 3))
        q1_c = np.empty((npos, nc, 3))
        q2_c = np.empty((nc, 3))
        cp_c = np.empty((npos, nc, 3))
        ct_c = np.empty((nc, 3))
        ut = np.empty(nc)
        stack_l = np.empty(2 * LHARD + 8, np.int64)
        stack_i = np.empty(2 * LHARD + 8, np.int64)

        for m in range(m_start, m_end):
            loc_r0 = (m - m_start) * nc
            for ia in range(na):
                r = r_act[ia]
                for j in range(3):
                    xm[j] = x[m, j] * refl[r, j]
                    for a in range(nc):
                        Tm[a, j] = T[m * tmul, a, j] * refl[r, j]

                min_rho2 = 1e300
                for i in range(npos):
                    px0, px1, px2 = xm[0] - mids[i, 0], xm[1] - mids[i, 1], xm[2] - mids[i, 2]
                    pt = px0 * th0 + px1 * th1 + px2 * th2
                    d0, d1, d2 = px0 - pt * th0, px1 - pt * th1, px2 - pt * th2
                    p_arr[i, 0], p_arr[i, 1], p_arr[i, 2] = d0, d1, d2
                    pt_arr[i] = pt
                    r2v = d0 * d0 + d1 * d1 + d2 * d2
                    rho2_arr[i] = r2v
                    if r2v < min_rho2:
                        min_rho2 = r2v

                for a in range(nc):
                    u0, u1, u2 = Tm[a, 0], Tm[a, 1], Tm[a, 2]
                    ut[a] = u0 * th0 + u1 * th1 + u2 * th2
                    for cc in range(ncomp):
                        f0, f1, f2 = rot[cc, 0], rot[cc, 1], rot[cc, 2]
                        D[a, cc] = u0 * f0 + u1 * f1 + u2 * f2
                        q2_c[a, cc] = ut[a] * ft[cc]
                        c0, c1, c2 = f1 * u2 - f2 * u1, f2 * u0 - f0 * u2, f0 * u1 - f1 * u0
                        ct_c[a, cc] = c0 * th0 + c1 * th1 + c2 * th2
                        for i in range(npos):
                            p0, p1, p2 = p_arr[i, 0], p_arr[i, 1], p_arr[i, 2]
                            upa = u0 * p0 + u1 * p1 + u2 * p2
                            fpc = f0 * p0 + f1 * p1 + f2 * p2
                            q0_c[i, a, cc] = upa * fpc
                            q1_c[i, a, cc] = upa * ft[cc] + ut[a] * fpc
                            cp_c[i, a, cc] = c0 * p0 + c1 * p1 + c2 * p2

                pt0 = pt_arr[0]
                I_tot[:] = 0.0
                sp, stack_l[0], stack_i[0] = 1, 0, 0
                while sp > 0:
                    sp -= 1
                    lev, idx = stack_l[sp], stack_i[sp]
                    wid = 2.0 / (1 << lev)
                    p = (1 << lev) - 1 + idx if lev <= ltab else -1
                    a0 = tlo[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * idx))
                    b0 = thi[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * (idx + 1)))
                    gap = max(a0 - pt0, 0.0, pt0 - b0)
                    if b0 - a0 > np.sqrt(min_rho2 + gap * gap) and lev < LHARD:
                        stack_l[sp], stack_i[sp], stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx, lev + 1, 2 * idx + 1
                        sp += 2
                        continue
                    if p >= 0:
                        q0 = off[p]
                        nq = off[p + 1] - q0
                        for iq in range(nq):
                            t_pan[iq] = tq[q0 + iq]
                    else:
                        q0 = -1
                        nq = g.shape[0]
                        u_mid = -1.0 + wid * (idx + 0.5)
                        u_half = 0.5 * wid
                        w_fac = u_half * H * 2.0 * lam / d
                        for iq in range(nq):
                            u = u_mid + u_half * g[iq]
                            ex = np.exp(lam * u)
                            inv_ex1 = 1.0 / (ex + 1.0)
                            t_pan[iq] = H * (ex - 1.0) * inv_ex1
                            w0 = w_fac * gw[iq] * ex * (inv_ex1 * inv_ex1)
                            W_loc[iq, 0] = w0
                            if n > 1:
                                W_loc[iq, 1] = w0 * u
                            for kk in range(2, n):
                                W_loc[iq, kk] = ((2 * kk - 1) * u * W_loc[iq, kk - 1] - (kk - 1) * W_loc[iq, kk - 2]) / kk

                    K_pan[:, :nq] = 0.0
                    for i in range(npos):
                        pt_i = pt_arr[i]
                        rho2_i = rho2_arr[i]
                        for iq in range(nq):
                            tau = t_pan[iq] - pt_i
                            tau2 = tau * tau
                            r2 = rho2_i + tau2
                            inv_r = 1.0 / np.sqrt(r2)
                            inv_r2 = inv_r * inv_r
                            kr = k * (r2 * inv_r)
                            ik = inv_r * inv_k
                            ik2 = ik * ik
                            phr = np.cos(kr) * inv_r
                            phi = np.sin(kr) * inv_r
                            c1 = 1.0 - ik2
                            c3 = 1.0 - 3.0 * ik2
                            ik3 = 3.0 * ik
                            fIr = c1 * phr - ik * phi
                            fIi = c1 * phi + ik * phr
                            fRr = (c3 * phr - ik3 * phi) * inv_r2
                            fRi = (c3 * phi + ik3 * phr) * inv_r2
                            fMr = (-ik * phr - phi) * inv_r
                            fMi = (phr - ik * phi) * inv_r
                            tfRr, tfRi = tau * fRr, tau * fRi
                            t2fRr, t2fRi = tau2 * fRr, tau2 * fRi
                            tfMr, tfMi = tau * fMr, tau * fMi
                            ch = 0
                            for s in range(ns):
                                for a in range(nc):
                                    for cc in range(ncomp):
                                        Dac = D[a, cc]
                                        q0v, q1v, q2v = q0_c[i, a, cc], q1_c[i, a, cc], q2_c[a, cc]
                                        cpv, ctv = cp_c[i, a, cc], ct_c[a, cc]
                                        er = Dac * fIr - q0v * fRr + q1v * tfRr - q2v * t2fRr
                                        ei = Dac * fIi - q0v * fRi + q1v * tfRi - q2v * t2fRi
                                        mr = cpv * fMr - ctv * tfMr
                                        mi = cpv * fMi - ctv * tfMi
                                        we = 2.0 * w[0, s, i, cc]
                                        wm = 2.0 * w[1, s, i, cc]
                                        K_pan[ch, iq] += we * er
                                        K_pan[ch + 1, iq] += we * ei
                                        K_pan[ch + 2, iq] += wm * mr
                                        K_pan[ch + 3, iq] += wm * mi
                                        ch += 4

                    if q0 >= 0:
                        for iq in range(nq):
                            wq_row = q0 + iq
                            for ch in range(nch):
                                sv = K_pan[ch, iq]
                                for kk in range(n):
                                    I_tot[ch, kk] += sv * W[wq_row, kk]
                    else:
                        for iq in range(nq):
                            for ch in range(nch):
                                sv = K_pan[ch, iq]
                                for kk in range(n):
                                    I_tot[ch, kk] += sv * W_loc[iq, kk]

                ch = 0
                for s in range(ns):
                    for a in range(nc):
                        for cc in range(ncomp):
                            for kk in range(n):
                                col_e = ((0 * ns + s) * n + kk) * 3 + cc
                                col_m = ((1 * ns + s) * n + kk) * 3 + cc
                                V_loc[col_e, a] = I_tot[ch, kk] + 1j * I_tot[ch + 1, kk]
                                V_loc[col_m, a] = I_tot[ch + 2, kk] + 1j * I_tot[ch + 3, kk]
                            ch += 4

                for cl in range(ncls):
                    if m < m_lim[cl]:
                        sg = sgn[cl, r]
                        na_cl = nact[cl]
                        for j_act in range(na_cl):
                            lc = lcols[cl, j_act]
                            for a in range(nc):
                                B_tile[cl, j_act, loc_r0 + a] += sg * V_loc[lc, a]

        row0 = m_start * nc
        for cl in range(ncls):
            if m_start < m_lim[cl]:
                r_cnt = (min(m_end, m_lim[cl]) - m_start) * nc
                na_cl = nact[cl]
                for j_act in range(na_cl):
                    gc = gcols[cl, j_act]
                    for rr in range(r_cnt):
                        A[cl, row0 + rr, gc] += B_tile[cl, j_act, rr]


def line_add(x, T, refl, r_act, t_hat, mids, lam, d, n, k, H, A, cmap, sgn, rot=None, w=None, m_lim=None):
    """General wrapper for parallel skeleton lines y_i(t) = mids[i] + t t_hat in frame rot"""
    mids = np.ascontiguousarray(np.atleast_2d(mids), dtype=float)
    w_arr = np.ones((2, 1, len(mids), 3)) if w is None else np.ascontiguousarray(w)
    cmap = np.ascontiguousarray(cmap)
    lcols, gcols, nact = _pack_cmap(cmap)
    _, _, ncomp = _mode_strides(cmap, w_arr.shape[1], n)
    tq, W, off, tlo, thi, ltab, g, gw = build_tree(lam, d, n, H, k)
    rot_arr = np.eye(3) if rot is None else np.ascontiguousarray(rot)
    x_c, T_c, r_c = np.ascontiguousarray(x), np.ascontiguousarray(T), np.ascontiguousarray(r_act, np.int64)
    th_c = np.ascontiguousarray(t_hat, dtype=float)
    ml = np.full(sgn.shape[0], len(x_c), np.int64) if m_lim is None else np.ascontiguousarray(m_lim, np.int64)
    if len(mids) == 1 and w_arr.shape[1] == 1 and (w is None or np.all(w_arr == 1.0)):
        _line_kernel_single(x_c, T_c, refl, r_c, th_c, mids[0], rot_arr, tq, W, off, tlo, thi,
                            ltab, g, gw, lam, d, n, k, H, A, lcols, gcols, nact, ml, sgn, ncomp)
    else:
        _line_kernel_multi(x_c, T_c, refl, r_c, th_c, mids, w_arr, rot_arr, tq, W, off, tlo, thi,
                           ltab, g, gw, lam, d, n, k, H, A, lcols, gcols, nact, ml, sgn, ncomp)


@njit(parallel=True, fastmath=True, cache=True)
def _points_kernel(x, T, refl, pe, pm, n_only_e, k, A, col0, m_lim, sgn):
    """Electric dipoles at pe and magnetic at pm, assembled columnwise into F-ordered A"""
    M = int(m_lim.max())
    nc, ne = T.shape[1], pe.shape[0]
    ncls, nr = sgn.shape
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    for j in prange(ne):
        has_m = j >= n_only_e
        y0, y1, y2 = pe[j, 0], pe[j, 1], pe[j, 2]
        acce = np.empty((ncls, nc, 3), np.complex128)
        accm = np.empty((ncls, nc, 3), np.complex128)
        for m in range(M):
            acce[:] = 0.0
            if has_m:
                accm[:] = 0.0
            for r in range(nr):
                e0 = x[m, 0] * refl[r, 0] - y0
                e1 = x[m, 1] * refl[r, 1] - y1
                e2 = x[m, 2] * refl[r, 2] - y2
                r2 = e0 * e0 + e1 * e1 + e2 * e2
                inv_r = 1.0 / np.sqrt(r2)
                inv_r2 = inv_r * inv_r
                kr = k * (r2 * inv_r)
                ik = inv_r * inv_k
                ik2 = ik * ik
                ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
                fI = (1.0 - ik2 + 1j * ik) * ph
                fR = (1.0 - 3.0 * ik2 + 3j * ik) * ph * inv_r2
                if not has_m:
                    for a in range(nc):
                        t0 = T[m * tmul, a, 0] * refl[r, 0]
                        t1 = T[m * tmul, a, 1] * refl[r, 1]
                        t2 = T[m * tmul, a, 2] * refl[r, 2]
                        td = fR * (t0 * e0 + t1 * e1 + t2 * e2)
                        ge0 = fI * t0 - td * e0
                        ge1 = fI * t1 - td * e1
                        ge2 = fI * t2 - td * e2
                        for s in range(ncls):
                            sg = sgn[s, r]
                            acce[s, a, 0] += sg * ge0
                            acce[s, a, 1] += sg * ge1
                            acce[s, a, 2] += sg * ge2
                else:
                    fM = (-ik + 1j) * ph * inv_r
                    for a in range(nc):
                        t0 = T[m * tmul, a, 0] * refl[r, 0]
                        t1 = T[m * tmul, a, 1] * refl[r, 1]
                        t2 = T[m * tmul, a, 2] * refl[r, 2]
                        td = fR * (t0 * e0 + t1 * e1 + t2 * e2)
                        ge0 = fI * t0 - td * e0
                        ge1 = fI * t1 - td * e1
                        ge2 = fI * t2 - td * e2
                        gm0 = fM * (t1 * e2 - t2 * e1)
                        gm1 = fM * (t2 * e0 - t0 * e2)
                        gm2 = fM * (t0 * e1 - t1 * e0)
                        for s in range(ncls):
                            sg = sgn[s, r]
                            acce[s, a, 0] += sg * ge0
                            acce[s, a, 1] += sg * ge1
                            acce[s, a, 2] += sg * ge2
                            accm[s, a, 0] += sg * gm0
                            accm[s, a, 1] += sg * gm1
                            accm[s, a, 2] += sg * gm2
            for s in range(ncls):
                if m < m_lim[s]:
                    ce = col0[s] + 3 * j
                    for a in range(nc):
                        A[s, m * nc + a, ce] += acce[s, a, 0]
                        A[s, m * nc + a, ce + 1] += acce[s, a, 1]
                        A[s, m * nc + a, ce + 2] += acce[s, a, 2]
                    if has_m:
                        cm = col0[s] + 3 * (ne + (j - n_only_e))
                        for a in range(nc):
                            A[s, m * nc + a, cm] += accm[s, a, 0]
                            A[s, m * nc + a, cm + 1] += accm[s, a, 1]
                            A[s, m * nc + a, cm + 2] += accm[s, a, 2]


def points_add(x, T, refl, pe, pm, k, A, col0, sgn, m_lim=None):
    """As _line_kernel, for electric point dipoles at pe and magnetic ones at pm (columns col0[s] + 3 j + c in class s)"""
    ne, nm = len(pe), len(pm)
    ml = np.full(sgn.shape[0], len(x), np.int64) if m_lim is None else np.ascontiguousarray(m_lim, np.int64)
    if nm > 0 and ne >= nm and np.array_equal(pe[ne - nm:], pm):
        _points_kernel(x, T, refl, pe, pm, ne - nm, k, A, col0, ml, sgn)
    else:
        _points_kernel(x, T, refl, pe, np.empty((0, 3)), ne, k, A, col0, ml, sgn)
        if nm > 0:
            _points_m_only(x, T, refl, pm, k, A, col0 + 3 * ne, ml, sgn)


@njit(parallel=True, fastmath=True, cache=True)
def _points_m_only(x, T, refl, pm, k, A, col0, m_lim, sgn):
    M = int(m_lim.max())
    nc, nm = T.shape[1], pm.shape[0]
    ncls, nr = sgn.shape
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    for j in prange(nm):
        y0, y1, y2 = pm[j, 0], pm[j, 1], pm[j, 2]
        accm = np.empty((ncls, nc, 3), np.complex128)
        for m in range(M):
            accm[:] = 0.0
            for r in range(nr):
                e0 = x[m, 0] * refl[r, 0] - y0
                e1 = x[m, 1] * refl[r, 1] - y1
                e2 = x[m, 2] * refl[r, 2] - y2
                r2 = e0 * e0 + e1 * e1 + e2 * e2
                inv_r = 1.0 / np.sqrt(r2)
                kr = k * (r2 * inv_r)
                ik = inv_r * inv_k
                ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
                fM = (-ik + 1j) * ph * inv_r
                for a in range(nc):
                    t0 = T[m * tmul, a, 0] * refl[r, 0]
                    t1 = T[m * tmul, a, 1] * refl[r, 1]
                    t2 = T[m * tmul, a, 2] * refl[r, 2]
                    gm0 = fM * (t1 * e2 - t2 * e1)
                    gm1 = fM * (t2 * e0 - t0 * e2)
                    gm2 = fM * (t0 * e1 - t1 * e0)
                    for s in range(ncls):
                        sg = sgn[s, r]
                        accm[s, a, 0] += sg * gm0
                        accm[s, a, 1] += sg * gm1
                        accm[s, a, 2] += sg * gm2
            for s in range(ncls):
                if m < m_lim[s]:
                    cm = col0[s] + 3 * j
                    for a in range(nc):
                        A[s, m * nc + a, cm] += accm[s, a, 0]
                        A[s, m * nc + a, cm + 1] += accm[s, a, 1]
                        A[s, m * nc + a, cm + 2] += accm[s, a, 2]


@njit(parallel=True, fastmath=True, cache=True)
def _radial_tails_kernel(x, T, refl, ray_pts, we, wm, k, A, col_e, col_m, m_lim, sgn):
    """Fields of radial corner/vertex tail rays (sharing electric and magnetic kernel evaluations)"""
    M = int(m_lim.max())
    nc = T.shape[1]
    ncls, nr = sgn.shape
    R_cnt, nstep = ray_pts.shape[0], ray_pts.shape[1]
    ne, nm = we.shape[0], wm.shape[0]
    tmul = 1 if T.shape[0] > 1 else 0
    inv_k = 1.0 / k
    for m in prange(M):
        acce = np.empty((ncls, ne, nc, 3), np.complex128)
        accm = np.empty((ncls, nm, nc, 3), np.complex128)
        Ge = np.empty((nc, 3), np.complex128)
        Gm = np.empty((nc, 3), np.complex128)
        Tm = np.empty((nc, 3))
        for ir in range(R_cnt):
            acce[:] = 0.0
            accm[:] = 0.0
            for r in range(nr):
                for a in range(nc):
                    for q in range(3):
                        Tm[a, q] = T[m * tmul, a, q] * refl[r, q]
                for t in range(nstep):
                    e0 = x[m, 0] * refl[r, 0] - ray_pts[ir, t, 0]
                    e1 = x[m, 1] * refl[r, 1] - ray_pts[ir, t, 1]
                    e2 = x[m, 2] * refl[r, 2] - ray_pts[ir, t, 2]
                    r2 = e0 * e0 + e1 * e1 + e2 * e2
                    inv_r = 1.0 / np.sqrt(r2)
                    inv_r2 = inv_r * inv_r
                    kr = k * (r2 * inv_r)
                    ik = inv_r * inv_k
                    ik2 = ik * ik
                    ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
                    fI = (1.0 - ik2 + 1j * ik) * ph
                    fR = (1.0 - 3.0 * ik2 + 3j * ik) * ph * inv_r2
                    fM = (-ik + 1j) * ph * inv_r
                    for a in range(nc):
                        t0, t1, t2 = Tm[a, 0], Tm[a, 1], Tm[a, 2]
                        td = fR * (t0 * e0 + t1 * e1 + t2 * e2)
                        Ge[a, 0] = fI * t0 - td * e0
                        Ge[a, 1] = fI * t1 - td * e1
                        Ge[a, 2] = fI * t2 - td * e2
                        Gm[a, 0] = fM * (t1 * e2 - t2 * e1)
                        Gm[a, 1] = fM * (t2 * e0 - t0 * e2)
                        Gm[a, 2] = fM * (t0 * e1 - t1 * e0)
                    for s in range(ncls):
                        sg = sgn[s, r]
                        for ie in range(ne):
                            wv = sg * we[ie, t]
                            for a in range(nc):
                                acce[s, ie, a, 0] += wv * Ge[a, 0]
                                acce[s, ie, a, 1] += wv * Ge[a, 1]
                                acce[s, ie, a, 2] += wv * Ge[a, 2]
                        for im in range(nm):
                            wv = sg * wm[im, t]
                            for a in range(nc):
                                accm[s, im, a, 0] += wv * Gm[a, 0]
                                accm[s, im, a, 1] += wv * Gm[a, 1]
                                accm[s, im, a, 2] += wv * Gm[a, 2]
            for s in range(ncls):
                if m < m_lim[s]:
                    be = col_e[s] + ir * ne * 3
                    for ie in range(ne):
                        for a in range(nc):
                            A[s, m * nc + a, be + ie * 3] += acce[s, ie, a, 0]
                            A[s, m * nc + a, be + ie * 3 + 1] += acce[s, ie, a, 1]
                            A[s, m * nc + a, be + ie * 3 + 2] += acce[s, ie, a, 2]
                    bm = col_m[s] + ir * nm * 3
                    for im in range(nm):
                        for a in range(nc):
                            A[s, m * nc + a, bm + im * 3] += accm[s, im, a, 0]
                            A[s, m * nc + a, bm + im * 3 + 1] += accm[s, im, a, 1]
                            A[s, m * nc + a, bm + im * 3 + 2] += accm[s, im, a, 2]


def radial_tails_add(x, T, refl, ray_pts, we, wm, k, A, col_e, col_m, sgn, m_lim=None):
    ml = np.full(sgn.shape[0], len(x), np.int64) if m_lim is None else np.ascontiguousarray(m_lim, np.int64)
    _radial_tails_kernel(x, T, refl, ray_pts, we, wm, k, A, col_e, col_m, ml, sgn)


# Direct matrix-free field and exact r -> infinity far-field evaluation from solved mode coefficients

@njit(parallel=True, fastmath=True, cache=True)
def eval_points_field(pts, refl, sgn, y, alpha, beta, k, out):
    """Adds Σ_{cl, r, j} sgn[cl, r] R_r E_dipole(R_r pts[m]; y[j], alpha[cl, j], beta[cl, j]) into out[cl, m, 3]"""
    P, ncls, nr, N = pts.shape[0], sgn.shape[0], refl.shape[0], y.shape[0]
    inv_k = 1.0 / k
    for m in prange(P):
        for r in range(nr):
            rx0, rx1, rx2 = refl[r, 0], refl[r, 1], refl[r, 2]
            xm0, xm1, xm2 = pts[m, 0] * rx0, pts[m, 1] * rx1, pts[m, 2] * rx2
            for j in range(N):
                d0 = xm0 - y[j, 0]
                d1 = xm1 - y[j, 1]
                d2 = xm2 - y[j, 2]
                r2 = d0 * d0 + d1 * d1 + d2 * d2
                inv_r = 1.0 / np.sqrt(r2)
                inv_r2 = inv_r * inv_r
                kr = k * (r2 * inv_r)
                ik = inv_r * inv_k
                ik2 = ik * ik
                ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
                fI = (1.0 - ik2 + 1j * ik) * ph
                fR = (1.0 - 3.0 * ik2 + 3j * ik) * ph * inv_r2
                fM = (-ik + 1j) * ph * inv_r
                for cl in range(ncls):
                    sg = sgn[cl, r]
                    a0, a1, a2 = alpha[cl, j, 0], alpha[cl, j, 1], alpha[cl, j, 2]
                    b0, b1, b2 = beta[cl, j, 0], beta[cl, j, 1], beta[cl, j, 2]
                    da = fR * (d0 * a0 + d1 * a1 + d2 * a2)
                    e0 = fI * a0 - da * d0 + fM * (d1 * b2 - d2 * b1)
                    e1 = fI * a1 - da * d1 + fM * (d2 * b0 - d0 * b2)
                    e2 = fI * a2 - da * d2 + fM * (d0 * b1 - d1 * b0)
                    out[cl, m, 0] += sg * rx0 * e0
                    out[cl, m, 1] += sg * rx1 * e1
                    out[cl, m, 2] += sg * rx2 * e2


@njit(parallel=True, fastmath=True, cache=True)
def _eval_line_field_kernel(pts, refl, r_act, sgn, t_hat, mids, tq, W, off, tlo, thi, ltab, g, gw,
                            lam, d, n, k, H, n_mom, M_mom, c_cart, out):
    """Adds field of a solved line/tail pair (pre-contracted on the 15 moment panels) into out[cl, m, 3]"""
    P, na, npos, ncls = pts.shape[0], r_act.shape[0], mids.shape[0], sgn.shape[0]
    inv_k = 1.0 / k
    th0, th1, th2 = t_hat[0], t_hat[1], t_hat[2]
    nch = min(64, max(1, P))
    for ch in prange(nch):
        px = np.empty((npos, 3))
        Wq = np.empty(n)
        mom = np.empty((ncls, 2, npos, 3), np.complex128)
        stack_l = np.empty(2 * LHARD + 8, np.int64)
        stack_i = np.empty(2 * LHARD + 8, np.int64)
        for m in range(ch, P, nch):
            for ia in range(na):
                r = r_act[ia]
                rx0, rx1, rx2 = refl[r, 0], refl[r, 1], refl[r, 2]
                xm0, xm1, xm2 = pts[m, 0] * rx0, pts[m, 1] * rx1, pts[m, 2] * rx2
                pt0 = (xm0 - mids[0, 0]) * th0 + (xm1 - mids[0, 1]) * th1 + (xm2 - mids[0, 2]) * th2
                rho2 = 1e300
                for i in range(npos):
                    px0 = xm0 - mids[i, 0]
                    px1 = xm1 - mids[i, 1]
                    px2 = xm2 - mids[i, 2]
                    px[i, 0], px[i, 1], px[i, 2] = px0, px1, px2
                    pt_i = px0 * th0 + px1 * th1 + px2 * th2
                    p0 = px0 - pt_i * th0
                    p1 = px1 - pt_i * th1
                    p2 = px2 - pt_i * th2
                    d2 = p0 * p0 + p1 * p1 + p2 * p2
                    if d2 < rho2:
                        rho2 = d2

                sp, stack_l[0], stack_i[0] = 1, 0, 0
                while sp > 0:
                    sp -= 1
                    lev, idx = stack_l[sp], stack_i[sp]
                    wid = 2.0 / (1 << lev)
                    p = (1 << lev) - 1 + idx if lev <= ltab else -1
                    a0 = tlo[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * idx))
                    b0 = thi[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * (idx + 1)))
                    gap = max(a0 - pt0, 0.0, pt0 - b0)
                    if b0 - a0 > np.sqrt(rho2 + gap * gap) and lev < LHARD:
                        stack_l[sp], stack_i[sp], stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx, lev + 1, 2 * idx + 1
                        sp += 2
                        continue
                    nq = off[p + 1] - off[p] if p >= 0 else g.shape[0]
                    for iq in range(nq):
                        q_idx = off[p] + iq if p >= 0 else -1
                        if 0 <= q_idx < n_mom:
                            t = tq[q_idx]
                        elif q_idx >= 0:
                            t = tq[q_idx]
                            for cl in range(ncls):
                                for kind in range(2):
                                    for i in range(npos):
                                        s0, s1, s2 = 0.0j, 0.0j, 0.0j
                                        for kk in range(n):
                                            wk = W[q_idx, kk]
                                            s0 += wk * c_cart[cl, kind, i, kk, 0]
                                            s1 += wk * c_cart[cl, kind, i, kk, 1]
                                            s2 += wk * c_cart[cl, kind, i, kk, 2]
                                        mom[cl, kind, i, 0] = s0
                                        mom[cl, kind, i, 1] = s1
                                        mom[cl, kind, i, 2] = s2
                        else:
                            u = -1.0 + wid * (idx + 0.5) + 0.5 * wid * g[iq]
                            ex = np.exp(lam * u)
                            t = H * (ex - 1.0) / (ex + 1.0)
                            Wq[0] = 0.5 * wid * gw[iq] * H * 2.0 * lam * ex / ((ex + 1.0) * (ex + 1.0)) / d
                            if n > 1:
                                Wq[1] = Wq[0] * u
                            for kk in range(2, n):
                                Wq[kk] = ((2 * kk - 1) * u * Wq[kk - 1] - (kk - 1) * Wq[kk - 2]) / kk
                            for cl in range(ncls):
                                for kind in range(2):
                                    for i in range(npos):
                                        s0, s1, s2 = 0.0j, 0.0j, 0.0j
                                        for kk in range(n):
                                            wk = Wq[kk]
                                            s0 += wk * c_cart[cl, kind, i, kk, 0]
                                            s1 += wk * c_cart[cl, kind, i, kk, 1]
                                            s2 += wk * c_cart[cl, kind, i, kk, 2]
                                        mom[cl, kind, i, 0] = s0
                                        mom[cl, kind, i, 1] = s1
                                        mom[cl, kind, i, 2] = s2

                        for i in range(npos):
                            d0 = px[i, 0] - t * th0
                            d1 = px[i, 1] - t * th1
                            d2 = px[i, 2] - t * th2
                            r2 = d0 * d0 + d1 * d1 + d2 * d2
                            inv_r = 1.0 / np.sqrt(r2)
                            inv_r2 = inv_r * inv_r
                            kr = k * (r2 * inv_r)
                            ik = inv_r * inv_k
                            ik2 = ik * ik
                            ph = (np.cos(kr) + 1j * np.sin(kr)) * inv_r
                            fI = (1.0 - ik2 + 1j * ik) * ph
                            fR = (1.0 - 3.0 * ik2 + 3j * ik) * ph * inv_r2
                            fM = (-ik + 1j) * ph * inv_r
                            for cl in range(ncls):
                                sg = sgn[cl, r]
                                if 0 <= q_idx < n_mom:
                                    a0 = M_mom[q_idx, cl, 0, i, 0]
                                    a1 = M_mom[q_idx, cl, 0, i, 1]
                                    a2 = M_mom[q_idx, cl, 0, i, 2]
                                    b0 = M_mom[q_idx, cl, 1, i, 0]
                                    b1 = M_mom[q_idx, cl, 1, i, 1]
                                    b2 = M_mom[q_idx, cl, 1, i, 2]
                                else:
                                    a0 = mom[cl, 0, i, 0]
                                    a1 = mom[cl, 0, i, 1]
                                    a2 = mom[cl, 0, i, 2]
                                    b0 = mom[cl, 1, i, 0]
                                    b1 = mom[cl, 1, i, 1]
                                    b2 = mom[cl, 1, i, 2]
                                da = fR * (d0 * a0 + d1 * a1 + d2 * a2)
                                e0 = fI * a0 - da * d0 + fM * (d1 * b2 - d2 * b1)
                                e1 = fI * a1 - da * d1 + fM * (d2 * b0 - d0 * b2)
                                e2 = fI * a2 - da * d2 + fM * (d0 * b1 - d1 * b0)
                                out[cl, m, 0] += sg * rx0 * e0
                                out[cl, m, 1] += sg * rx1 * e1
                                out[cl, m, 2] += sg * rx2 * e2


def eval_line_field(pts, refl, r_act, sgn, t_hat, mids, lam, d, n, k, H, c_pair, out, rot=None, w=None):
    """Pre-contracts mode coefficients on moment panels and adds line/tail fields into out[cl, m, 3]"""
    mids = np.ascontiguousarray(np.atleast_2d(mids), dtype=float)
    w = np.ones((2, 1, len(mids), 3)) if w is None else np.ascontiguousarray(w)
    rot_arr = np.eye(3) if rot is None else np.ascontiguousarray(rot)
    tq, W, off, tlo, thi, ltab, g, gw = build_tree(lam, d, n, H, k)
    c_cart = np.ascontiguousarray(2.0 * np.einsum("ksic,lksmc,cq->lkimq", w, c_pair, rot_arr))
    n_mom = int(off[min(511, len(off) - 1)])
    M_mom = np.ascontiguousarray(np.einsum("pm,lkimq->plkiq", W[:n_mom], c_cart))
    _eval_line_field_kernel(pts, refl, np.ascontiguousarray(r_act, np.int64), sgn,
                            np.ascontiguousarray(t_hat, dtype=float), mids, tq, W, off, tlo, thi,
                            ltab, g, gw, lam, d, n, k, H, n_mom, M_mom, c_cart, out)


@njit(parallel=True, fastmath=True, cache=True)
def eval_far_dipoles(dirs, refl, r_act, sgn, y, alpha, beta, k, out):
    """Exact r -> infinity far-field pattern F(dirs) added into out[cl, D, 3] over active images r_act"""
    D_cnt, na, N, ncls = dirs.shape[0], r_act.shape[0], y.shape[0], sgn.shape[0]
    for m in prange(D_cnt):
        for ia in range(na):
            r = r_act[ia]
            rx0, rx1, rx2 = refl[r, 0], refl[r, 1], refl[r, 2]
            xh0, xh1, xh2 = dirs[m, 0] * rx0, dirs[m, 1] * rx1, dirs[m, 2] * rx2
            for j in range(N):
                kr = -k * (xh0 * y[j, 0] + xh1 * y[j, 1] + xh2 * y[j, 2])
                ph = np.cos(kr) + 1j * np.sin(kr)
                mph = 1j * ph
                for cl in range(ncls):
                    sg = sgn[cl, r]
                    a0, a1, a2 = alpha[cl, j, 0], alpha[cl, j, 1], alpha[cl, j, 2]
                    b0, b1, b2 = beta[cl, j, 0], beta[cl, j, 1], beta[cl, j, 2]
                    xa = xh0 * a0 + xh1 * a1 + xh2 * a2
                    f0 = ph * (a0 - xa * xh0) + mph * (xh1 * b2 - xh2 * b1)
                    f1 = ph * (a1 - xa * xh1) + mph * (xh2 * b0 - xh0 * b2)
                    f2 = ph * (a2 - xa * xh2) + mph * (xh0 * b1 - xh1 * b0)
                    out[cl, m, 0] += sg * rx0 * f0
                    out[cl, m, 1] += sg * rx1 * f1
                    out[cl, m, 2] += sg * rx2 * f2


@njit(parallel=True, cache=True)
def build_c3_tri(A_tri, ic, tri_idx, hi_pts, hi_tri, lc_yz, sg_yz, c0, c1, c2, w1, w2, out):
    """Builds a C_3 Fourier block on Face 0 directly from the 1/48-cube triangle matrix A_tri"""
    nt, nh, nc = len(tri_idx), len(hi_pts), len(c0)
    for k in prange(nc):
        j0, j1, j2 = c0[k], c1[k], c2[k]
        y0, y1, y2 = lc_yz[j0], lc_yz[j1], lc_yz[j2]
        s0, s1, s2 = sg_yz[j0], sg_yz[j1], sg_yz[j2]
        for t in range(nt):
            m = tri_idx[t]
            r0, r1 = 2 * t, 2 * t + 1
            out[2 * m, k] = A_tri[ic, r0, j0] + w1 * A_tri[ic, r0, j1] + w2 * A_tri[ic, r0, j2]
            out[2 * m + 1, k] = A_tri[ic, r1, j0] + w1 * A_tri[ic, r1, j1] + w2 * A_tri[ic, r1, j2]
        for h in range(nh):
            m = hi_pts[h]
            r0, r1 = 2 * hi_tri[h] + 1, 2 * hi_tri[h]
            out[2 * m, k] = s0 * A_tri[ic, r0, y0] + w1 * s1 * A_tri[ic, r0, y1] + w2 * s2 * A_tri[ic, r0, y2]
            out[2 * m + 1, k] = s0 * A_tri[ic, r1, y0] + w1 * s1 * A_tri[ic, r1, y1] + w2 * s2 * A_tri[ic, r1, y2]


@njit(parallel=True, cache=True)
def build_z2_tri(A_tri, ic_a, ic_b, ic_byz, m_row_tri, m_row_ab, lc_b, lc_byz, sg_byz,
                 r_pair, r_fix, c_pair, col_p, col_s, c_fix, sign_eig, out):
    """Builds a Z_2 (y <-> z) symmetric/antisymmetric block directly from the 1/48-cube triangle matrix A_tri"""
    np_r, nf_r, np_c, nf_c = len(r_pair), len(r_fix), len(c_pair), len(c_fix)
    rt2, s2 = np.sqrt(2.0), np.sqrt(0.5)
    for k in prange(np_c):
        c0 = c_pair[k]
        c1 = col_p[c0]
        sc = sign_eig * col_s[c0]
        b0, b1 = lc_b[c0], lc_b[c1]
        y0, y1 = lc_byz[c0], lc_byz[c1]
        sy0, sy1 = sg_byz[c0], sg_byz[c1]
        for i in range(np_r):
            r = r_pair[i]
            rt = m_row_tri[r]
            kd = m_row_ab[r]
            if kd == 0:
                out[i, k] = A_tri[ic_a, rt, c0] + sc * A_tri[ic_a, rt, c1]
            elif kd == 1:
                out[i, k] = A_tri[ic_b, rt, b0] + sc * A_tri[ic_b, rt, b1]
            else:
                out[i, k] = sy0 * A_tri[ic_byz, rt, y0] + sc * sy1 * A_tri[ic_byz, rt, y1]
        for i in range(nf_r):
            rt = m_row_tri[r_fix[i]]
            out[np_r + i, k] = s2 * (A_tri[ic_a, rt, c0] + sc * A_tri[ic_a, rt, c1])
    for k in prange(nf_c):
        c0 = c_fix[k]
        b0, y0, sy0 = lc_b[c0], lc_byz[c0], sg_byz[c0]
        col_out = np_c + k
        for i in range(np_r):
            r = r_pair[i]
            rt = m_row_tri[r]
            kd = m_row_ab[r]
            if kd == 0:
                v = A_tri[ic_a, rt, c0]
            elif kd == 1:
                v = A_tri[ic_b, rt, b0]
            else:
                v = sy0 * A_tri[ic_byz, rt, y0]
            out[i, col_out] = rt2 * v
        for i in range(nf_r):
            out[np_r + i, col_out] = A_tri[ic_a, m_row_tri[r_fix[i]], c0]


@njit(parallel=True, fastmath=True, cache=True)
def eval_points_adj(pts, refl, sgn, u, y, k, out_a, out_b):
    """Exact conjugate transpose of eval_points_field for boundary dipoles u[cl, m, 3] = Σ_a Tw[m, a] conj(v[cl, m, a])"""
    P, ncls, nr, N = pts.shape[0], sgn.shape[0], refl.shape[0], y.shape[0]
    inv_k = 1.0 / k
    for j in prange(N):
        y0, y1, y2 = y[j, 0], y[j, 1], y[j, 2]
        for m in range(P):
            for r in range(nr):
                rx0, rx1, rx2 = refl[r, 0], refl[r, 1], refl[r, 2]
                d0 = pts[m, 0] * rx0 - y0
                d1 = pts[m, 1] * rx1 - y1
                d2 = pts[m, 2] * rx2 - y2
                r2 = d0 * d0 + d1 * d1 + d2 * d2
                inv_r = 1.0 / np.sqrt(r2)
                inv_r2 = inv_r * inv_r
                kr = k * (r2 * inv_r)
                ik = inv_r * inv_k
                ik2 = ik * ik
                phc = (np.cos(kr) - 1j * np.sin(kr)) * inv_r
                fIc = (1.0 - ik2 - 1j * ik) * phc
                fRc = (1.0 - 3.0 * ik2 - 3j * ik) * phc * inv_r2
                fMc = (-ik - 1j) * phc * inv_r
                for cl in range(ncls):
                    sg = sgn[cl, r]
                    u0 = sg * rx0 * np.conj(u[cl, m, 0])
                    u1 = sg * rx1 * np.conj(u[cl, m, 1])
                    u2 = sg * rx2 * np.conj(u[cl, m, 2])
                    du = fRc * (d0 * u0 + d1 * u1 + d2 * u2)
                    out_a[cl, j, 0] += fIc * u0 - du * d0
                    out_a[cl, j, 1] += fIc * u1 - du * d1
                    out_a[cl, j, 2] += fIc * u2 - du * d2
                    out_b[cl, j, 0] += fMc * (u1 * d2 - u2 * d1)
                    out_b[cl, j, 1] += fMc * (u2 * d0 - u0 * d2)
                    out_b[cl, j, 2] += fMc * (u0 * d1 - u1 * d0)


@njit(parallel=True, fastmath=True, cache=True)
def _eval_line_adj_kernel(pts, refl, r_act, sgn, t_hat, mids, tq, W, off, tlo, thi, ltab, g, gw,
                          lam, d, n, k, H, n_mom, u, G_mom, Z_cart):
    """Exact conjugate transpose of _eval_line_field_kernel using moment-panel accumulation on levels < LMOM"""
    P, na, npos, ncls = pts.shape[0], r_act.shape[0], mids.shape[0], sgn.shape[0]
    inv_k = 1.0 / k
    th0, th1, th2 = t_hat[0], t_hat[1], t_hat[2]
    nch = G_mom.shape[0]
    for ch in prange(nch):
        px = np.empty((npos, 3))
        Wq = np.empty(n)
        stack_l = np.empty(2 * LHARD + 8, np.int64)
        stack_i = np.empty(2 * LHARD + 8, np.int64)
        for m in range(ch, P, nch):
            for ia in range(na):
                r = r_act[ia]
                rx0, rx1, rx2 = refl[r, 0], refl[r, 1], refl[r, 2]
                xm0, xm1, xm2 = pts[m, 0] * rx0, pts[m, 1] * rx1, pts[m, 2] * rx2
                pt0 = (xm0 - mids[0, 0]) * th0 + (xm1 - mids[0, 1]) * th1 + (xm2 - mids[0, 2]) * th2
                rho2 = 1e300
                for i in range(npos):
                    px0, px1, px2 = xm0 - mids[i, 0], xm1 - mids[i, 1], xm2 - mids[i, 2]
                    px[i, 0], px[i, 1], px[i, 2] = px0, px1, px2
                    pt_i = px0 * th0 + px1 * th1 + px2 * th2
                    p0, p1, p2 = px0 - pt_i * th0, px1 - pt_i * th1, px2 - pt_i * th2
                    d2 = p0 * p0 + p1 * p1 + p2 * p2
                    if d2 < rho2:
                        rho2 = d2

                sp, stack_l[0], stack_i[0] = 1, 0, 0
                while sp > 0:
                    sp -= 1
                    lev, idx = stack_l[sp], stack_i[sp]
                    wid = 2.0 / (1 << lev)
                    p = (1 << lev) - 1 + idx if lev <= ltab else -1
                    a0 = tlo[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * idx))
                    b0 = thi[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * (idx + 1)))
                    gap = max(a0 - pt0, 0.0, pt0 - b0)
                    if b0 - a0 > np.sqrt(rho2 + gap * gap) and lev < LHARD:
                        stack_l[sp], stack_i[sp], stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx, lev + 1, 2 * idx + 1
                        sp += 2
                        continue
                    nq = off[p + 1] - off[p] if p >= 0 else g.shape[0]
                    for iq in range(nq):
                        q_idx = off[p] + iq if p >= 0 else -1
                        if q_idx >= 0:
                            t = tq[q_idx]
                        else:
                            u_loc = -1.0 + wid * (idx + 0.5) + 0.5 * wid * g[iq]
                            ex = np.exp(lam * u_loc)
                            t = H * (ex - 1.0) / (ex + 1.0)
                            Wq[0] = 0.5 * wid * gw[iq] * H * 2.0 * lam * ex / ((ex + 1.0) * (ex + 1.0)) / d
                            if n > 1:
                                Wq[1] = Wq[0] * u_loc
                            for kk in range(2, n):
                                Wq[kk] = ((2 * kk - 1) * u_loc * Wq[kk - 1] - (kk - 1) * Wq[kk - 2]) / kk

                        for i in range(npos):
                            d0, d1, d2 = px[i, 0] - t * th0, px[i, 1] - t * th1, px[i, 2] - t * th2
                            r2 = d0 * d0 + d1 * d1 + d2 * d2
                            inv_r = 1.0 / np.sqrt(r2)
                            inv_r2 = inv_r * inv_r
                            kr = k * (r2 * inv_r)
                            ik = inv_r * inv_k
                            ik2 = ik * ik
                            phc = (np.cos(kr) - 1j * np.sin(kr)) * inv_r
                            fIc = (1.0 - ik2 - 1j * ik) * phc
                            fRc = (1.0 - 3.0 * ik2 - 3j * ik) * phc * inv_r2
                            fMc = (-ik - 1j) * phc * inv_r
                            for cl in range(ncls):
                                sg = sgn[cl, r]
                                u0 = sg * rx0 * np.conj(u[cl, m, 0])
                                u1 = sg * rx1 * np.conj(u[cl, m, 1])
                                u2 = sg * rx2 * np.conj(u[cl, m, 2])
                                du = fRc * (d0 * u0 + d1 * u1 + d2 * u2)
                                ge0, ge1, ge2 = fIc * u0 - du * d0, fIc * u1 - du * d1, fIc * u2 - du * d2
                                gm0 = fMc * (u1 * d2 - u2 * d1)
                                gm1 = fMc * (u2 * d0 - u0 * d2)
                                gm2 = fMc * (u0 * d1 - u1 * d0)
                                if 0 <= q_idx < n_mom:
                                    G_mom[ch, q_idx, cl, 0, i, 0] += ge0
                                    G_mom[ch, q_idx, cl, 0, i, 1] += ge1
                                    G_mom[ch, q_idx, cl, 0, i, 2] += ge2
                                    G_mom[ch, q_idx, cl, 1, i, 0] += gm0
                                    G_mom[ch, q_idx, cl, 1, i, 1] += gm1
                                    G_mom[ch, q_idx, cl, 1, i, 2] += gm2
                                elif q_idx >= 0:
                                    for kk in range(n):
                                        wk = W[q_idx, kk]
                                        Z_cart[ch, cl, 0, i, kk, 0] += wk * ge0
                                        Z_cart[ch, cl, 0, i, kk, 1] += wk * ge1
                                        Z_cart[ch, cl, 0, i, kk, 2] += wk * ge2
                                        Z_cart[ch, cl, 1, i, kk, 0] += wk * gm0
                                        Z_cart[ch, cl, 1, i, kk, 1] += wk * gm1
                                        Z_cart[ch, cl, 1, i, kk, 2] += wk * gm2
                                else:
                                    for kk in range(n):
                                        wk = Wq[kk]
                                        Z_cart[ch, cl, 0, i, kk, 0] += wk * ge0
                                        Z_cart[ch, cl, 0, i, kk, 1] += wk * ge1
                                        Z_cart[ch, cl, 0, i, kk, 2] += wk * ge2
                                        Z_cart[ch, cl, 1, i, kk, 0] += wk * gm0
                                        Z_cart[ch, cl, 1, i, kk, 1] += wk * gm1
                                        Z_cart[ch, cl, 1, i, kk, 2] += wk * gm2


def eval_line_adj(pts, refl, r_act, sgn, t_hat, mids, lam, d, n, k, H, u, out_pair, rot=None, w=None):
    """Exact conjugate transpose of eval_line_field, adding mode sensitivities into out_pair[cl, 2, ns, n, 3]"""
    mids = np.ascontiguousarray(np.atleast_2d(mids), dtype=float)
    w = np.ones((2, 1, len(mids), 3)) if w is None else np.ascontiguousarray(w)
    rot_arr = np.eye(3) if rot is None else np.ascontiguousarray(rot)
    tq, W, off, tlo, thi, ltab, g, gw = build_tree(lam, d, n, H, k)
    n_mom = int(off[min(511, len(off) - 1)])
    ncls, npos = sgn.shape[0], len(mids)
    nch = min(2, max(1, len(pts)))
    G_mom = np.zeros((nch, n_mom, ncls, 2, npos, 3), np.complex128)
    Z_cart = np.zeros((nch, ncls, 2, npos, n, 3), np.complex128)
    _eval_line_adj_kernel(pts, refl, np.ascontiguousarray(r_act, np.int64), sgn,
                          np.ascontiguousarray(t_hat, dtype=float), mids, tq, W, off, tlo, thi,
                          ltab, g, gw, lam, d, n, k, H, n_mom, u, G_mom, Z_cart)
    Z_tot = Z_cart.sum(axis=0) + np.einsum("pm,plkiq->lkimq", W[:n_mom], G_mom.sum(axis=0))
    out_pair += 2.0 * np.einsum("ksic,lkimq,cq->lksmc", w, Z_tot, rot_arr)

