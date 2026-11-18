"""
Fields of the line sources of the cube (Section 5.1 and Appendix B of the
paper), assembled straight into the matrices of the parity classes.

A line at depth d carries a dipole density that is a polynomial in the
stretched variable u in [-1, 1], with t(u) = H tanh(λu/2) along the edge and
λ = ln((2H - d)/d). Each matrix entry is the integral of a dipole field
against one Legendre mode. Most of the work is shared between targets:

  * [-1, 1] is cut into a binary tree of panels with fixed nodes and, for
    each node, weights against all modes, computed once per line. The four
    coarsest levels use PM Gauss nodes in t and exact moments of the Lagrange
    polynomials, so a distant target needs PM kernel values for the whole
    line; the finer levels use NG Gauss nodes in u.
  * a target descends the tree until a panel is no longer than its distance
    from the target, and uses the nodes of the panels it accepts.
  * the electric and the magnetic line at the same position share the kernel
    values, and the nine lines of a singular tail are summed with their
    weights before the modes are accumulated.
  * the eight mirror images of a target are summed with the signs of each
    parity class before anything is written.
"""
import numpy as np
from numba import njit, prange

LMOM = 4          # levels 0..3 use moments (panel widths 2, 1, 0.5, 0.25 in u)
PM = 12           # nodes of a moment panel for a line of one wavelength (20 gives the same cube fields to 1e-8)
NG = 8            # Gauss nodes of a deeper panel, likewise (10: the same to 1e-8)
ALPHA_M = 1.0     # accept a moment panel if its length <= ALPHA_M * distance
ALPHA_G = 1.0     # the same for Gauss panels
LHARD = 40        # absolute depth limit of the descent

_TREES = {}


def _lagrange(tn, t):
    """Values of the Lagrange polynomials on nodes tn at points t: (len(t), len(tn))."""
    w = np.array([1.0 / np.prod(tn[j] - np.delete(tn, j)) for j in range(len(tn))])
    diff = t[:, None] - tn[None, :]
    exact = diff == 0.0
    diff[exact] = 1.0
    L = w[None, :] / diff
    L /= L.sum(axis=1, keepdims=True)
    rows = exact.any(axis=1)
    L[rows] = exact[rows].astype(float)
    return L


def nodes(k, H):
    """Nodes per panel (moment panels, Gauss panels) for wavenumber k and
    lines of length up to 2 H. A panel is accepted when its length is at most
    its distance from the target, so the kernel phase varies along the coarsest
    panels in proportion to k H, and the rule must grow with it: 12 / 8 nodes
    are converged (1e-8 in the cube fields) at one wavelength, but give an
    error of 5e-3 at three wavelengths, where 20 / 10 are needed."""
    s = max(0.0, k * 2 * H / (2 * np.pi) - 1.0)
    return int(np.ceil(PM + 8 * s)), int(np.ceil(NG + 2 * s))


def build_tree(lam, d, n, H, ltab, k=2 * np.pi):
    """Node positions t, weights W (nodes, n) and panel table of one line
    parametrisation. Panels are stored level by level (heap order)."""
    PM, NG = nodes(k, H)
    key = (round(lam, 12), round(d, 15), n, ltab, PM, NG)
    if key in _TREES:
        return _TREES[key]
    leg = np.polynomial.legendre
    gm, gwm = leg.leggauss(PM)
    gg, gwg = leg.leggauss(NG)
    gf, gwf = leg.leggauss(16)
    npan = 2 ** (ltab + 1) - 1
    off = np.zeros(npan + 1, np.int64)
    tlo, thi = np.empty(npan), np.empty(npan)
    ts, Ws = [], []
    p = 0
    for lev in range(ltab + 1):
        m = 2 ** lev
        ulo = -1.0 + 2.0 * np.arange(m) / m
        uhi = ulo + 2.0 / m
        a, b = H * np.tanh(0.5 * lam * ulo), H * np.tanh(0.5 * lam * uhi)
        tlo[p:p + m], thi[p:p + m] = a, b
        if lev < LMOM:
            for i in range(m):
                tn = 0.5 * (a[i] + b[i]) + 0.5 * (b[i] - a[i]) * gm
                nsub = max(8, int(np.ceil((uhi[i] - ulo[i]) / 0.03125)))
                e = np.linspace(ulo[i], uhi[i], nsub + 1)
                u = (0.5 * (e[:-1] + e[1:])[:, None] + 0.5 * (e[1:] - e[:-1])[:, None] * gf).ravel()
                w = (0.5 * (e[1:] - e[:-1])[:, None] * gwf).ravel()
                w = w * H * 0.5 * lam / np.cosh(0.5 * lam * u) ** 2 / d
                L = _lagrange(tn, H * np.tanh(0.5 * lam * u))                 # (fine, PM)
                V = leg.legvander(u, n - 1)                                   # (fine, n)
                ts.append(tn)
                Ws.append(L.T @ (V * w[:, None]))
                off[p + i + 1] = off[p + i] + PM
        else:
            u = 0.5 * (ulo + uhi)[:, None] + 0.5 * (uhi - ulo)[:, None] * gg  # (m, NG)
            w = 0.5 * (uhi - ulo)[:, None] * gwg * H * 0.5 * lam / np.cosh(0.5 * lam * u) ** 2 / d
            ts.append((H * np.tanh(0.5 * lam * u)).ravel())
            Ws.append(leg.legvander(u.ravel(), n - 1) * w.ravel()[:, None])
            off[p + 1:p + m + 1] = off[p] + NG * np.arange(1, m + 1)
        p += m
    out = (np.ascontiguousarray(np.concatenate(ts)), np.ascontiguousarray(np.vstack(Ws)), off, tlo, thi, gg, gwg)
    _TREES[key] = out
    return out


@njit(parallel=True, fastmath=True, cache=True)
def _tree_kernel(x, T, refl, ax, b, c, pos, rot, we, wm, do_e, do_m, tq, W, off, tlo, thi, ltab, lmom,
                 alpha_m, alpha_g, lhard, g, gw, lam, d, n, k, H, A_e, cmap_e, A_m, cmap_m, sgn):
    """For every target x[m], every image R = refl[r] (a vector of signs) and
    every class s, adds
        sgn[s, r] * (row vectors T[m] * R) . (field at x[m] * R of column j)
    to A[s, m * nc + component, cmap[s, j]], for the columns
    j = (slot * n + mode) * 3 + source component with cmap[s, j] >= 0. The
    images are summed before anything is written, so each entry of A is
    touched once. T is (M, nc, 3) or (1, nc, 3) (the same for all targets)."""
    M = x.shape[0]
    npos = pos.shape[0]
    nr = refl.shape[0]
    nse, nsm = we.shape[0], wm.shape[0]
    nc = T.shape[1]
    tmul = 1 if T.shape[0] > 1 else 0
    ncls = sgn.shape[0]
    ng = g.shape[0]
    deep = np.zeros(M, np.int64)
    for m in prange(M):
        mt = m * tmul
        ar_e = np.zeros((nr, nse, nc * 3, n))
        ai_e = np.zeros((nr, nse, nc * 3, n))
        ar_m = np.zeros((nr, nsm, nc * 3, n))
        ai_m = np.zeros((nr, nsm, nc * 3, n))
        Se = np.zeros((nse, nc, 3), np.complex128)
        Sm = np.zeros((nsm, nc, 3), np.complex128)
        Ge = np.zeros((3, 3), np.complex128)
        Gm = np.zeros((3, 3), np.complex128)
        Tm = np.zeros((3, 3))
        Wq = np.empty(n)
        dv = np.empty(3)
        xm = np.empty(3)
        stack_l = np.empty(2 * lhard + 8, np.int64)
        stack_i = np.empty(2 * lhard + 8, np.int64)
        for r in range(nr):
            for j in range(3):
                xm[j] = x[m, j] * refl[r, j]
                for a in range(nc):
                    Tm[a, j] = T[mt, a, j] * refl[r, j]
            xa, xb, xc = xm[ax], xm[b], xm[c]
            rho2 = 1e300
            for i in range(npos):
                r2 = (xb - pos[i, 0]) ** 2 + (xc - pos[i, 1]) ** 2
                if r2 < rho2:
                    rho2 = r2
            sp = 1
            stack_l[0] = 0
            stack_i[0] = 0
            while sp > 0:
                sp -= 1
                lev, idx = stack_l[sp], stack_i[sp]
                wid = 2.0 / (1 << lev)
                if lev <= ltab:
                    p = (1 << lev) - 1 + idx
                    a0, b0 = tlo[p], thi[p]
                else:
                    p = -1
                    u_lo = -1.0 + wid * idx
                    a0, b0 = H * np.tanh(0.5 * lam * u_lo), H * np.tanh(0.5 * lam * (u_lo + wid))
                gap = max(a0 - xa, 0.0, xa - b0)
                dist = np.sqrt(rho2 + gap * gap)
                if (b0 - a0) > (alpha_m if lev < lmom else alpha_g) * dist and lev < lhard:
                    stack_l[sp], stack_i[sp] = lev + 1, 2 * idx
                    stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx + 1
                    sp += 2
                    continue
                if p >= 0:
                    q0, nq = off[p], off[p + 1] - off[p]
                else:
                    q0, nq = 0, ng
                    if lev > deep[m]:
                        deep[m] = lev
                for iq in range(nq):
                    if p >= 0:
                        t = tq[q0 + iq]
                        for kk in range(n):
                            Wq[kk] = W[q0 + iq, kk]
                    else:
                        u = -1.0 + wid * (idx + 0.5) + 0.5 * wid * g[iq]
                        ex = np.exp(lam * u)
                        t = H * (ex - 1.0) / (ex + 1.0)
                        wt = 0.5 * wid * gw[iq] * H * 2.0 * lam * ex / ((ex + 1.0) * (ex + 1.0)) / d
                        Wq[0] = wt
                        if n > 1:
                            Wq[1] = wt * u
                        for kk in range(2, n):
                            Wq[kk] = ((2 * kk - 1) * u * Wq[kk - 1] - (kk - 1) * Wq[kk - 2]) / kk
                    for ie in range(nse):
                        for a in range(nc):
                            for cc in range(3):
                                Se[ie, a, cc] = 0.0
                    for ie in range(nsm):
                        for a in range(nc):
                            for cc in range(3):
                                Sm[ie, a, cc] = 0.0
                    dv[ax] = xa - t
                    for i in range(npos):
                        dv[b] = xb - pos[i, 0]
                        dv[c] = xc - pos[i, 1]
                        rr = np.sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2])
                        kr = k * rr
                        ir = 1.0 / rr
                        ikr = 1.0 / kr
                        ph = (np.cos(kr) + 1j * np.sin(kr)) * ir
                        e0, e1, e2 = dv[0] * ir, dv[1] * ir, dv[2] * ir
                        if do_e:
                            fI = (1.0 + 1j * ikr - ikr * ikr) * ph
                            fR = (1.0 + 3j * ikr - 3.0 * ikr * ikr) * ph
                            Ge[0, 0] = fI - fR * e0 * e0; Ge[0, 1] = -fR * e0 * e1; Ge[0, 2] = -fR * e0 * e2
                            Ge[1, 0] = Ge[0, 1]; Ge[1, 1] = fI - fR * e1 * e1; Ge[1, 2] = -fR * e1 * e2
                            Ge[2, 0] = Ge[0, 2]; Ge[2, 1] = Ge[1, 2]; Ge[2, 2] = fI - fR * e2 * e2
                            for a in range(nc):
                                g0 = Tm[a, 0] * Ge[0, 0] + Tm[a, 1] * Ge[1, 0] + Tm[a, 2] * Ge[2, 0]
                                g1 = Tm[a, 0] * Ge[0, 1] + Tm[a, 1] * Ge[1, 1] + Tm[a, 2] * Ge[2, 1]
                                g2 = Tm[a, 0] * Ge[0, 2] + Tm[a, 1] * Ge[1, 2] + Tm[a, 2] * Ge[2, 2]
                                for cc in range(3):
                                    loc = g0 * rot[cc, 0] + g1 * rot[cc, 1] + g2 * rot[cc, 2]
                                    for ie in range(nse):
                                        Se[ie, a, cc] += we[ie, i, cc] * loc
                        if do_m:
                            f = (1j - ikr) * ph
                            Gm[0, 0] = 0.0; Gm[0, 1] = -f * e2; Gm[0, 2] = f * e1
                            Gm[1, 0] = f * e2; Gm[1, 1] = 0.0; Gm[1, 2] = -f * e0
                            Gm[2, 0] = -f * e1; Gm[2, 1] = f * e0; Gm[2, 2] = 0.0
                            for a in range(nc):
                                g0 = Tm[a, 0] * Gm[0, 0] + Tm[a, 1] * Gm[1, 0] + Tm[a, 2] * Gm[2, 0]
                                g1 = Tm[a, 0] * Gm[0, 1] + Tm[a, 1] * Gm[1, 1] + Tm[a, 2] * Gm[2, 1]
                                g2 = Tm[a, 0] * Gm[0, 2] + Tm[a, 1] * Gm[1, 2] + Tm[a, 2] * Gm[2, 2]
                                for cc in range(3):
                                    loc = g0 * rot[cc, 0] + g1 * rot[cc, 1] + g2 * rot[cc, 2]
                                    for ie in range(nsm):
                                        Sm[ie, a, cc] += wm[ie, i, cc] * loc
                    if do_e:
                        for ie in range(nse):
                            for a in range(nc):
                                for cc in range(3):
                                    sr, si = Se[ie, a, cc].real, Se[ie, a, cc].imag
                                    row = a * 3 + cc
                                    for kk in range(n):
                                        ar_e[r, ie, row, kk] += sr * Wq[kk]
                                        ai_e[r, ie, row, kk] += si * Wq[kk]
                    if do_m:
                        for ie in range(nsm):
                            for a in range(nc):
                                for cc in range(3):
                                    sr, si = Sm[ie, a, cc].real, Sm[ie, a, cc].imag
                                    row = a * 3 + cc
                                    for kk in range(n):
                                        ar_m[r, ie, row, kk] += sr * Wq[kk]
                                        ai_m[r, ie, row, kk] += si * Wq[kk]
        for cl in range(ncls):
            if do_e:
                for ie in range(nse):
                    for kk in range(n):
                        for cc in range(3):
                            col = cmap_e[cl, (ie * n + kk) * 3 + cc]
                            if col >= 0:
                                for a in range(nc):
                                    vr, vi = 0.0, 0.0
                                    for r in range(nr):
                                        vr += sgn[cl, r] * ar_e[r, ie, a * 3 + cc, kk]
                                        vi += sgn[cl, r] * ai_e[r, ie, a * 3 + cc, kk]
                                    A_e[cl, m * nc + a, col] += vr + 1j * vi
            if do_m:
                for ie in range(nsm):
                    for kk in range(n):
                        for cc in range(3):
                            col = cmap_m[cl, (ie * n + kk) * 3 + cc]
                            if col >= 0:
                                for a in range(nc):
                                    vr, vi = 0.0, 0.0
                                    for r in range(nr):
                                        vr += sgn[cl, r] * ar_m[r, ie, a * 3 + cc, kk]
                                        vi += sgn[cl, r] * ai_m[r, ie, a * 3 + cc, kk]
                                    A_m[cl, m * nc + a, col] += vr + 1j * vi
    return deep.max()


_EYE = np.eye(3)
_ONE = np.ones((1, 1, 3))
_NOMAP = np.full((1, 1), -1, np.int64)


def _ltab(lam, d):
    """Tabulate down to panels of length about d in the middle of the line."""
    return int(min(max(LMOM + 2, np.ceil(np.log2(lam / (2.0 * d)))), 12))


def tree_add(x, T, refl, ax, b, c, pos, lam, d, n, k, H, A, cmap_e, cmap_m, sgn, rot=None, we=None, wm=None):
    """Adds the fields of the electric (cmap_e given) and/or magnetic (cmap_m
    given) line distributions at the positions `pos` ((npos, 2): coordinates
    b, c), all with the parametrisation (lam, d, n), to the class matrices
    A (classes, rows, columns; row = target * nc + component), summed over
    the mirror images `refl` ((images, 3) signs) with the signs `sgn`
    (classes, images). `cmap_*` ((classes, slots * n * 3), -1 = column not
    kept) say where each column goes; `we`, `wm`: (slots, npos, 3) weights of
    the positions per exponent slot and local source component; `rot` maps
    Cartesian source components to the local frame. T: (M, nc, 3) or
    (1, nc, 3) row vectors."""
    pos = np.ascontiguousarray(np.atleast_2d(pos), dtype=float)
    ltab = _ltab(lam, d)
    tq, W, off, tlo, thi, g, gw = build_tree(lam, d, n, H, ltab, k)
    do_e, do_m = cmap_e is not None, cmap_m is not None
    one = np.ones((1, pos.shape[0], 3))
    return _tree_kernel(np.ascontiguousarray(x), np.ascontiguousarray(T), refl, ax, b, c, pos,
                        _EYE if rot is None else np.ascontiguousarray(rot),
                        np.ascontiguousarray(we) if we is not None else one,
                        np.ascontiguousarray(wm) if wm is not None else one,
                        do_e, do_m, tq, W, off, tlo, thi, ltab, LMOM, ALPHA_M, ALPHA_G, LHARD, g, gw, lam, d, n, k, H,
                        A, cmap_e if do_e else _NOMAP, A, cmap_m if do_m else _NOMAP, sgn)


@njit(parallel=True, fastmath=True, cache=True)
def points_add(x, T, refl, pe, pm, k, A, col0, sgn):
    """As _tree_kernel, for electric point dipoles at pe and magnetic point
    dipoles at pm (three components each; all components kept, columns
    col0[s] + 3 j + c in class s)."""
    M, nc = x.shape[0], T.shape[1]
    tmul = 1 if T.shape[0] > 1 else 0
    ne, nm = pe.shape[0], pm.shape[0]
    ncls, nr = sgn.shape
    for m in prange(M):
        mt = m * tmul
        G = np.zeros((3, 3), np.complex128)
        acc = np.zeros((ncls, nc, 3), np.complex128)
        for j in range(ne + nm):
            for s in range(ncls):
                for a in range(nc):
                    for cc in range(3):
                        acc[s, a, cc] = 0.0
            for r in range(nr):
                if j < ne:
                    y0, y1, y2 = pe[j, 0], pe[j, 1], pe[j, 2]
                else:
                    y0, y1, y2 = pm[j - ne, 0], pm[j - ne, 1], pm[j - ne, 2]
                d0, d1, d2 = x[m, 0] * refl[r, 0] - y0, x[m, 1] * refl[r, 1] - y1, x[m, 2] * refl[r, 2] - y2
                rr = np.sqrt(d0 * d0 + d1 * d1 + d2 * d2)
                kr = k * rr
                ir, ikr = 1.0 / rr, 1.0 / kr
                ph = (np.cos(kr) + 1j * np.sin(kr)) * ir
                e0, e1, e2 = d0 * ir, d1 * ir, d2 * ir
                if j < ne:
                    fI = (1.0 + 1j * ikr - ikr * ikr) * ph
                    fR = (1.0 + 3j * ikr - 3.0 * ikr * ikr) * ph
                    G[0, 0] = fI - fR * e0 * e0; G[0, 1] = -fR * e0 * e1; G[0, 2] = -fR * e0 * e2
                    G[1, 0] = G[0, 1]; G[1, 1] = fI - fR * e1 * e1; G[1, 2] = -fR * e1 * e2
                    G[2, 0] = G[0, 2]; G[2, 1] = G[1, 2]; G[2, 2] = fI - fR * e2 * e2
                else:
                    f = (1j - ikr) * ph
                    G[0, 0] = 0.0; G[0, 1] = -f * e2; G[0, 2] = f * e1
                    G[1, 0] = f * e2; G[1, 1] = 0.0; G[1, 2] = -f * e0
                    G[2, 0] = -f * e1; G[2, 1] = f * e0; G[2, 2] = 0.0
                for a in range(nc):
                    for cc in range(3):
                        v = (T[mt, a, 0] * refl[r, 0] * G[0, cc] + T[mt, a, 1] * refl[r, 1] * G[1, cc]
                             + T[mt, a, 2] * refl[r, 2] * G[2, cc])
                        for s in range(ncls):
                            acc[s, a, cc] += sgn[s, r] * v
            for s in range(ncls):
                for a in range(nc):
                    for cc in range(3):
                        A[s, m * nc + a, col0[s] + 3 * j + cc] += acc[s, a, cc]
