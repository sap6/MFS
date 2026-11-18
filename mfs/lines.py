"""
Fields of the line sources of the cube (Section 5.1 and Appendix B of the
paper), assembled straight into the matrices of the parity classes.

A line at depth d carries a dipole density that is a polynomial in the stretched
variable u in [-1, 1], with t(u) = H tanh(λu/2) along the edge and λ = ln((2H - d)/d).
Each matrix entry is the integral of a dipole field against one Legendre mode, and
most of the work is shared between targets:
  * [-1, 1] is cut into a binary tree of panels whose nodes and weights against all
    modes are computed once per line (exact moments in t on the four coarsest levels,
    so a distant target needs only PM kernel values for the whole line);
  * a target descends the tree until a panel is no longer than its distance from it;
  * the electric and magnetic line share the kernel values, the nine lines of a tail
    are summed with their weights, and the eight mirror images of a target are summed
    with the signs of each parity class before anything is written.
"""
import numpy as np
from numba import njit, prange
from scipy.interpolate import BarycentricInterpolator

LMOM = 4          # the four coarsest levels use moments (panel widths 2, 1, 1/2, 1/4 in u)
PM, NG = 12, 8    # nodes of a moment panel and of a Gauss panel, for a line of one wavelength
LHARD = 40        # deepest level of the descent

_TREES = {}       # the tables of each line, built once by build_tree


def nodes(k, H):
    """
    Nodes per panel (moment panels, Gauss panels) for wavenumber k and lines of length up to 2H.
    The kernel phase along the coarsest panels grows like kH: 12 / 8 nodes at one wavelength
    (the cube fields change by 1e-8 with more), 20 / 10 at two and 28 / 12 at three (Appendix B)
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
    gm, _ = leg.leggauss(pm)
    gg, gwg = leg.leggauss(ng)
    gf, gwf = leg.leggauss(16)                           # fine rule for the exact moments
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
            for i in range(m):
                tn = 0.5 * (a[i] + b[i]) + 0.5 * (b[i] - a[i]) * gm
                nsub = max(8, int(np.ceil((uhi[i] - ulo[i]) / 0.03125)))   # sub-intervals of width at most 1/32 in u
                e = np.linspace(ulo[i], uhi[i], nsub + 1)
                u = (0.5 * (e[:-1] + e[1:])[:, None] + 0.5 * (e[1:] - e[:-1])[:, None] * gf).ravel()
                w = (0.5 * (e[1:] - e[:-1])[:, None] * gwf).ravel()
                w = w * H * 0.5 * lam / np.cosh(0.5 * lam * u) ** 2 / d
                L = BarycentricInterpolator(tn, np.eye(pm))(H * np.tanh(0.5 * lam * u))   # Lagrange polynomials, (fine, pm)
                V = leg.legvander(u, n - 1)                                   # (fine, n)
                ts.append(tn)
                Ws.append(L.T @ (V * w[:, None]))
                off[p + i + 1] = off[p + i] + pm
        else:
            # Gauss panels: ng nodes in u, weights = the modes times dt/du
            u = 0.5 * (ulo + uhi)[:, None] + 0.5 * (uhi - ulo)[:, None] * gg  # (m, ng)
            w = 0.5 * (uhi - ulo)[:, None] * gwg * H * 0.5 * lam / np.cosh(0.5 * lam * u) ** 2 / d
            ts.append((H * np.tanh(0.5 * lam * u)).ravel())
            Ws.append(leg.legvander(u.ravel(), n - 1) * w.ravel()[:, None])
            off[p + 1:p + m + 1] = off[p] + ng * np.arange(1, m + 1)
        p += m
    out = (np.ascontiguousarray(np.concatenate(ts)), np.ascontiguousarray(np.vstack(Ws)), off, tlo, thi, ltab, gg, gwg)
    _TREES[key] = out
    return out


@njit(fastmath=True, cache=True, inline="always")
def _dipoles(Tm, e, k, r, out_e, out_m):
    """
    Tested fields of unit electric and magnetic dipoles at distance r in the direction e:
      out_e[a, c] = Tm[a] · G e_c   with G = fI I - fR e e^T,
      out_m[a, c] = Tm[a] · (fm e × e_c)
    """
    kr = k * r
    ph = (np.cos(kr) + 1j * np.sin(kr)) / r
    fI = (1.0 + 1j / kr - 1.0 / (kr * kr)) * ph
    fR = (1.0 + 3j / kr - 3.0 / (kr * kr)) * ph
    fm = (1j - 1.0 / kr) * ph
    for a in range(Tm.shape[0]):
        te = Tm[a, 0] * e[0] + Tm[a, 1] * e[1] + Tm[a, 2] * e[2]
        for c in range(3):
            out_e[a, c] = fI * Tm[a, c] - fR * te * e[c]
        out_m[a, 0] = fm * (Tm[a, 1] * e[2] - Tm[a, 2] * e[1])
        out_m[a, 1] = fm * (Tm[a, 2] * e[0] - Tm[a, 0] * e[2])
        out_m[a, 2] = fm * (Tm[a, 0] * e[1] - Tm[a, 1] * e[0])


@njit(parallel=True, fastmath=True, cache=True)
def _tree_kernel(x, T, refl, ax, b, c, pos, rot, w, tq, W, off, tlo, thi, ltab, g, gw, lam, d, n, k, H, A, cmap, sgn):
    """
    For every target x[m], image R = refl[r] and class s, adds
      sgn[s, r] (rows T[m] R) · (field at R x[m] of column j)
    to A[s, m nc + a, cmap[s, j]] for every column j = ((kind ns + slot) n + mode) 3 + component
    with cmap[s, j] >= 0 (kind 0: electric, 1: magnetic; components in the frame rot)
    """
    M, nr, npos, ns, nc = x.shape[0], refl.shape[0], pos.shape[0], w.shape[1], T.shape[1]
    tmul = 1 if T.shape[0] > 1 else 0
    for m in prange(M):
        accr, acci = np.zeros((nr, 2, ns, nc, 3, n)), np.zeros((nr, 2, ns, nc, 3, n))   # real and imaginary parts
        S = np.zeros((2, ns, nc, 3), np.complex128)
        Ge, Gm = np.empty((nc, 3), np.complex128), np.empty((nc, 3), np.complex128)
        Tm, xm, e, Wq = np.empty((nc, 3)), np.empty(3), np.empty(3), np.empty(n)
        stack_l, stack_i = np.empty(2 * LHARD + 8, np.int64), np.empty(2 * LHARD + 8, np.int64)
        for r in range(nr):
            # The image of the target and of its rows
            for j in range(3):
                xm[j] = x[m, j] * refl[r, j]
                for a in range(nc):
                    Tm[a, j] = T[m * tmul, a, j] * refl[r, j]
            rho2 = 1e300                                       # squared distance to the nearest line, across the edge
            for i in range(npos):
                rho2 = min(rho2, (xm[b] - pos[i, 0]) ** 2 + (xm[c] - pos[i, 1]) ** 2)

            # Descend the tree: accept a panel when it is no longer than its distance from the target
            sp, stack_l[0], stack_i[0] = 1, 0, 0
            while sp > 0:
                sp -= 1
                lev, idx = stack_l[sp], stack_i[sp]
                wid = 2.0 / (1 << lev)
                p = (1 << lev) - 1 + idx if lev <= ltab else -1
                a0 = tlo[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * idx))
                b0 = thi[p] if p >= 0 else H * np.tanh(0.5 * lam * (-1.0 + wid * (idx + 1)))
                gap = max(a0 - xm[ax], 0.0, xm[ax] - b0)
                if b0 - a0 > np.sqrt(rho2 + gap * gap) and lev < LHARD:
                    stack_l[sp], stack_i[sp], stack_l[sp + 1], stack_i[sp + 1] = lev + 1, 2 * idx, lev + 1, 2 * idx + 1
                    sp += 2
                    continue
                for iq in range(off[p + 1] - off[p] if p >= 0 else g.shape[0]):
                    # The node t and its weights against the modes: tabulated, or (deep panels) by recurrence
                    if p >= 0:
                        t = tq[off[p] + iq]
                        for kk in range(n):
                            Wq[kk] = W[off[p] + iq, kk]
                    else:
                        u = -1.0 + wid * (idx + 0.5) + 0.5 * wid * g[iq]
                        ex = np.exp(lam * u)
                        t = H * (ex - 1.0) / (ex + 1.0)
                        Wq[0] = 0.5 * wid * gw[iq] * H * 2.0 * lam * ex / ((ex + 1.0) * (ex + 1.0)) / d
                        if n > 1:
                            Wq[1] = Wq[0] * u
                        for kk in range(2, n):
                            Wq[kk] = ((2 * kk - 1) * u * Wq[kk - 1] - (kk - 1) * Wq[kk - 2]) / kk

                    # Fields of the dipoles at t on every line, in the frame rot, summed with the weights w
                    S[:] = 0.0
                    for i in range(npos):
                        e[ax], e[b], e[c] = xm[ax] - t, xm[b] - pos[i, 0], xm[c] - pos[i, 1]
                        rr = np.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
                        e /= rr
                        _dipoles(Tm, e, k, rr, Ge, Gm)
                        for a in range(nc):
                            for cc in range(3):
                                le = Ge[a, 0] * rot[cc, 0] + Ge[a, 1] * rot[cc, 1] + Ge[a, 2] * rot[cc, 2]
                                lm = Gm[a, 0] * rot[cc, 0] + Gm[a, 1] * rot[cc, 1] + Gm[a, 2] * rot[cc, 2]
                                for s in range(ns):
                                    S[0, s, a, cc] += w[0, s, i, cc] * le
                                    S[1, s, a, cc] += w[1, s, i, cc] * lm
                    for kind in range(2):
                        for s in range(ns):
                            for a in range(nc):
                                for cc in range(3):
                                    sr, si = S[kind, s, a, cc].real, S[kind, s, a, cc].imag
                                    for kk in range(n):
                                        accr[r, kind, s, a, cc, kk] += sr * Wq[kk]
                                        acci[r, kind, s, a, cc, kk] += si * Wq[kk]

        # Sum the images with the signs of each class, and write each entry once
        for cl in range(sgn.shape[0]):
            for kind in range(2):
                for s in range(ns):
                    for kk in range(n):
                        for cc in range(3):
                            col = cmap[cl, ((kind * ns + s) * n + kk) * 3 + cc]
                            if col >= 0:
                                for a in range(nc):
                                    vr, vi = 0.0, 0.0
                                    for r in range(nr):
                                        vr += sgn[cl, r] * accr[r, kind, s, a, cc, kk]
                                        vi += sgn[cl, r] * acci[r, kind, s, a, cc, kk]
                                    A[cl, m * nc + a, col] += vr + 1j * vi


def tree_add(x, T, refl, ax, b, c, pos, lam, d, n, k, H, A, cmap, sgn, rot=None, w=None):
    """
    Adds the fields of an electric and a magnetic line density (or tail) to the class matrices A
      x, T : targets (M, 3) and their rows (M, nc, 3), or (1, nc, 3) for the same rows everywhere
      refl : the mirror images (images, 3); sgn: their signs in each class (classes, images)
      ax, b, c: the axis of the edge and the two others; pos: (lines, 2) their coordinates b, c
      lam, d, n: the stretching, the depth and the number of modes
      cmap : (classes, 2 slots n 3) the column of each unknown in each class, -1 if not kept
      rot  : the frame of the source components; w: (2, slots, lines, 3) weights of the lines
    """
    pos = np.ascontiguousarray(np.atleast_2d(pos), dtype=float)
    w = np.ones((2, 1, len(pos), 3)) if w is None else np.ascontiguousarray(w)
    tq, W, off, tlo, thi, ltab, g, gw = build_tree(lam, d, n, H, k)
    _tree_kernel(np.ascontiguousarray(x), np.ascontiguousarray(T), refl, ax, b, c, pos,
                 np.eye(3) if rot is None else np.ascontiguousarray(rot), w, tq, W, off, tlo, thi, ltab, g, gw,
                 lam, d, n, k, H, A, np.ascontiguousarray(cmap), sgn)


@njit(parallel=True, fastmath=True, cache=True)
def points_add(x, T, refl, pe, pm, k, A, col0, sgn):
    """As _tree_kernel, for electric point dipoles at pe and magnetic ones at pm (columns col0[s] + 3 j + c in class s)"""
    M, nc, ne = x.shape[0], T.shape[1], pe.shape[0]
    ncls, nr = sgn.shape
    tmul = 1 if T.shape[0] > 1 else 0
    for m in prange(M):
        acc = np.zeros((ncls, nc, 3), np.complex128)
        Ge, Gm = np.empty((nc, 3), np.complex128), np.empty((nc, 3), np.complex128)
        Tm, e = np.empty((nc, 3)), np.empty(3)
        for j in range(ne + pm.shape[0]):
            acc[:] = 0.0
            for r in range(nr):
                y = pe[j] if j < ne else pm[j - ne]
                for q in range(3):
                    e[q] = x[m, q] * refl[r, q] - y[q]
                    for a in range(nc):
                        Tm[a, q] = T[m * tmul, a, q] * refl[r, q]
                rr = np.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
                e /= rr
                _dipoles(Tm, e, k, rr, Ge, Gm)
                for s in range(ncls):
                    for a in range(nc):
                        for cc in range(3):
                            acc[s, a, cc] += sgn[s, r] * (Ge[a, cc] if j < ne else Gm[a, cc])
            for s in range(ncls):
                for a in range(nc):
                    for cc in range(3):
                        A[s, m * nc + a, col0[s] + 3 * j + cc] += acc[s, a, cc]
