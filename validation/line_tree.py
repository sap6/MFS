"""Fast assembly of the field of line dipole distributions (cube edges).

Same integrals as line_kernel.line_block, organised so that almost all of the
work is shared between targets:

  * the stretched variable u in [-1, 1] of a line is cut into a binary tree of
    panels. Every panel carries a fixed set of nodes and, for each node, the
    weights against all Legendre modes of the dipole density, computed once
    per (depth, number of modes):
      - levels 0 .. LMOM-1 (panels wide in u, where the modes and the map
        t(u) are not smooth on the scale of the panel): PM Gauss nodes in the
        physical variable t and exact moments of the Lagrange polynomials
        against the modes, so that a distant target needs PM kernel values
        for the whole line;
      - levels LMOM .. ltab: NG Gauss nodes in u;
      - deeper levels (very close targets only) are generated on the fly.
  * a target descends the tree until a panel is short compared with its
    distance from the target (length <= alpha * distance), and uses the nodes
    of the panels it accepts;
  * the electric and the magnetic line at the same position share the kernel
    values, and the lines of a singular tail (same parametrisation, different
    depths) are summed with their weights before the modes are accumulated.

`tree_block` returns blocks identical in layout to Line.block /
TailLine.block of pec_cube_lightning3d.py; test_line_tree.py compares the two.
"""
import numpy as np
from numba import njit, prange

LMOM = 4          # levels 0..3 use moments (panel widths 2, 1, 0.5, 0.25 in u)
PM = 20           # nodes of a moment panel
NG = 10           # Gauss nodes of a deeper panel
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


def build_tree(lam, d, n, H, ltab):
    """Node positions t, weights W (nodes, n) and panel table of one line
    parametrisation. Panels are stored level by level (heap order)."""
    key = (round(lam, 12), round(d, 15), n, ltab)
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
def _tree_kernel(x, T, ax, b, c, pos, rot, we, wm, do_e, do_m, tq, W, off, tlo, thi, ltab, lmom,
                 alpha_m, alpha_g, lhard, g, gw, lam, d, n, k, H, A_e, cmap_e, A_m, cmap_m, sgn):
    """Adds sgn[s] * (field of column j at row (target, component)) to
    A[target * nc + component, cmap[s, j]] for every class s and every column
    j = (slot * n + mode) * 3 + source component with cmap[s, j] >= 0. T is
    (M, nc, 3) (one matrix of row vectors per target) or (1, nc, 3) (the same
    for all targets)."""
    M = x.shape[0]
    npos = pos.shape[0]
    nse, nsm = we.shape[0], wm.shape[0]
    nc = T.shape[1]
    tmul = 1 if T.shape[0] > 1 else 0
    ncls = sgn.shape[0]
    ng = g.shape[0]
    deep = np.zeros(M, np.int64)
    for m in prange(M):
        xa, xb, xc = x[m, ax], x[m, b], x[m, c]
        rho2 = 1e300
        for i in range(npos):
            r2 = (xb - pos[i, 0]) ** 2 + (xc - pos[i, 1]) ** 2
            if r2 < rho2:
                rho2 = r2
        Tm = np.zeros((3, 3))
        mt = m * tmul
        for a in range(nc):
            for j in range(3):
                Tm[a, j] = T[mt, a, j]
        ar_e = np.zeros((nse, nc * 3, n))
        ai_e = np.zeros((nse, nc * 3, n))
        ar_m = np.zeros((nsm, nc * 3, n))
        ai_m = np.zeros((nsm, nc * 3, n))
        Se = np.zeros((nse, nc, 3), np.complex128)
        Sm = np.zeros((nsm, nc, 3), np.complex128)
        Ge = np.zeros((3, 3), np.complex128)
        Gm = np.zeros((3, 3), np.complex128)
        Wq = np.empty(n)
        dv = np.empty(3)
        stack_l = np.empty(2 * lhard + 8, np.int64)
        stack_i = np.empty(2 * lhard + 8, np.int64)
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
                                    ar_e[ie, row, kk] += sr * Wq[kk]
                                    ai_e[ie, row, kk] += si * Wq[kk]
                if do_m:
                    for ie in range(nsm):
                        for a in range(nc):
                            for cc in range(3):
                                sr, si = Sm[ie, a, cc].real, Sm[ie, a, cc].imag
                                row = a * 3 + cc
                                for kk in range(n):
                                    ar_m[ie, row, kk] += sr * Wq[kk]
                                    ai_m[ie, row, kk] += si * Wq[kk]
        for cl in range(ncls):
            sg = sgn[cl]
            if do_e:
                for ie in range(nse):
                    for kk in range(n):
                        for cc in range(3):
                            col = cmap_e[cl, (ie * n + kk) * 3 + cc]
                            if col >= 0:
                                for a in range(nc):
                                    A_e[m * nc + a, col] += sg * (ar_e[ie, a * 3 + cc, kk] + 1j * ai_e[ie, a * 3 + cc, kk])
            if do_m:
                for ie in range(nsm):
                    for kk in range(n):
                        for cc in range(3):
                            col = cmap_m[cl, (ie * n + kk) * 3 + cc]
                            if col >= 0:
                                for a in range(nc):
                                    A_m[m * nc + a, col] += sg * (ar_m[ie, a * 3 + cc, kk] + 1j * ai_m[ie, a * 3 + cc, kk])
    return deep.max()


_EYE = np.eye(3)
_ONE = np.ones((1, 1, 3))
_DUMMY = np.zeros((1, 1), complex)
_NOMAP = np.full((1, 1), -1, np.int64)
_PLUS = np.ones(1)


def _ltab(lam, d):
    """Tabulate down to panels of length about d in the middle of the line."""
    return int(min(max(LMOM + 2, np.ceil(np.log2(lam / (2.0 * d)))), 12))


def tree_add(x, T, ax, b, c, pos, lam, d, n, k, H, A, cmap_e, cmap_m, sgn, rot=None, we=None, wm=None):
    """Adds the fields of the electric (cmap_e given) and/or magnetic (cmap_m
    given) line distributions at the positions `pos` ((npos, 2): coordinates
    b, c), all with the parametrisation (lam, d, n), to the matrix A (rows:
    target * nc + component). `cmap_*` ((classes, slots * n * 3), -1 = column
    not kept) and `sgn` (classes,) say where each column goes and with which
    sign; `we`, `wm`: (slots, npos, 3) weights of the positions per exponent
    slot and local source component; `rot` maps Cartesian source components
    to the local frame. T: (M, nc, 3) or (1, nc, 3) row vectors."""
    pos = np.ascontiguousarray(np.atleast_2d(pos), dtype=float)
    ltab = _ltab(lam, d)
    tq, W, off, tlo, thi, g, gw = build_tree(lam, d, n, H, ltab)
    do_e, do_m = cmap_e is not None, cmap_m is not None
    one = np.ones((1, pos.shape[0], 3))
    return _tree_kernel(np.ascontiguousarray(x), np.ascontiguousarray(T), ax, b, c, pos,
                        _EYE if rot is None else np.ascontiguousarray(rot),
                        np.ascontiguousarray(we) if we is not None else one,
                        np.ascontiguousarray(wm) if wm is not None else one,
                        do_e, do_m, tq, W, off, tlo, thi, ltab, LMOM, ALPHA_M, ALPHA_G, LHARD, g, gw, lam, d, n, k, H,
                        A, cmap_e if do_e else _NOMAP, A, cmap_m if do_m else _NOMAP, sgn)


def tree_block(x, T, ax, b, c, pos, lam, d, n, k, H, rot=None, we=None, wm=None):
    """Field blocks of the electric (we given) and/or magnetic (wm given)
    line distributions; see tree_add. Returns (out_e, out_m), each
    (M, nc, slots * n, 3) or None."""
    pos = np.ascontiguousarray(np.atleast_2d(pos), dtype=float)
    ltab = _ltab(lam, d)
    tq, W, off, tlo, thi, g, gw = build_tree(lam, d, n, H, ltab)
    Tm = np.eye(3)[None] if T is None else np.ascontiguousarray(T)
    nc = Tm.shape[1]
    do_e, do_m = we is not None, wm is not None
    out_e = np.zeros((len(x) * nc, we.shape[0] * n * 3), complex) if do_e else _DUMMY
    out_m = np.zeros((len(x) * nc, wm.shape[0] * n * 3), complex) if do_m else _DUMMY
    _tree_kernel(np.ascontiguousarray(x), Tm, ax, b, c, pos, _EYE if rot is None else np.ascontiguousarray(rot),
                 np.ascontiguousarray(we) if do_e else _ONE, np.ascontiguousarray(wm) if do_m else _ONE,
                 do_e, do_m, tq, W, off, tlo, thi, ltab, LMOM, ALPHA_M, ALPHA_G, LHARD, g, gw, lam, d, n, k, H,
                 out_e, np.arange(out_e.shape[1])[None, :] if do_e else _NOMAP,
                 out_m, np.arange(out_m.shape[1])[None, :] if do_m else _NOMAP, _PLUS)
    return (out_e.reshape(len(x), nc, -1, 3) if do_e else None), (out_m.reshape(len(x), nc, -1, 3) if do_m else None)
