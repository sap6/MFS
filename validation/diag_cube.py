"""Where does the cube boundary residual sit? Solves one configuration per
argument group and prints, on an independent test grid, the RMS residual by
distance to the nearest edge (rows) and distance to the nearest vertex along
that edge (columns), plus the optical-theorem mismatch.

Usage: python3 diag_cube.py J "key=value ..." ["key=value ..." ...]
Keys: any Sources keyword (gamma, n0, n_base, vrays, tail_n, q, ...), col_p, col_q, col_extra, beta.
"""
import sys
import time

import numpy as np

import pec_cube_lightning3d as m

RZ = [0, 1e-3, 1e-2, 3e-2, 0.1, 0.5]
TZ = [0, 1e-2, 3e-2, 0.1, 0.25, 0.5]


def run(J, **kw):
    col = {k: kw.pop(k) for k in ("col_p", "col_q", "col_extra", "beta") if k in kw}
    col_p, col_q, col_extra, beta = (col.get("col_p", 3), col.get("col_q", 0.5), col.get("col_extra", 1),
                                     col.get("beta", 0.5))
    src = m.SymSources(J, q=kw.pop("q", 0.5), kinds="em", **kw)
    lv = int(np.ceil(np.log(src.depths[-1] / (2 * m.H)) / np.log(col_q))) + col_extra
    c, info = m.solve(src, lv, int(col_p), col_q, beta)
    ext, sca = m.cross_sections(src, c)
    bem = m.load_bem()
    E = (src.matrix(bem["pts"]) @ c).reshape(-1, 3)
    x, w, T = m.cube_surface(int(lv * np.log(col_q) / np.log(0.45)) + 1, 4, 0.45)
    sel = m.quarter(x)
    x, w, T = x[sel], w[sel], T[sel]
    rhs = -np.einsum("mab,mb->ma", T, m.incident(x))
    m2 = np.sum(np.abs((src.matrix(x, T) @ c).reshape(-1, 2) - rhs) ** 2, axis=1)
    d = np.sort(m.H - np.abs(x), axis=1)            # 0 (on the face), edge distance, along-edge distance
    rho, tau = d[:, 1], d[:, 2]
    print(f"unk={info['unknowns']} rows={info['rows']} asm={info['t_assemble']:.0f}s solve={info['t_solve']:.0f}s "
          f"vs_bem={np.linalg.norm(E - bem['E']) / np.linalg.norm(bem['E']):.2e} "
          f"sigma_ext={ext:.5f} sigma_sca={sca:.5f} optical={abs(ext - sca) / sca:.1e}")
    print("   rho \\ tau  " + " ".join(f"<{t:<7.2g}" for t in TZ[1:]) + "  | weighted share")
    tot = np.sum(w * rho * m2)
    for lo, hi in zip(RZ[:-1], RZ[1:]):
        line = f"   <{hi:<8.2g} "
        for tl, th in zip(TZ[:-1], TZ[1:]):
            s = (rho >= lo) & (rho < hi) & (tau >= tl) & (tau < th)
            line += f" {np.sqrt(np.sum(w[s] * m2[s]) / np.sum(w[s])):7.1e}" if s.any() else "    -   "
        s = (rho >= lo) & (rho < hi)
        print(line + f"  | {np.sum(w[s] * rho[s] * m2[s]) / tot:5.2f}", flush=True)
    # Rough far-field effect of the residual in each zone: the residual
    # tangential field acts as a magnetic current M = E_res x n on the
    # surface; near a flat conductor it radiates like 2M in free space.
    res = (src.matrix(x, T) @ c).reshape(-1, 2) - rhs
    Ev = np.einsum("ma,mac->mc", res, T)                       # residual as a 3-vector
    nrm = np.cross(T[:, 0], T[:, 1])
    nrm *= np.sign(np.einsum("mc,mc->m", nrm, x))[:, None]      # outward
    M = 2 * np.cross(Ev, nrm) * w[:, None]
    dJ = src.depths[-1]
    zones = {"vertex (tau<4d)": tau < 4 * dJ, "edge, rho<d": (tau >= 4 * dJ) & (rho < dJ),
             "edge, d<rho<0.1": (tau >= 4 * dJ) & (rho >= dJ) & (rho < 0.1), "interior": (tau >= 4 * dJ) & (rho >= 0.1),
             "all": np.ones(len(x), bool)}
    pts = bem["pts"]
    out = "   far-field effect of the residual by zone (relative to |E_sc|): "
    for name, sel in zones.items():
        Ez = 0
        for R, sig in m.SymSources.GROUP:                       # images of the quarter surface
            G = m.magnetic_pairs(pts[:, None, :] - (x[sel] * R)[None, :, :])
            Ez = Ez + m.K / (4 * np.pi) * np.einsum("pmab,mb->pa", G, sig * M[sel] * R * np.prod(R))
        out += f"{name}: {np.linalg.norm(Ez) / np.linalg.norm(bem['E']):.1e}  "
    print(out, flush=True)
    return E


if __name__ == "__main__":
    J = int(sys.argv[1])
    prev = None
    for grp in sys.argv[2:]:
        kw = {}
        for a in grp.split():
            k, v = a.split("=")
            kw[k] = float(v) if ("." in v or "e" in v) else int(v)
        print(f"--- J={J} {grp}")
        t0 = time.time()
        E = run(J, **kw)
        if prev is not None:
            print(f"   field change from previous config: {np.linalg.norm(E - prev) / np.linalg.norm(E):.2e}")
        prev = E
