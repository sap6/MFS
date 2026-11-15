"""Where does the cube boundary residual sit, on the fitting grid and on an
independent test grid? Usage: python diag_cube.py J [col_p] [col_q]"""
import sys
import numpy as np
import pec_cube_lightning3d as m

J = int(sys.argv[1]); col_p = int(sys.argv[2]) if len(sys.argv) > 2 else 3
col_q = float(sys.argv[3]) if len(sys.argv) > 3 else 1 / 3
src = m.SymSources(J, q=0.5, gamma=1.0, kinds="em")
d_min = src.depths[-1]
lv = int(np.ceil(np.log(d_min / (2 * m.H)) / np.log(col_q))) + 1
c, info = m.solve(src, lv, col_p, col_q, 0.5)
bem = m.load_bem()
E = (src.matrix(bem["pts"]) @ c).reshape(-1, 3)
print(f"J={J} col_p={col_p} col_q={col_q:.2f} unk={info['unknowns']} rows={info['rows']} "
      f"vs_bem={np.linalg.norm(E - bem['E']) / np.linalg.norm(bem['E']):.3e}", flush=True)
edges = [0, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.5]
for name, (lvl, p, q) in (("fit grid ", (lv, col_p, col_q)), ("test grid", (int(lv * np.log(col_q) / np.log(0.45)) + 1, 4, 0.45))):
    x, w, T = m.cube_surface(lvl, p, q)
    sel = m.quarter(x); x, w, T = x[sel], w[sel], T[sel]
    rhs = -np.einsum("mab,mb->ma", T, m.incident(x))
    res = (src.matrix(x, T) @ c).reshape(-1, 2) - rhs
    m2 = np.sum(np.abs(res) ** 2, axis=1)
    rho = m.edge_dist(x)
    vd = np.linalg.norm(m.H - np.abs(x), axis=1)
    ww = w * rho
    tot = np.sum(ww * np.sum(np.abs(rhs) ** 2, axis=1))
    print(f"{name}: {len(x)} pts, weighted (rho^1) rel residual {np.sqrt(np.sum(ww * m2) / tot):.3e}")
    line = "   rms by edge distance: "
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (rho >= lo) & (rho < hi)
        if s.any():
            line += f"[{lo:.0e},{hi:.0e}) {np.sqrt(np.sum(w[s] * m2[s]) / np.sum(w[s])):.1e}  "
    print(line)
    line = "   rms by vertex distance (rho<0.01): "
    near = rho < 0.01
    for lo, hi in zip([0, 1e-4, 1e-3, 1e-2, 0.1], [1e-4, 1e-3, 1e-2, 0.1, 1]):
        s = near & (vd >= lo) & (vd < hi)
        if s.any():
            line += f"[{lo:.0e},{hi:.0e}) {np.sqrt(np.sum(w[s] * m2[s]) / np.sum(w[s])):.1e}  "
    print(line, flush=True)
