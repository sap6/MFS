"""Where does the cube boundary residual sit? Bins by distance to the nearest
edge and to the nearest vertex. Usage: python diag_cube.py J [col_p]"""
import sys
import numpy as np
import pec_cube_lightning3d as m

J = int(sys.argv[1]); col_p = int(sys.argv[2]) if len(sys.argv) > 2 else 2
src = m.Sources(J)
d_min = src.depths[-1]
lv = int(np.ceil(np.log(d_min / (2 * m.H)) / np.log(1 / 3))) + 1
c, info = m.solve(src, lv, col_p, 1 / 3, 1 / 3)
bem = m.load_bem()
E = (src.matrix(bem["pts"]) @ c).reshape(-1, 3)
print(f"J={J} col_p={col_p} unk={info['unknowns']} rows={info['rows']} "
      f"vs_bem={np.linalg.norm(E - bem['E']) / np.linalg.norm(bem['E']):.3e}", flush=True)
x, w, T = m.cube_surface(lv + 1, 3, 0.4)          # independent test grid
rhs = -np.einsum("mab,mb->ma", T, m.incident(x))
res = (src.matrix(x, T) @ c).reshape(-1, 2) - rhs
m2 = np.sum(np.abs(res) ** 2, axis=1)
rho = m.edge_dist(x)
vd = np.linalg.norm(m.H - np.abs(x), axis=1)
tot = np.sum(w * np.sum(np.abs(rhs) ** 2, axis=1))
print("unweighted rel L2 residual:", np.sqrt(np.sum(w * m2) / tot))
edges = [0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.5]
print("by edge distance: [lo, hi)  rms residual   share of L2^2   share with rho^(2/3) weight")
ww = w * rho ** (2 / 3)
for lo, hi in zip(edges[:-1], edges[1:]):
    s = (rho >= lo) & (rho < hi)
    if s.any():
        print(f"  [{lo:.0e},{hi:.0e})  {np.sqrt(np.sum(w[s] * m2[s]) / np.sum(w[s])):.2e}   "
              f"{np.sum(w[s] * m2[s]) / np.sum(w * m2):.2f}   {np.sum(ww[s] * m2[s]) / np.sum(ww * m2):.2f}")
print("by vertex distance (points within 0.03 of an edge only):")
near = rho < 0.03
for lo, hi in zip([0, 0.01, 0.03, 0.1, 0.3], [0.01, 0.03, 0.1, 0.3, 1]):
    s = near & (vd >= lo) & (vd < hi)
    if s.any():
        print(f"  [{lo},{hi})  rms {np.sqrt(np.sum(w[s] * m2[s]) / np.sum(w[s])):.2e}")
