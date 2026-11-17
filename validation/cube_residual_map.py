"""Where is the boundary residual of the cube solution? Relative tangential
residual at points that were not used in the fit, binned by the distance to
the nearest edge and to the nearest vertex (standard wave).

    python3 cube_residual_map.py J
"""
import sys

import numpy as np

import cube_mfs as cm
import pec_cube_lightning3d as m

if __name__ == "__main__":
    J = int(sys.argv[1])
    d, p = m.D, m.P
    cube = cm.Cube(J, classes=cm.classes_of(d, p))
    coef, res = cube.solve(d, p)
    lv = int(np.ceil(np.log(cube.depths[-1] / (2 * m.H)) / np.log(0.45))) + 2
    x, w, T = m.cube_surface(lv, 4, 0.45, 0.07)                      # another graded grid: ratio 0.45, 4 nodes per panel
    sel = np.all(x > 0, axis=1)
    x, w, T = x[sel], w[sel], T[sel]
    A = cube.matrix(x, T, list(coef))
    E = sum(A[i][:, :cube.unknowns[c]] @ coef[c] for i, c in enumerate(coef)).reshape(-1, 2)
    inc = np.einsum("mab,mb->ma", T, m.incident(x))
    r2 = np.sum(np.abs(E + inc) ** 2, axis=1)
    i2 = np.sum(np.abs(inc) ** 2, axis=1)
    rho = m.edge_dist(x)
    rv = np.linalg.norm(x - m.H, axis=1)                             # distance to the vertex (+, +, +)
    print(f"J = {J}: residual of the fit {res:.2e}; at {len(x)} test points, unweighted L2: {np.sqrt(np.sum(w * r2) / np.sum(w * i2)):.2e}, "
          f"with the fit's weight rho^(2 beta): {np.sqrt(np.sum(w * rho * r2) / np.sum(w * rho * i2)):.2e}")
    edges = [0, 1e-5, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.5]
    tot = np.sum(w * rho * r2)
    print("distance to the nearest edge     points   rms residual   share of the weighted squared residual")
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (rho >= lo) & (rho < hi)
        if s.any():
            print(f"  [{lo:7.0e}, {hi:7.0e})   {s.sum():8d}   {np.sqrt(np.sum(w[s] * r2[s]) / np.sum(w[s])):10.2e}   {np.sum(w[s] * rho[s] * r2[s]) / tot:8.3f}")
    print("distance to the vertex           points   rms residual   share of the weighted squared residual")
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (rv >= lo) & (rv < hi)
        if s.any():
            print(f"  [{lo:7.0e}, {hi:7.0e})   {s.sum():8d}   {np.sqrt(np.sum(w[s] * r2[s]) / np.sum(w[s])):10.2e}   {np.sum(w[s] * rho[s] * r2[s]) / tot:8.3f}")
