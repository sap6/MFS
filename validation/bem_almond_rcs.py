"""Bempp reference (EFIE, RWG, dense LU) for the monostatic RCS of the NASA
almond in its own plane, HH and VV: one factorisation, all look angles.
Conventions as in almond_rcs.py (azimuth 0 = nose-on).

    python3 bem_almond_rcs.py FREQ_GHz n_t n_s [step_deg]
Writes bem_almond_rcs_FREQGHz_NDOF.npz and prints sigma in m^2.
"""
import sys
import time

import bempp_cl.api as bempp
import numpy as np
from scipy.linalg import lu_factor, lu_solve

from pec_almond import mesh

C0 = 299792458.0


@bempp.callable(complex=True, parameterized=True)
def trace(x, n, domain_index, result, par):
    ph = np.exp(1j * par[6] * (x[0] * par[0] + x[1] * par[1] + x[2] * par[2]))
    result[0] = ph * (par[4] * n[2] - par[5] * n[1])
    result[1] = ph * (par[5] * n[0] - par[3] * n[2])
    result[2] = ph * (par[3] * n[1] - par[4] * n[0])


if __name__ == "__main__":
    freq, n_t, n_s = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    step = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    k = 2 * np.pi * freq * 1e9 / C0
    t0 = time.time()
    v, e = mesh(n_t, n_s)
    grid = bempp.Grid(v.T, e.T.astype("uint32"))
    div = bempp.function_space(grid, "RWG", 0)
    curl = bempp.function_space(grid, "SNC", 0)
    A = bempp.as_matrix(bempp.operators.boundary.maxwell.electric_field(div, div, curl, k).weak_form())
    lu = lu_factor(A, overwrite_a=True)
    phi = np.deg2rad(np.arange(0.0, 180.0 + 1e-9, step))
    d = -np.column_stack((np.cos(phi), np.sin(phi), np.zeros_like(phi)))
    rfar = 1.0e5
    pot = bempp.operators.potential.maxwell.electric_field(div, (-rfar * d).T, k)
    out = {"phi_deg": np.rad2deg(phi), "freq": freq, "ndof": div.global_dof_count}
    for name in ("HH", "VV"):
        sig = []
        for i, ph in enumerate(phi):
            p = np.array([-np.sin(ph), np.cos(ph), 0.0]) if name == "HH" else np.array([0.0, 0.0, 1.0])
            rhs = bempp.GridFunction(div, fun=trace, dual_space=curl,
                                     function_parameters=np.concatenate((d[i], p, [k]))).projections(curl)
            lam = bempp.GridFunction(div, coefficients=lu_solve(lu, rhs))
            F = -(pot * lam)[:, i] * rfar * np.exp(-1j * k * rfar)
            sig.append(4 * np.pi * abs(p @ F) ** 2)
        out["sigma_" + name] = np.array(sig)
    np.savez(f"bem_almond_rcs_{freq}GHz_{out['ndof']}.npz", **out)
    print(f"almond {freq} GHz, Bempp dofs={out['ndof']}, {time.time() - t0:.0f}s")
    print("azimuth  sigma_HH[m^2]  sigma_VV[m^2]")
    for a_, hh, vv in zip(out["phi_deg"], out["sigma_HH"], out["sigma_VV"]):
        print(f"{a_:6.1f} {hh:13.5e} {vv:13.5e}")
