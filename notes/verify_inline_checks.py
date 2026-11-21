"""
Standalone verification of the inline sensitivity, ablation and cross-validation
numbers stated in Sections 5, 6 and Appendices A--C of the paper that are not
stored directly in data/*.npz.

    mfs-safe-run python notes/verify_inline_checks.py [START_SECTION [END_SECTION]]
"""
import gc
import os
import sys
import numpy as np

from mfs import almond, core, cube, lines, tetrahedron

START = int(sys.argv[1]) if len(sys.argv) > 1 else 1
END = int(sys.argv[2]) if len(sys.argv) > 2 else 7
D = os.path.join(os.path.dirname(__file__), "..", "data", "")
K = 2 * np.pi
d0, p0 = core.wave()
dw, pw = core.wave(60, 25, 35)
rel = core.rel
cc = np.load(D + "cube_convergence.npz")
pts_c, ref_c = cc["pts"], cc["ref_E"]

if START <= 1 <= END:
    print("== 1. Tail unknown counts (Sections 5.2, 5.6, 6.1)", flush=True)
    c8_t = cube.Cube(8, classes=[]).unknowns[(1, -1, 1)]
    c8_nt = cube.Cube(8, tails=False, classes=[]).unknowns[(1, -1, 1)]
    print(f"Cube J=8 unknowns/class: with tails {c8_t}, without {c8_nt} (+{(c8_t - c8_nt) / c8_nt:.1%})", flush=True)
    t8_t = tetrahedron.Tetrahedron(8, classes=[]).unknowns[(1, 1, 1)]
    t8_nt = tetrahedron.Tetrahedron(8, tails=False, classes=[]).unknowns[(1, 1, 1)]
    print(f"Tetrahedron J=8 unknowns/class: with tails {t8_t}, without {t8_nt} (+{t8_t - t8_nt} cols, incl. 90 corner-ray cols)", flush=True)
    ya, ta = almond.source_rows(6, 3, s=0.5)
    print(f"Almond 6:3 (10 rows) tail unknowns: even-z={almond.tail_weights(ta, (1, 1)).shape[1]}, "
          f"odd-z={almond.tail_weights(ta, (1, -1)).shape[1]}", flush=True)

if START <= 2 <= END:
    print("\n== 2. Cube 1-depth/octave at J = 10..12 & J = 9 refinements (Section 5.4)", flush=True)
    for J in (10, 11, 12):
        gc.collect()
        Ej, _, inf = cube.solve(J, d0, p0, pts_c, K, q=0.5)
        print(f"Cube 1-depth (q=0.5) J={J}: err vs J=12 ref = {rel(Ej, ref_c):.2e}, opt = {inf['optical']:.2e}", flush=True)

    E9 = cc["dense_E9"]
    for label, kw in (("+4 extra modes", dict(extra_modes=4)), ("12x12 inner base", dict(n_base=12))):
        gc.collect()
        Er, _, _ = cube.solve(9, d0, p0, pts_c, K, **kw)
        print(f"Cube J=9 shift ({label}): {rel(Er, E9):.2e}", flush=True)

    orig_tail_e = dict(cube.TAIL_E)
    try:
        cube.TAIL_E = {"e": (2 / 3, 5 / 3 - 1.0), "m": (4 / 3, 7 / 3 - 1.0)}
        for J in (6, 8):
            gc.collect()
            Eb, _, _ = cube.solve(J, d0, p0, pts_c, K)
            base_err = rel(cc[f"dense_E{J}"], ref_c)
            bad_err = rel(Eb, ref_c)
            print(f"Cube J={J} bisector exponents -1: err {bad_err:.2e} vs {base_err:.2e} ({bad_err / base_err:.2f}x)", flush=True)
    finally:
        cube.TAIL_E = orig_tail_e

if START <= 3 <= END:
    print("\n== 3. Cube multi-wavelength mode refinements (Section 5.5)", flush=True)
    lines._TREES.clear()
    gc.collect()
    E8_k3_def = cube.solve(8, d0, p0, pts_c, 3 * K)[0]
    gc.collect()
    E8_k3_8 = cube.solve(8, d0, p0, pts_c, 3 * K, extra_modes=8)[0]
    print(f"Cube k=6pi J=8 (extra_modes 8 vs 12): shift = {rel(E8_k3_8, E8_k3_def):.2e}", flush=True)
    del E8_k3_def, E8_k3_8

    lines._TREES.clear()
    gc.collect()
    E7_k3_def = cube.solve(7, d0, p0, pts_c, 3 * K)[0]
    gc.collect()
    E7_k3_16 = cube.solve(7, d0, p0, pts_c, 3 * K, extra_modes=16)[0]
    print(f"Cube k=6pi J=7 (extra_modes 16 vs 12): shift = {rel(E7_k3_16, E7_k3_def):.2e}", flush=True)
    del E7_k3_def, E7_k3_16

    lines._TREES.clear()
    gc.collect()
    E8_k2_def = cube.solve(8, d0, p0, pts_c, 2 * K)[0]
    gc.collect()
    E8_k2_dbl = cube.solve(8, d0, p0, pts_c, 2 * K, extra_modes=8)[0]
    print(f"Cube k=4pi J=8 (extra_modes 8 vs 4): shift = {rel(E8_k2_dbl, E8_k2_def):.2e}", flush=True)
    del E8_k2_def, E8_k2_dbl

if START <= 4 <= END:
    print("\n== 4. Regular tetrahedron checks & second corner exponent (Section 5.6 & Appendix A)", flush=True)
    tc = np.load(D + "tetrahedron_convergence.npz")
    pts_t, ref_t = tc["pts"], tc["ref_E"]
    print(f"Tet Bempp n=44 (8/10) vs MFS J=10: {rel(tc['bem_E44_q8'], ref_t):.2e}, vs J=9: {rel(tc['bem_E44_q8'], tc['dense_E9']):.2e}, "
          f"n=32 vs n=44: {rel(tc['bem_E32_q8'], tc['bem_E44_q8']):.2e}", flush=True)

    orig_nu = tetrahedron.NU
    try:
        tetrahedron.NU = (0.32221,)
        gc.collect()
        Et1, _, _ = tetrahedron.solve(7, d0, p0, pts_t, K)
        print(f"Tet J=7 single electric corner exponent vs two: err {rel(Et1, ref_t):.2e} vs {rel(tc['dense_E7'], ref_t):.2e}", flush=True)
    finally:
        tetrahedron.NU = orig_nu

if START <= 5 <= END:
    print("\n== 5. Almond (+1,+1) tip-tail residuals & 10.25 GHz 35.5 deg null flank (Sections 6.1--6.2)", flush=True)
    ka = 2 * np.pi * 1.19e9 / almond.C0
    phi = np.radians(np.arange(0.0, 180.0 + 1e-9, 5.0))
    d_sw = -np.column_stack((np.cos(phi), np.sin(phi), 0 * phi))
    p_sw = np.column_stack((-np.sin(phi), np.cos(phi), 0 * phi))

    def max_rel_inc_res(s, tail_flag, s_col=0.3):
        y, tail = almond.source_rows(6, 3, s=s)
        hmax = min(0.012, 2 * np.pi / ka / 5)
        sc = s_col
        m0 = almond._count_collocation(sc, hmax)
        if 2 * m0 < 2.5 * 3 * len(y):
            sc *= np.sqrt(2 * m0 / (2.5 * 3 * len(y)))
        x, w, nrm = almond.collocation(sc, hmax)
        Tw = core.tangent_basis(nrm) * np.sqrt(w / w.sum())[:, None, None]
        b_cls = -np.einsum("mab,mwb->maw", Tw, core.class_part(x, ka, d_sw, p_sw, almond.REFL, (1, 1, 1))).reshape(-1, len(d_sw))
        inc_w = np.array([core.plane_wave(x, ka, dv, pv) for dv, pv in zip(d_sw, p_sw)])
        b_inc = np.einsum("mab,wmb->maw", Tw, inc_w).reshape(-1, len(d_sw))
        _, r_cls, _ = almond.solve_class(ka, y, (1, 1), d_sw, p_sw, tail=(tail if tail_flag else None), s_col=s_col)
        return float(np.max(r_cls * np.linalg.norm(b_cls, axis=0) / np.linalg.norm(b_inc, axis=0)))

    print(f"Almond 1.19 GHz (+1,+1) max residual rel E_inc: no tail = {max_rel_inc_res(0.5, False):.1%}, "
          f"with tail (s=0.5) = {max_rel_inc_res(0.5, True):.2%}, with tail (s=0.4) = {max_rel_inc_res(0.4, True):.3%}", flush=True)

    a10 = np.load(D + "almond_rcs_10.25GHz.npz")
    ref10 = np.loadtxt(D + "austin_almond/ref_rcs.III.A.s1.fx4.H.txt")[:, 3]
    mfs10 = 10 * np.log10(a10["sigma_HH"])
    th10 = ref10.max() - 80.0
    diff10 = np.abs(np.maximum(mfs10, th10) - np.maximum(ref10, th10))
    imax = int(np.argmax(diff10))
    print(f"Almond 10.25 GHz HH peak diff: {diff10[imax]:.2f} dB at {a10['phi_deg'][imax]:.1f} deg "
          f"(MFS {mfs10[imax]:.2f} vs ref {ref10[imax]:.2f} dBsm)", flush=True)

if START <= 6 <= END:
    print("\n== 6. Line quadrature refinement (Appendix B)", flush=True)
    orig_nodes = lines.nodes
    try:
        lines._TREES.clear()
        lines.nodes = lambda k_, H_: (20, 10)
        gc.collect()
        E6_q20 = cube.solve(6, d0, p0, pts_c, K)[0]
        print(f"Cube k=2pi J=6 line quad (20,10) vs (12,8): shift = {rel(E6_q20, cc['dense_E6']):.2e}", flush=True)

        lines._TREES.clear()
        lines.nodes = orig_nodes
        gc.collect()
        E6_k2_base = cube.solve(6, d0, p0, pts_c, 2 * K)[0]
        lines._TREES.clear()
        lines.nodes = lambda k_, H_: (32, 14)
        gc.collect()
        E6_k2_fine = cube.solve(6, d0, p0, pts_c, 2 * K)[0]
        print(f"Cube k=4pi J=6 line quad (32,14) vs (20,10): shift = {rel(E6_k2_fine, E6_k2_base):.2e}", flush=True)
    finally:
        lines.nodes = orig_nodes
        lines._TREES.clear()

if START <= 7 <= END:
    print("\n== 7. Bempp symmetry & quadrature checks (Appendix C)", flush=True)
    from mfs import bem
    bempp = bem.bempp
    pts50 = 2.0 * core.fibonacci(50)
    for reg, sing in ((4, 6), (8, 10)):
        bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = reg, sing
        Es, _ = bem.cube(8, 2.0, K, dw, pw, pts50)
        Ef, _ = bem.efie(bem.cube_grid(8, 2.0)[0], K, dw, pw, pts50)
        print(f"Bempp cube n=8 g=2.0 mirror fold vs full ({reg}/{sing}): {rel(Es, Ef):.2e}", flush=True)

    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = 6, 8
    Ets, _ = bem.tetrahedron(6, 1.8, K, dw, pw, pts50)
    Etf, _ = bem.efie(bem.tetrahedron_grid(6, 1.8)[0], K, dw, pw, pts50)
    bempp.GLOBAL_PARAMETERS.quadrature.regular, bempp.GLOBAL_PARAMETERS.quadrature.singular = 4, 4
    Eos, _ = bem.torus(8, K, d0, p0, pts50)
    Eof, _ = bem.efie(bem.torus_grid(1.0, 0.5, 8)[0], K, d0, p0, pts50)
    print(f"Bempp rotation fold vs full: tetrahedron = {rel(Ets, Etf):.2e}, torus = {rel(Eos, Eof):.2e}", flush=True)
