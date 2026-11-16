# Project notes

Last updated: 2026-11-15. All scripts are in `validation/`. Problem throughout:
exterior scattering of a plane wave by a perfectly conducting (PEC) body,
tangential E_sc = -tangential E_inc on the surface, least-squares collocation.

## 0. Scope (from the author)

- The paper should be clearly **derived from the MSci thesis**: same topic
  (MFS for 3D Maxwell), same core method and fundamental solution, and the
  thesis's example shapes, including the **cube** and the **almond**. It should
  state that it grew out of the thesis and cite it.
- Beyond that, changes are **not limited**. There is no need to stay close to
  the thesis's text, structure, parameters, experiments or figures. Improve,
  replace or add whatever makes the paper better.
- It should not turn into an unrelated paper: the thesis's question and method
  remain the starting point, and new techniques are developments of them.

### Agreed paper shape ("Option A")

Working title: stability and convergence of the MFS for 3D Maxwell scattering.
The author liked this outline as derived from the thesis, extended, and suited
to JCP. Status of each element:

| Element | Status |
|---|---|
| Tangential-trace formulation | done in all current solvers |
| Mie-validated convergence for the PEC sphere in N and source radius | done |
| Almond at benchmark frequencies against published RCS data or Bempp | done against Bempp at 1.19, 7, 9.92 GHz, all look angles, HH and VV (T4); overlay of the measured curves (Woo et al.) left to the author |
| Cube result | honest account of where plain MFS degrades is done; edge sources with a singular tail reach about 0.15% (Bempp-limited reference) and converge with depth (T1, in progress) |
| Cost comparison with a boundary element solver | cube: table done (`fig_cost_cube.py`), MFS 2x to 10x faster than Bempp between 0.7% and 0.2% and still converging; sphere and almond tables and all plots open (T5) |
| Helmholtz material | shrinks to about a page (write-up) |

Aim: the cube and almond should work well with the method, not only be shown
as limitations. If the cube cannot be pushed much below 0.5%, present what it
reaches honestly alongside the explanation of why plain MFS fails there.

## 1. Results so far

| Shape | Method | Best result | Reference |
|---|---|---|---|
| Unit sphere, k = 2pi | dipoles on inner sphere R = 0.5 | 2e-9 (N = 400, <1 s); 2e-14 (N = 1200) | Mie |
| Unit sphere, k = 4pi | same | 1.3e-12 (N = 1200) | Mie |
| Unit cube, k = 2pi | uniform dipoles on inner cube | stalls at 25% | Bempp |
| Unit cube, k = 2pi | lightning-type edge sources, no tail | 0.5% (J = 10); 0.8% at J = 12 (stall) | Bempp h = 0.05 |
| Unit cube, k = 2pi | edge sources + singular tail, balanced refinement (section 2, items 8, 12) | 0.93%, 0.33%, 0.13%, 0.049%, 0.016% for J = 4..8 against J = 9; J = 9 is 0.05% from extrapolated Bempp; optical mismatch down to 4.6e-5 | MFS J = 9, extrapolated Bempp |
| NASA almond, 1.19 GHz, nose-on | confocal / shrunk source surface | unstable or 20%+ | - |
| NASA almond, 1.19 / 7 / 9.92 GHz, monostatic RCS, HH and VV | skeleton sources + tip tail, four parity classes | within 0.03 / 0.04 / 0.12 dB rms of Bempp; self-converged to 0.001 dB (1.19 GHz) and 0.01 dB rms (9.92 GHz) | Bempp 9558 / 13632 dofs |
| NASA almond, 1.19 GHz, nose-on | skeleton sources on mid-plane | self-converged to 5e-5, optical mismatch 1e-4; Bempp approaches it like 1/N (0.67%, 0.29%, 0.20%) | Bempp 4212, 9558, 13632 dofs |

Bempp itself: sphere 0.4% with 9750 dofs (307 s); cube at h = 0.05 (8514 dofs,
334 s in the container, 156 s on the 16-core laptop) is about 0.15-0.2% from the
extrapolated limit, not 0.12% (that was only the change from h = 0.07); almond 98 s (4212 dofs) and 275 s (9558 dofs), still
converging towards the MFS field.

## 2. What we learned (do not re-test)

1. The thesis's cube and almond results (1e-9, 1e-8) are interior problems
   whose solution is the plane wave itself. They are not scattering and hide
   every difficulty. The exterior sphere in the thesis has 2M equations for 3N
   unknowns with M = N (underdetermined).
2. Sphere: convergence is exponential in sqrt(N), not N. Small R is faster
   until the coefficient norm grows and the far field loses digits; best R is
   about 0.35 (k = 2pi) to 0.45 (k = 4pi) at N = 400.
3. Edges: uniform or shrunk-surface sources cannot converge (edge singularity,
   field ~ rho^(-1/3)). Clustering sources towards the edge works only if
   **all** of these hold (found on the 2D prototypes):
   - each clustered source has electric **and** magnetic dipoles (2D: monopole
     and dipole). Electric only stalls at 5% (cube) / 19% (2D TE);
   - successive depths differ by a factor of about 0.5. A factor 0.25 stalls;
   - sample points across the edge comfortably outnumber the unknowns per
     scale (ratio 2 or more). Ratio 1.25 gives silent garbage. Use
     `col_q=0.5 col_p=3`;
   - rows weighted by sqrt(area) * (distance to edge)^beta with beta = 0.5
     (1/3 to 0.75 similar; 0 and 1 are much worse).
4. In 3D the edge sources must be **line** dipole distributions along the edge
   (density polynomial in a log-stretched variable), not rows of points; the
   line integrals use a target-adapted quadrature (`line_kernel.py`, Numba).
5. Thin bodies (almond): the scattered field's singularities sit on the focal
   curve in the mid-plane, within 6% of the rim, so a shrunk copy of the surface
   cannot be the source surface. Point dipoles (e_x, e_y, m_z by symmetry) on
   the mid-plane, in rows clustered towards the focal curve with spacing
   proportional to the distance to the surface, converge quickly.
6. Mirror symmetry cuts unknowns and rows by 4 (solve time by ~60).
7. Did not help on the cube: bigger base grid, cell-averaged ("moment")
   testing, more modes beyond gamma = 2, extra corner rays (0.58% -> 0.51%).

8. **Singular tail (the fix for the cube stall).** The error of the edge
   sources was the cut-off of the edge singularity at the deepest free depth
   d_J (error ~ d_J^0.55), not conditioning or sampling. Below d_J, sources at
   depths d_J * 2^-i (i = 1..9, remainder lumped into the last) are tied
   together with fixed weights (d_i/d_J)^e, one unknown set per exponent e.
   The exponents must be the physical ones, in the frame of the edge (axis,
   radial = along the bisector, normal), sorted by mirror symmetry about the
   bisector plane: electric axis 2/3, 1; electric radial 5/3, 2; magnetic axis
   1, 4/3; magnetic radial 2, 7/3; no normal components. Giving every
   component 2/3 (axis) and 5/3 (transverse) doubles the error; exponents one
   too low for the transverse components give garbage. 2D prototype
   (`proto2d_square_tail.py`): TE 1e-2 -> 8e-5 and TM 1e-3 -> 2e-7 at J = 4 with
   three tail unknowns per corner; in TE only the normal dipole with 5/3
   matters, monopole tails do nothing (parity).
9. With the tail, did not matter (J = 4 or 5): more modes (gamma 3), bigger
   base grid, 4 vertex rays, col_p 4, deeper collocation, beta 0.25 (beta 0.01
   is worse), denser depths (q = 0.71), a third exponent, tail length 14
   instead of 9, tapering the tail to zero at the line ends (slightly worse).
10. Error measures that need no reference: the optical theorem
    (sigma_ext = sigma_sca, `cross_sections()`); it is not automatic for
    collocation. It is useless for Bempp, where Galerkin makes it hold to 1e-7
    on any mesh. Bempp cube cross-sections: 2.91892, 2.92320, 2.92473 for
    h = 0.1, 0.07, 0.05 (quadrature order does not matter), extrapolating to
    2.9255-2.9275; MFS J = 6: sigma_ext 2.92686, sigma_sca 2.92576.
11. The Bempp mesh sequence converges towards the MFS field: MFS (J = 6) is
    0.63%, 0.36%, 0.25% from Bempp h = 0.1, 0.07, 0.05.

12. **Refine everything together.** With the tail, depth J alone converges
    only to about 0.1%: other resolution limits then show (variants at J = 6
    differ by 0.03-0.12%). The balanced family (`balanced()` in
    `fig_cost_cube.py`: n0 = 2J - 4 modes offset, base grid 4/6/8 with
    interior panels 0.25/0.15/0.1) converges cleanly: J = 5, 6, 7, 8 change by
    0.196%, 0.085%, 0.033%; optical mismatch 8.2e-4, 3.3e-4, 1.5e-4, 7.9e-5;
    distance to extrapolated Bempp 0.33%, 0.14%, 0.07%, 0.05%.
    A fine base grid needs fine interior collocation: n_base = 8 with panels
    of 0.25 gives garbage (95%) without any warning except the optical
    mismatch (0.5).
13. Timing on the same machine (container, 4 cores): Bempp h = 0.1, 0.07,
    0.05, 0.04 takes 57, 166, 323, 689 s and is 0.68%, 0.40%, 0.25%, 0.18% from
    the converged MFS field. The MFS uses the two mirror symmetries of this
    incidence; Bempp does not (say so in the paper).
14. LAPACK: plain QR is twice as fast as gelsy and equally accurate here
    (rank was always full). The laptop's 16 cores speed up Bempp and the
    Numba assembly about 2x but its dense QR was slower than the container's.

15. **Cube cost table** (`python3 fig_cost_cube.py`, container, 4 cores, error
    against MFS J = 9): MFS + tail J = 3..8: 16.5, 25, 42, 77, 154, 220 s for
    2.7%, 0.93%, 0.33%, 0.13%, 0.049%, 0.016%. Bempp h = 0.1, 0.07, 0.05, 0.04:
    57, 166, 323, 689 s for 0.68%, 0.39%, 0.26%, 0.19%. Edge sources without
    tail, J = 4, 6, 8: 9, 37, 125 s for 5.4%, 2.5%, 1.1%. Plain MFS, 72 to 450
    unknowns: 30% to 21%, no convergence.

16. **Almond, all look angles** (`almond_rcs.py`). Any plane wave splits into
    four parity classes under y -> -y and z -> -z; each is solved with a
    quarter of the sources, and one QR per class serves every angle. HH uses
    the two classes even in z (dipoles e_x, e_y, m_z on the mid-plane), VV the
    two odd in z (e_z, m_x, m_y).
17. **Tip tail (almond).** The class even in both mirrors carries current
    along the axis into the tip (an elliptic cone, half-angles 34 and 12.7
    degrees), where the field is strongly singular. The nose-on case never
    excited it. Free sources cut off 1 mm from the tip leave 19% residual in
    that class, half of it within 1 mm of the tip, and HH errors up to 1.2 dB;
    more sources do not help. Fix (`tip_tail`, `tail_matrix`): continue every
    source row to 1e-4 mm from the tip and tie the continued dipoles together
    with weights (r / r_cut)^e, e in {2.2, 2.4, 2.7, 3.0} (dipole exponent =
    field exponent nu + 2; the true nu of the elliptic cone is not known, the
    fit combines the set). Electric dipoles only; adding magnetic ones, more
    or lower exponents, a longer tail or collocation closer to the tip changes
    nothing. Close exponents make the coefficient norm huge (1e7 to 1e12)
    without hurting the result; the set of four keeps it at 1e3.
18. **Row spacing.** For that class the spacing along a source row must be
    0.5 x (distance to the surface) or less: residual 3.1%, 0.23%, 0.034% for
    0.7, 0.5, 0.4. More rows help far less. The other three classes were
    already fine at 0.7.
19. Almond RCS at 1.19 GHz, 19 angles: spacing 0.5 against 0.4 differs by
    0.0011 dB (HH) and 0.0002 dB (VV); against Bempp 9558 dofs 0.023 and
    0.029 dB rms. Bempp references for 1.19, 7, 9.92 GHz are stored in
    `bem_almond_rcs_ref.json` (from `bem_almond_rcs.py`, laptop).

20. **Almond cost table** (`python3 fig_cost_almond.py`, 1.19 GHz nose-on,
    container, 4 cores, error against the finest MFS layout): MFS with 795,
    1377, 2103, 2937 unknowns: 2.4, 4.1, 7.9, 13.5 s for 6.7%, 0.40%, 0.029%,
    0.0056%. Bempp with 2106, 4212, 9558 dofs: 28, 89, 412 s for 1.3%, 0.67%,
    0.30% (compilation excluded by a warm-up solve).

## 3. Current best commands

```
python3 pec_sphere_convergence.py && python3 plot_pec_sphere.py
python3 fig_cost_cube.py ; python3 fig_cost_almond.py ; python3 fig_cost_sphere.py ; python3 fig_almond_rcs.py   # tables and figures from saved data (add `run` to recompute)
python3 bem_reference.py cube                      # writes bem_cube_h*.npz
python3 pec_cube_lightning3d.py 4,5,6 tag=NAME kinds=em q=0.5 gamma=2.0 col_p=3 col_q=0.5 beta=0.5 sym=1 tail=9
python3 diag_cube.py 5 "tail_n=9"                 # residual map by zone, optical theorem
python3 bem_optical.py 0.1,0.07,0.05              # Bempp cross-sections and far field
python3 pec_almond_skeleton.py 1.19 --bem
python3 almond_rcs.py 1.19 6:3 10 --s=0.5         # monostatic RCS, HH and VV, compares with stored Bempp
python3 proto2d_square.py TE ; python3 proto2d_ellipse.py     # 2D prototypes, seconds
python3 proto2d_square_tail.py TE                 # 2D singular tail
```

Cube cost with tail (container, 4 cores): J = 4: 1470 unknowns, 37 + 7 s;
J = 5: 1800, 49 + 13 s; J = 6: 2124, 62 + 22 s. Without tail at J = 10: 3126
unknowns, 35k rows, ~90 + 85 s for 0.58%.
Almond: 979 sources -> 21 s solve; 1597 -> 50 s.

Files: `mie.py` (reference), `bem_reference.py` (Bempp EFIE),
`pec_cube_lightning3d.py` + `line_kernel.py` (current cube solver),
`diag_cube.py` (residual map), `bem_optical.py` (Bempp cross-sections, timed),
`pec_almond_skeleton.py` (current almond solver), `pec_almond.py` (geometry, mesh,
and the confocal solver kept as the record of the failed layout).

## 4. Task list (one per session)

- [x] **T0 Clean-up.** Done 2026-11-16: one-off cube result files and the
      superseded solvers `pec_cube_edge.py`, `pec_cube_lightning.py`,
      `check_cube_dipole.py` removed; their helpers now live in
      `pec_cube_lightning3d.py`. Kept on purpose: `pec_almond.py` (geometry,
      mesh, and the confocal solver that documents why a shrunk source surface
      fails) and `results_pec_cube.json` (uniform sources, 25%). Plain MFS on
      the cube is reproduced by `pec_cube_lightning3d.py -1 n_base=8`.
- [x] **T1 Cube stall.** Solved 2026-11-16 by the singular tail (section 2,
      items 8-11). Self-convergence 0.21%, 0.083%, 0.050% for J = 6, 7, 8;
      optical-theorem mismatch 3e-4; cross-section agrees with h^(4/3)-
      extrapolated Bempp to 2e-4. Open: the optical mismatch stops falling at
      3e-4 from J = 7 (next error source not identified; vertices are the
      candidate, as nothing singular is built in there).
- [x] **T2 Almond robustness.** Done 2026-11-16. `solve()` now refines the
      collocation until rows/unknowns >= 2.5. Layouts with 979, 1597, 2246
      sources: field changes 5.5e-5 and 4.0e-5, optical-theorem mismatch
      0.8e-4 to 1.5e-4, cross-section 5.7452e-3 to 5.7455e-3 m^2. Bempp with
      4212, 9558, 13632 dofs differs from the MFS by 0.667%, 0.292%, 0.204%:
      proportional to 1/N with the MFS as the limit. The boundary residual
      stays at 2-4e-3 while the field converges (it does not radiate; where it
      sits was not checked).
- [x] **T3 Almond, general incidence.** Covered by `almond_rcs.py`: any
      plane wave is handled through the four parity classes. Only in-plane
      look angles were run; the thesis case (d = -z, p = y) goes through the
      same `solve_class` but was not run or checked.
- [x] **T4 Almond at benchmark frequencies.** Done 2026-11-16 against Bempp;
      **not** against the measured data of Woo et al. 1993, which could not be
      fetched from the cloud session (the author should overlay digitised
      curves). Monostatic RCS, HH and VV, azimuth 0-180 degrees, layout 6:3
      with row spacing 0.5: about 4200 unknowns per class and 2 minutes per
      frequency at 1.19, 7 and 9.92 GHz alike (the source count is set by the
      geometry, not the wavelength, up to 8.35 wavelengths). MFS against Bempp
      (9558 dofs at 1.19 GHz, 13632 at 7 and 9.92): rms 0.023 / 0.029 dB (HH /
      VV) at 1.19, 0.044 / 0.026 dB at 7, 0.115 / 0.088 dB at 9.92 GHz (max
      0.46 dB at a null). MFS self-convergence at 9.92 GHz: row spacing 0.4,
      collocation 0.22, or layout 8:5 change HH by at most 0.036, 0.004, 0.034
      dB and VV by at most 0.0024 dB. Figure: `fig_almond_rcs.py`.
- [ ] **T5 Error-versus-time plots,** MFS against Bempp, for sphere, cube and
      almond. One script per figure.
- [ ] **T6 Literature check.** Is a 3D lightning-type method for Maxwell new?
      Must cite and compare with: method of auxiliary sources (Fikioris et
      al.), discrete sources method (Eremin, Doicu, Wriedt; sources on the
      axis/skeleton are close to our almond layout), generalized multipole
      technique, Gopal & Trefethen lightning solvers, 2023 regularized MFS for
      PEC, arXiv 2605.15527, Barnett & Betcke 2008. Read full texts.
      *Search-level pass, 2026-11-16 (abstracts and snippets only; arxiv.org is
      blocked from the cloud session, so no full text was read and nothing
      here is verified):*
      - No 3D lightning-type or singular-source MFS for Maxwell or Helmholtz
        on a polyhedron turned up. The lightning Stokes paper says a 3D
        extension is "undoubtedly possible, but not straightforward".
      - Closest precedents for building the edge singularity into a source
        model, all 2D: Leviatan et al., IEE Proc. 1997 (fields near edges
        expanded in infinite-wedge solutions, current filaments elsewhere;
        triangular cylinder); Valtchev et al. (MFS enriched with
        corner-adapted functions, J. Comput. Appl. Math. 470); Barnett &
        Betcke 2010 (MFS plus Fourier-Bessel corner functions for polygons).
        Our tail differs: it stays inside the MFS (fixed-weight sums of
        fundamental solutions), which is what carries over to 3D edges.
      - To read: lightning Helmholtz solver (arXiv 2310.01665; sources on the
        interior bisectors, the 2D version of our layout); Tsitsas, Zouros,
        Fikioris, Leviatan, IEEE TAP 2018 (MAS currents diverge when the
        auxiliary curve does not enclose the singularities; supports our
        almond and cube explanations); a 2026 MFS paper on 3D harmonic
        problems on the unit cube (arXiv 2603.09926); multibody MFS (arXiv
        2603.19113); "Lightning-fast method of fundamental solutions"
        (graphics, 3D Laplace); arXiv 2605.15527 (not yet seen).
- [ ] **T7 Paper outline** (separate write-up session, see section 6).

## 5. Experiment log (newest last; one line each)

- Sphere sweep, k = 2pi and 4pi, N to 1200, R 0.1-0.95: see table; sqrt(N) rate.
- Cube, uniform sources n = 4..11: residual 41% -> 32%, far field 25% off.
- Cube, point clusters at edges (electric only): 32% -> 13% residual, not converged.
- Cube, electric line sources, q = 0.25: 25% -> 5%, stalls at J = 3.
- Cube, three lines per depth (electric): 2.9% at J = 2, stalls at J = 3.
- 2D square: monopoles only stall (TM 0.5%, TE 19%); with dipoles TM 1e-8, TE 5e-6.
- Cube, e+m, q = 0.5, col_q = 1/3: 5.4%, 2.6%, 1.3% (J = 4, 6, 8), then worse.
- Cube, e+m, q = 0.5, col_q = 0.5, gamma = 2: 2.5%, 1.1%, 0.58%, 0.80% (J = 6..12).
- Cube, q = 0.7 with col_q = 1/3: blows up (128%), under-sampled.
- Cube, vrays = 4 at J = 10: 0.51%.
- Almond, confocal sources, first layout: unstable (known-solution check failed).
- Almond, confocal, rim-resolved: stable, 22% residual with 1243 sources.
- 2D ellipse: confocal monopoles and skeleton mono+dipole both reach 1e-12.
- Almond, skeleton sources, nose-on, 1.19 GHz: residual 3.7%, 0.31%, 0.24%; field
  change 5e-5; 0.30% from Bempp (9558 dofs); 2246 sources under-sampled.

## 6. Write-up session

- Attach both `sap6/MFS` and `sap6/mfs-paper`. Read this file and the
  manuscript; do not re-run experiments.
- Target: Journal of Computational Physics; post to arXiv first. Fallback
  venue: Engineering Analysis with Boundary Elements.
- Proposed story: (1) MFS for 3D Maxwell with a Mie-validated stability and
  convergence study on the sphere; (2) why plain MFS fails on the thesis's cube
  and almond (singularities of the continued field); (3) source layouts that
  fix it (edge line sources; skeleton sources). Keep the thesis shapes.
- Known problems in the current manuscript: thesis wording ("this project",
  "chapter"); 2D Helmholtz material too long (cut to about a page); interior
  problems labelled as scattering; no Mie comparison; claim of exponential
  convergence in N; GMRES / block-structure paragraph does not match the code;
  conclusion over-claims; 12 references, none on MAS / GMT / DSM.
- The author writes the prose from their own bullet points. Claude gives
  referee-style critique, checks every number against script output, verifies
  each citation by DOI, and formats. AI use will be declared (Elsevier policy:
  writing help in a declaration, research use in the methods).
- Overleaf routine: push Overleaf -> GitHub before a session, work on a
  branch, merge, then pull GitHub -> Overleaf. Ask before pushing.
