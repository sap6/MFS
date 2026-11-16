# Project notes

Last updated: 2026-11-16. All scripts are in `validation/`. Problem throughout:
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
| Cube result | done: plain MFS stalls at 21-30%; edge sources with the singular tail converge to 0.016% (J = 8) and agree with extrapolated Bempp to 0.05% (T1) |
| Cost comparison with a boundary element solver | done for sphere, cube and almond, tables and figures (T5) |
| Helmholtz material | shrinks to about a page (write-up) |

Aim (met on 2026-11-16): the cube and almond work well with the method and
beat Bempp in accuracy and time, with the plain-MFS failures kept as the
motivation. The one new idea behind both is the same: where the field is
singular (cube edges, almond tip), sources below the cut-off depth are tied
together by power laws instead of being left free ("singular tail").

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

Bempp itself (container, compilation excluded): sphere 0.37% with 9750 dofs
(411 s); cube 0.26% at h = 0.05 (8514 dofs, 323 s) and 0.19% at h = 0.04 (13302
dofs, 689 s); almond 0.67% (4212 dofs, 89 s) and 0.30% (9558 dofs, 412 s). In
every case Bempp converges towards the MFS field.

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

21. **Sphere cost table** (`python3 fig_cost_sphere.py`, k = 2 pi, against Mie,
    container, no symmetry used): MFS with 300, 600, 1200, 2400 unknowns: 0.06,
    0.16, 1.0, 5.4 s for 2.7e-3, 2.7e-6, 2.0e-9, 1.8e-13. Bempp with 1269,
    2733, 4824, 9750 dofs: 15, 46, 112, 411 s for 3.2%, 1.4%, 0.76%, 0.37%.

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

Done on 2026-11-16 (details in section 2, items 8-21): T0 clean-up, T1 cube
stall, T2 almond robustness, T3 almond general incidence, T4 almond at 7 and
9.92 GHz, T5 error-versus-time tables and figures.

T3 check of the thesis case (wave along -z, polarised along y, 1.19 GHz):
MFS backscatter 0.110432 m^2 (row spacing 0.5) and 0.110436 m^2 (0.4); Bempp
0.109478 (4212 dofs) and 0.110008 m^2 (9558 dofs), converging towards the MFS.

Open, in the order I would take them:

- [ ] **T6 Literature check** (needs arXiv access, see section 7). Is a 3D
      lightning-type or singular-source MFS for Maxwell new? Must cite and
      compare with: method of auxiliary sources (Fikioris et al.), discrete
      sources method (Eremin, Doicu, Wriedt; sources on the axis/skeleton are
      close to our almond layout), generalized multipole technique, Gopal &
      Trefethen lightning solvers, 2023 regularized MFS for PEC, arXiv
      2605.15527 (could not be identified by search), Barnett & Betcke. Read
      full texts and verify every citation by DOI.
      *Search-level pass only so far (abstracts and snippets, nothing
      verified):*
      - No 3D lightning-type or singular-source MFS for Maxwell or Helmholtz
        on a polyhedron turned up. The lightning Stokes paper says a 3D
        extension is "undoubtedly possible, but not straightforward".
      - Closest precedents for building the edge singularity into a source
        model, all 2D: Leviatan et al., IEE Proc. 1997 (fields near edges
        expanded in infinite-wedge solutions, current filaments elsewhere;
        triangular cylinder); Valtchev et al. (MFS enriched with
        corner-adapted functions, J. Comput. Appl. Math. 470); Barnett &
        Betcke 2010 (MFS plus Fourier-Bessel corner functions for polygons).
        Our tails differ: they stay inside the MFS (fixed-weight sums of
        fundamental solutions), which is what carries over to 3D edges and to
        the almond tip.
      - To read: lightning Helmholtz solver (arXiv 2310.01665; sources on the
        interior bisectors, the 2D version of our layout); Tsitsas, Zouros,
        Fikioris, Leviatan, IEEE TAP 2018 (MAS currents diverge when the
        auxiliary curve does not enclose the singularities; supports our
        almond and cube explanations); a 2026 MFS paper on 3D harmonic
        problems on the unit cube (arXiv 2603.09926); multibody MFS (arXiv
        2603.19113); "Lightning-fast method of fundamental solutions"
        (graphics, 3D Laplace).
      - Also needed for the new material: the tip-singularity exponent of an
        elliptic cone (to replace the fitted exponent set by a computed one,
        or at least to cite), and the corner exponents of the cube.
- [ ] **T7 Paper outline and critique** (separate write-up session, section 6).
- [ ] **T8 Measured almond data.** Overlay the monostatic curves of Woo et
      al. 1993 (IEEE AP Magazine) on `fig_almond_rcs`. The author has to
      supply digitised curves; nothing was compared with measurements yet.
- [ ] **T9 Cube beyond 0.02% (optional).** The balanced family is still
      converging at J = 9 (optical mismatch 4.6e-5). Nothing singular is built
      in at the vertices; if a referee asks for more digits, a vertex tail
      (same idea as the almond tip, exponents fitted from a set) is the
      candidate. Not needed for the paper as it stands.
- [ ] **T10 General incidence for the cube (optional).** The cube solver is
      hard-wired to this incidence (even in x, odd in y). The almond's
      parity-class machinery (`parity_class`, `class_part` in the almond
      scripts) shows how to do any incidence with 4 (or, with the third
      mirror, 8) class solves that share the assembly.
- [ ] **T11 Tidy the solvers (optional).** `pec_almond_skeleton.py` main
      still runs the nose-on case with row spacing 0.7 and no tip tail (fine
      for that class only); `almond_rcs.py` is the general, current solver.
      The cube keeps the name `pec_cube_lightning3d.py`.

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

- 2D square with tail (3 unknowns per corner): TE 1e-2 -> 7.8e-5, TM 1e-3 -> 2.4e-7 at J = 4.
- Cube, first 3D tail (same exponents for all components): 1.9%, 0.85%, 0.46% (J = 4, 5, 6); wrong exponents for transverse components: garbage (65%).
- Cube, tail with physical exponents in the edge frame: 0.87%, 0.35%, 0.25% vs Bempp h = 0.05.
- Cube, one knob at a time at J = 4-6 (modes, base grid, vertex rays, collocation, weights, depth spacing, taper): no single knob matters; n_base = 8 with coarse collocation: 95% error.
- Cube, balanced refinement J = 5..9: changes 0.196%, 0.085%, 0.033%, 0.016%; optical 8.2e-4 .. 4.6e-5; 0.05% from extrapolated Bempp.
- Bempp cube h = 0.1 .. 0.04: cross-section 2.91892, 2.92320, 2.92473, 2.92560; 57 .. 689 s; quadrature order irrelevant.
- Almond 1.19 GHz nose-on, collocation following the sources: layouts agree to 5e-5; Bempp 0.667%, 0.292%, 0.204% (1/N).
- Almond RCS, first run: HH off by 20-25% at oblique angles; residual of the (+,+) class 19%, half of it within 1 mm of the tip.
- Almond tip tail: (+,+) residual 19% -> 2.8%; with row spacing 0.5: 0.23%; 0.4: 0.034%. Without the tail at spacing 0.4: HH wrong by up to 1.2 dB.
- Almond RCS vs Bempp, rms dB (HH / VV): 0.023 / 0.029 (1.19 GHz), 0.044 / 0.026 (7 GHz), 0.115 / 0.088 (9.92 GHz).
- Cost tables: sphere, cube, almond (section 2, items 15, 20, 21).
- Almond, thesis broadside case: MFS 0.110432 m^2, Bempp 0.110008 m^2 (9558 dofs).

## 6. Write-up session

- Attach both `sap6/MFS` and `sap6/mfs-paper`. Read this file and the
  manuscript; do not re-run experiments.
- Target: Journal of Computational Physics; post to arXiv first. Fallback
  venue: Engineering Analysis with Boundary Elements.
- Proposed story: (1) MFS for 3D Maxwell with a Mie-validated stability and
  convergence study on the sphere; (2) why plain MFS fails on the thesis's cube
  and almond (singularities of the continued field); (3) source layouts that
  fix it: sources on the skeleton (cube: bisector planes of the edges; almond:
  mid-plane up to the focal curve), clustered towards the singular set, plus
  the singular tail at edges and tip; (4) accuracy and cost against Bempp, and
  the almond benchmark RCS. Keep the thesis shapes.
- Numbers the paper can state, each with its script: section 1 table and
  section 2 items 15, 19, 20, 21 (`fig_cost_*.py`, `fig_almond_rcs.py`,
  `almond_rcs.py`). Figures exist as PDF and PNG in `validation/`.
- What the paper must not claim, or must qualify:
  - No comparison with measured almond data was made (T8).
  - The cube and almond timings use the two mirror symmetries of the chosen
    incidence; Bempp uses none. For a general incidence the MFS needs up to
    four class solves (almond: implemented; cube: not implemented, T10).
  - The cube reference below 0.2% is the MFS itself; the independent check
    is the h^(4/3)-extrapolated Bempp field (agreement 0.05%) and the
    cross-section (2.9281 extrapolated Bempp against 2.9278-2.9284 MFS).
  - The almond tip exponents are a fitted set, not derived; the edge
    exponents of the cube are the standard wedge values, sorted by a parity
    argument that is in these notes but has not been written out as a proof.
  - "Skeleton" completeness (the field of a polyhedron is a layer potential
    on the bisector planes, by reflecting the exterior field across each
    face) is my own argument from this session, unchecked against the
    literature; treat it as a conjecture to verify or cite.
  - Novelty of the 3D tail is unverified until T6 is done.
  - Bempp was used with default settings, dense LU, no preconditioning or
    compression; say so, since a tuned BEM could be faster.
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

## 7. Working environment (read before running anything)

- **Where to run.** The cloud container (4 cores, 15 GB) is the primary copy
  and is fast enough for everything except large Bempp meshes. Its dense QR
  is faster than the laptop's. Rough times here: cube J = 6: 77 s, J = 8:
  220 s; almond RCS, all angles: 2 min per frequency; Bempp 13k dofs: 11 min.
- **The author's laptop** (`flowx13`, 16 cores, 15 GB, Linux) can be attached
  through Remote Control (folder `Projects`). It has a clone at
  `~/Projects/MFS` with a virtual environment `.venv` (numpy, scipy, numba,
  bempp-cl, gmsh). That clone is at an **old commit (efac0ae) plus hand-made
  edits** to `pec_almond_skeleton.py` and two extra scripts
  (`bem_optical.py`, `bem_almond_rcs.py`); it is only good for Bempp
  reference runs. To use it for anything else, bring it to the current
  branch first (fetch, show the author the diff, then check out). Use it for:
  Bempp with more than about 10k dofs, where it is about twice as fast.
  Run commands there as `/home/sap/Projects/MFS/.venv/bin/python3 -W ignore
  script.py ...` from `~/Projects/MFS/validation`; background jobs there do
  not notify the session, so read their output file to check.
- **Session restarts.** The cloud worker was restarted several times in one
  day. Files in the repository survive; the scratch folder and running jobs
  do not. Put anything worth keeping in `validation/` and commit often.
  Long jobs: start them in the background and write results to a file.
- **Bempp timings** must exclude its one-off compilation: do a small warm-up
  solve first (the `fig_cost_*` scripts do).
- **Network.** WebSearch works. arxiv.org (and mirrors) cannot be fetched
  from the cloud session, neither by the fetch tool nor from the container.
  For T6 the author has to allow it or supply the PDFs.
- **Scratch outputs** named `*_t_*.npy/json` and `almond_rcs_*.npz` are
  ignored by git on purpose.
- **Figures** follow one style (colour plus marker plus dash per method, so
  they print in greyscale); see `STYLE` in the `fig_*.py` scripts.
