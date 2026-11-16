# Project notes

Last updated: 2026-11-16 (second session of that day). All scripts are in
`validation/`. Companion files: `DERIVATIONS.md` (proofs and exponent tables),
`LITERATURE.md` (T6), `PAPER_OUTLINE.md` (T7). Problem throughout:
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
| Almond at benchmark frequencies against published RCS data or Bempp | done: against Bempp at 1.19, 7, 9.92 GHz (T4) and against the public Austin benchmark (measurements and reference simulation) at 3.5, 5.125, 7, 10.25 GHz, all look angles, HH and VV (T8) |
| Cube result | done: plain MFS stalls at 21-30%; edge sources with edge and vertex tails converge to 0.011% (J = 8) and agree with graded Bempp to 0.022%, with its extrapolation to 0.016% (T1, T9); any incident wave (T10) |
| Cost comparison with a boundary element solver | done for sphere, cube and almond, tables and figures (T5); cube also against Bempp on edge-graded meshes |
| Helmholtz material | shrinks to about a page (write-up) |

Aim (met on 2026-11-16): the cube and almond work well with the method, with
the plain-MFS failures kept as the motivation. The one new idea behind both
is the same: where the field is singular (cube edges and vertices, almond
tip), sources below the cut-off depth are tied together by power laws instead
of being left free ("singular tail"); the exponents are computed, not fitted.
Against Bempp: far ahead on the sphere and the almond; on the cube far ahead
of quasi-uniform meshes but only level with edge-graded meshes around 0.15%
(ahead below 0.1%), see section 2, item 27. Do not write "beats Bempp" for
the cube without that qualification.

## 1. Results so far

| Shape | Method | Best result | Reference |
|---|---|---|---|
| Unit sphere, k = 2pi | dipoles on inner sphere R = 0.5 | 2e-9 (N = 400, <1 s); 2e-14 (N = 1200) | Mie |
| Unit sphere, k = 4pi | same | 1.3e-12 (N = 1200) | Mie |
| Unit cube, k = 2pi | uniform dipoles on inner cube | stalls at 25% | Bempp |
| Unit cube, k = 2pi | lightning-type edge sources, no tail | 0.5% (J = 10); 0.8% at J = 12 (stall) | Bempp h = 0.05 |
| Unit cube, k = 2pi | edge sources + edge and vertex tails, balanced refinement (section 2, items 8, 12, 22, 25) | 1.27%, 0.52%, 0.20%, 0.089%, 0.033%, 0.011% for J = 3..8 against J = 9; optical mismatch down to 4.0e-5 | MFS J = 9; graded Bempp (10368 dofs) is 0.022% from J = 9, its extrapolation 0.016% |
| Unit cube, oblique wave (from theta 60, phi 25, polarisation angle 35 degrees) | the same, four parity classes (`cube_general.py`) | see item 26 | graded Bempp, 10368 dofs |
| NASA almond, 1.19 GHz, nose-on | confocal / shrunk source surface | unstable or 20%+ | - |
| NASA almond, 1.19 / 7 / 9.92 GHz, monostatic RCS, HH and VV | skeleton sources + tip tail, four parity classes | within 0.03 / 0.04 / 0.12 dB rms of Bempp; self-converged to 0.001 dB (1.19 GHz) and 0.01 dB rms (9.92 GHz) | Bempp 9558 / 13632 dofs |
| NASA almond, 1.19 GHz, nose-on | skeleton sources on mid-plane | self-converged to 5e-5, optical mismatch 1e-4; Bempp approaches it like 1/N (0.67%, 0.29%, 0.20%) | Bempp 4212, 9558, 13632 dofs |
| NASA almond, 3.5 / 5.125 / 7 / 10.25 GHz, monostatic RCS, HH and VV | the same solver as the RCS row above, tip exponents from the cone | mean difference from the Austin reference simulation 0.04-0.14 dB (HH), 0.01-0.02 dB (VV); from the measurements 0.6-1.8 dB, the same as the reference simulation's own distance from them | Austin RCS benchmark III-A (Kelley et al. 2020) |

Bempp itself (container, compilation excluded): sphere 0.37% with 9750 dofs
(411 s); cube, quasi-uniform: 0.26% at h = 0.05 (8514 dofs, 393 s) and 0.19% at
h = 0.04 (13302 dofs, 832 s); cube, graded to the edges: 0.105% with 4608 dofs
(161 s); almond 0.67% (4212 dofs, 89 s) and 0.30% (9558 dofs, 412 s). In every
case Bempp converges towards the MFS field. (Cube times are from the second
session's container, about 20% slower than the first; sphere and almond times
are from the first.)

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

8. **Singular tail (the fix for the cube stall).** *(Exponent list superseded
   by item 22: one derived exponent per component.)* The error of the edge
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

15. **Cube cost table, first version** *(superseded by item 27)*
    (`python3 fig_cost_cube.py`, container, 4 cores, error
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
17. **Tip tail (almond).** *(The fitted exponent set is superseded by item 24.)* The class even in both mirrors carries current
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


### Second session, 2026-11-16 (items 22-32)

22. **Edge tail exponents are derived** (`DERIVATIONS.md` B). The tail stands
    for the surface current that the continued field leaves on the bisector
    half-plane of the edge. Parity under the mirror in that plane sorts the
    Meixner terms: even terms give only electric current, odd terms only
    magnetic, nothing normal to the plane; integer-power terms give no
    current at all. One exponent per component: electric axis 2/3, electric
    radial 5/3, magnetic axis 4/3, magnetic radial 7/3 (now the defaults of
    `Sources`). Checked with `check_tail_exponents.py`: same error as the old
    sets with integer second exponents (J = 5: 0.314% against 0.327%; J = 7:
    0.049% both), 8-10% fewer unknowns; the derived second terms (8/3, 11/3)
    change nothing; the J = 9 field moved by 1e-5.
23. **Skeleton theorem** (`DERIVATIONS.md` A). For a convex PEC polyhedron
    the scattered field is exactly a layer of electric and magnetic current
    on the interior medial surface (image principle across each face, patched
    over the nearest-face cells). It explains why skeleton sources are
    complete, why a shrunk surface fails, and why both dipole kinds are
    needed. Not found in the literature; convex bodies only.
24. **Tip exponents are computed** (`cone_exponents.py`, `DERIVATIONS.md` C).
    Method of particular solutions on the sphere with poles on the focal
    segment of the cone's cross-section; exact to 8 digits for a circular
    cone, stable to 5 digits from 12 to 36 poles. Almond tip (half-angles
    34.09 and 12.71 degrees): the only exponent below 1 is nu = 0.29959
    (E ~ r^-0.70), in the class even under both mirrors. One electric
    exponent e = nu + 2 = 2.2996 replaces the fitted set {2.2, 2.4, 2.7, 3.0}:
    same residual (0.24% at spacing 0.5), HH closer to the fine layout
    (0.0001 against 0.0003 dB rms). `almond_rcs.rcs` now takes the exponents
    per class from the cone (two electric, one magnetic); more exponents
    change nothing.
25. **Cube vertices** (`python3 cone_exponents.py cube`, `DERIVATIONS.md` D).
    Exponents outside a cube corner: nu = 0.45417 (the known value), 1.23087;
    tau = 0.8400, 1.20635. A vertex tail on the body diagonal (`VertexTail`,
    e = nu + 2 and tau + 2, 12 unknowns) cuts the error by a third to 40% at
    every depth: J = 5: 0.31% -> 0.20%, J = 6: 0.133% -> 0.086%, J = 7: 0.049%
    -> 0.029% (against J = 9), and 0.059% -> 0.042% against graded Bempp. A
    second exponent changes nothing. It is now part of the balanced family.
26. **Cube, any incident wave** (`cube_general.py`, T10). Four parity classes
    under x -> -x, y -> -y; the class matrices are signed sums of the same
    four reflected blocks, so the assembly is shared and each class has its
    own QR. It reproduces the symmetric solver to 2e-15. Oblique wave (from
    theta 60, phi 25, polarisation angle 35 degrees), J = 5, 6, 7:
    86, 176, 370 s (assembly 47, 84, 169 s; four factorisations 25, 70, 171 s)
    for 0.27%, 0.112%, 0.055% against graded Bempp (10368 dofs, quadrature
    6/8, itself 0.016% from its 7200-dof mesh); optical mismatch 8.8e-4,
    4.0e-4, 2.1e-4. The symmetric case costs 53, 95, 199 s, so a general wave
    costs 1.6-1.9 times as much. Graded Bempp (no symmetry needed) gives
    0.15% in 76 s and 0.105% in 161 s for any wave: for a general incidence
    it is faster than the MFS at 0.1-0.4% and slower below about 0.06%.
27. **Cube cost table, final** (`python3 fig_cost_cube.py`, second container,
    4 cores, everything timed in one session, error against MFS J = 9 with
    5082 unknowns and 419 s):
    - MFS, edge and vertex tails, J = 3..8: 22, 32, 53, 95, 199, 284 s for
      1.27%, 0.52%, 0.20%, 0.089%, 0.033%, 0.011%; optical mismatch 3.0e-3 to
      6.4e-5.
    - Edge sources without tails, J = 4, 6, 8: 14, 54, 180 s for 5.4%, 2.5%,
      1.1%. Plain MFS, 72 to 450 unknowns: 30% to 21%.
    - Bempp, quasi-uniform, h = 0.1, 0.07, 0.05, 0.04 (2214 to 13302 dofs):
      63, 174, 393, 832 s for 0.68%, 0.39%, 0.26%, 0.19%.
    - Bempp, graded to the edges (`bem_cube_graded.py`, grading exponent 2,
      default quadrature), 1152, 2592, 4608 dofs: 32, 76, 161 s for 0.44%,
      0.15%, 0.105%. With quadrature orders 6/8, 4608 dofs: 555 s for 0.060%.
    - Reading: against quasi-uniform Bempp the MFS is 10-15 times faster at
      0.2%; against graded Bempp it is level at 0.15-0.2%, 1.7 times faster
      at 0.1%, about 4 times at 0.06%, and graded Bempp was not taken below
      0.02% (25 min of assembly on the laptop for 0.022%).
    - Bempp's time is assembly: 158 of 161 s at 4608 dofs, 1521 of 1540 s at
      10368 dofs (laptop). A faster linear solver would not change the
      comparison; matrix compression or higher-order elements might, and
      were not tried.
28. **Graded Bempp meshes** (`bem_cube_graded.py`). Tensor grid per face with
    nodes 0.5 (1 - (1 - |s|)^g), each rectangle cut in two triangles
    (anisotropic near the edges). g = 2 is best; g = 1.5 similar at 4608 dofs
    (0.113%); g = 2.5 worse (0.20% against 0.13% at 2592 dofs); g = 3 fails
    (4.8%) with default quadrature. The default quadrature (4 regular, 6
    singular) limits the graded meshes at 0.1%: 4608 dofs give 0.104%, 0.060%,
    0.057% with orders 4/6, 6/8, 8/10, and the gap grows with the mesh
    (7200 dofs: 0.038% against 0.032%). Use 8/10 for reference runs.
29. **Independent reference for the cube.** Graded Bempp with quadrature
    8/10 and 4608, 7200, 10368 dofs is 0.058%, 0.032%, 0.022% from MFS J = 9;
    extrapolated from these three (fitted rate n^-2.9) 0.016%; cross-section
    2.92710, 2.92745, 2.92762, extrapolated 2.92785. Extrapolated
    quasi-uniform Bempp stays 0.05% away (unchanged). MFS self-convergence:
    0.011% at J = 8. So the cube field is confirmed by a second method to
    about 0.02%; below that the reference is the MFS itself.
30. **Almond against the Austin RCS benchmark** (`fig_almond_benchmark.py`,
    data in `validation/austin_almond/`, T8). Same almond (9.936 in), same
    angle and polarisation conventions. Mean |difference| in dB over 361
    angles, MFS against the suite's reference simulation (ARCHIE-AIM, finest
    mesh), HH / VV: 3.5 GHz 0.043 / 0.011; 5.125 GHz 0.055 / 0.013; 7 GHz
    0.072 / 0.016; 10.25 GHz 0.140 / 0.024. MFS against the two measurements:
    0.64-1.83 dB, within 0.07 dB of the reference simulation's own distance
    from them in every case. Largest pointwise differences from the reference
    simulation are in deep HH nulls: 1.5 dB (3.5 GHz, 130 deg), 0.9, 0.7 and
    4.3 dB (10.25 GHz, 35.5 deg, a null at -66.7 dBsm). There the MFS is
    converged: layouts 6:3 and 8:5 differ by 0.0006 dB on average and 0.024
    dB at that null, so the 4.3 dB is not an MFS resolution error (the
    reference's mesh, or the 0.01% difference in length, are the candidates;
    not investigated). About 3 min per frequency for all angles.
31. **Never run two jobs at once.** Two Numba/BLAS jobs side by side made four
    QR factorisations take 453 s instead of 22 s and spoiled every timing of
    that batch. Queue jobs strictly one after another.
32. **Laptop.** Bempp is 3-4 times faster there than in the container (2592
    dofs: 22 s against 75 s). A separate worktree of the working branch is at
    `~/Projects/MFS-ccr` (branch `laptop-bem`); the old clone `~/Projects/MFS`
    was not touched. Results were committed there and pushed to the working
    branch.

## 3. Current best commands

```
python3 pec_sphere_convergence.py && python3 plot_pec_sphere.py
python3 fig_cost_cube.py ; python3 fig_cost_almond.py ; python3 fig_cost_sphere.py   # tables and figures from saved data (add `run` to recompute)
python3 fig_almond_rcs.py ; python3 fig_almond_benchmark.py                         # almond RCS against Bempp / against the Austin benchmark
python3 cube_general.py 6 60 25 35                # cube, any plane wave (from theta, phi; polarisation angle), four parity classes
python3 check_tail_exponents.py 5                 # edge and vertex tail variants at one depth
python3 cone_exponents.py ; python3 cone_exponents.py cube   # tip exponents of the almond cone; vertex exponents of the cube
python3 bem_optical.py 0.1,0.07,0.05,0.04         # Bempp, quasi-uniform cube meshes: field, cross-sections, time
python3 bem_cube_graded.py 8,12,16 2              # Bempp, meshes graded to the edges; add "8 10" for reference quadrature, --wave=60,25,35 for an oblique wave
python3 pec_cube_lightning3d.py 4,5,6 tag=NAME kinds=em q=0.5 gamma=2.0 col_p=3 col_q=0.5 beta=0.5 sym=1 tail=9   # single-knob runs (no vertex tail, no balanced refinement)
python3 diag_cube.py 5 "tail_n=9"                 # residual map by zone, optical theorem
python3 pec_almond_skeleton.py 1.19 --bem         # nose-on case only (see its header)
python3 almond_rcs.py 1.19 6:3 10 --s=0.5         # monostatic RCS, HH and VV, compares with stored Bempp
python3 proto2d_square.py TE ; python3 proto2d_ellipse.py ; python3 proto2d_square_tail.py TE   # 2D prototypes, seconds
```

Cube cost (second container, 4 cores): section 2, item 27. Almond: about
3 min per frequency for all look angles and both polarisations (6:3 layout);
the 8:5 layout takes 7 min.

Files: `mie.py` (reference); `pec_cube_lightning3d.py` + `line_kernel.py`
(cube library and solver for the standard wave; `Sources`, `SymSources`,
`TailLine`, `VertexTail`); `fig_cost_cube.py` (balanced refinement family
`balanced()`, cost table); `cube_general.py` (any wave); `diag_cube.py`
(residual map); `pec_almond_skeleton.py` (almond library; its main is the
nose-on case); `almond_rcs.py` (general almond solver); `pec_almond.py`
(geometry, mesh, and the confocal solver kept as the record of the failed
layout); `cone_exponents.py` (singularity exponents of cones and of the cube
vertex); `bem_reference.py`, `bem_optical.py`, `bem_cube_graded.py`,
`bem_almond_rcs.py` (Bempp); `austin_almond/` (third-party benchmark data,
own licence, see its README).

## 4. Task list (one per session)

Done in the first session of 2026-11-16 (section 2, items 8-21): T0 clean-up,
T1 cube stall, T2 almond robustness, T3 almond general incidence, T4 almond
at 7 and 9.92 GHz, T5 error-versus-time tables and figures.

Done in the second session of 2026-11-16 (items 22-32), at the author's
request to finish the list in one go:

- **T6 Literature check**: `LITERATURE.md`. arXiv full texts read; everything
  else checked against Crossref records or abstracts only, because publisher
  sites are blocked. No 3D edge-singular MFS for Maxwell or Helmholtz found;
  Kong and Shen (arXiv:2610.26540, October 2026) call the 3D edge and vertex
  case open. arXiv 2605.15527 is MFS for Maxwell in bi-periodic layered media
  (smooth interfaces), not a competitor.
- **T7 Outline and critique**: `PAPER_OUTLINE.md` (the manuscript itself was
  only read, at commit 4f8731f).
- **T8 Measured almond data**: done with the public Austin benchmark instead
  of digitised curves from Woo et al. (item 30).
- **T9 Cube vertices**: exponents computed, vertex tail implemented and
  adopted (item 25).
- **T10 Cube, general incidence**: `cube_general.py` (item 26).
- **T11 Tidy**: headers of both solver modules now say what is library, what
  is the standard case and where the general solver is. Files were *not*
  renamed (every script imports them by name).

T3 check of the thesis case (wave along -z, polarised along y, 1.19 GHz):
MFS backscatter 0.110432 m^2 (row spacing 0.5) and 0.110436 m^2 (0.4); Bempp
0.109478 (4212 dofs) and 0.110008 m^2 (9558 dofs), converging towards the MFS.

Open, in the order I would take them:

- [ ] **T12 Read what could not be fetched** (author or a session with
      publisher access): the list at the end of `LITERATURE.md`. Most urgent:
      Ye et al. 2023 (regularised MFS for 3D PEC), De Smedt and Van Bladel
      1986 (compare their exponent tables with ours), Kelley et al. 2020
      (measurement uncertainty to quote), Liu's thesis.
- [ ] **T13 Write-up** (author writes; section 6 and `PAPER_OUTLINE.md`).
      Still missing for it: sketches of the two source layouts (cube
      cross-section at an edge; almond mid-plane), and a figure of the
      boundary residual near an edge with and without tail.
- [ ] **T14 Re-time the sphere and almond cost tables (optional).** They were
      timed in the first container; the cube table in the second, which is
      about 20% slower. Each table is consistent in itself. Re-run
      `fig_cost_sphere.py run` and `fig_cost_almond.py run` if one machine
      for all tables is wanted.
- [ ] **T15 Vertex treatment beyond the diagonal (optional).** The vertex
      tail carries only the radial law; a tail on the three bisector planes
      near the vertex (with the angular structure of the vertex singularity)
      might change the rate rather than a constant. Not needed for the paper.
- [ ] **T16 A compressed or higher-order BEM as competitor (optional).**
      Bempp was used with dense assembly; item 27 says what that does and
      does not show.
- [ ] **T17 Non-convex polyhedra (optional, new work).** The skeleton
      theorem and the layouts cover convex bodies only.
- [ ] **T18 Rename the solver files (optional)**, e.g. `cube_mfs.py`,
      `almond_mfs.py`, once nothing else is in flight.

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

- Edge tail exponents, derived set against the old one: J = 5 0.314% / 0.327%, J = 7 0.049% / 0.049%; derived second terms: no change.
- Tip cone exponents (almond): nu = 0.29959 the only one below 1; single exponent 2.2996 equals the fitted set of four (residual 0.24%, RCS within 0.0002 dB).
- Cube vertex exponents: nu = 0.45417 (known value reproduced), tau = 0.8400. Vertex tail: 0.31 / 0.133 / 0.049% -> 0.20 / 0.086 / 0.029% (J = 5, 6, 7).
- Graded Bempp, grading scan at 2592 dofs: g = 1, 2, 2.5 (q 6/8), 3: 0.63%, 0.15%, 0.20%, 4.8%. Quadrature 4/6, 6/8, 8/10 at 4608 dofs: 0.104%, 0.060%, 0.057%.
- Graded Bempp reference (q 8/10), 4608 / 7200 / 10368 dofs: 0.058%, 0.032%, 0.022% from MFS J = 9; extrapolated 0.016%.
- Cube cost table re-timed in one container: MFS J = 3..8 22-284 s for 1.27%-0.011%; Bempp uniform 63-832 s for 0.68%-0.19%; graded 32-161 s for 0.44%-0.105%.
- Cube, oblique wave, J = 5, 6, 7: 86, 176, 370 s for 0.27%, 0.112%, 0.055% against graded Bempp; standard wave reproduced to 2e-15.
- Almond against the Austin benchmark, 3.5-10.25 GHz: 0.04-0.14 dB (HH), 0.01-0.02 dB (VV) from the reference simulation; about 1 dB from measurements, as the reference is.
- Almond 10.25 GHz, layouts 6:3 and 8:5: 0.0006 dB mean difference, 0.024 dB in the -66.7 dBsm null where the reference simulation differs by 4.3 dB.
- Two jobs at once: QR 20 times slower, timings of that batch discarded.

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
  section 2 items 19-21 and 24-30 (`fig_cost_*.py`, `fig_almond_rcs.py`,
  `fig_almond_benchmark.py`, `cube_general.py`, `cone_exponents.py`,
  `check_tail_exponents.py`). Figures exist as PDF and PNG in `validation/`.
  Outline, critique and likely referee questions: `PAPER_OUTLINE.md`.
- What the paper must not claim, or must qualify (state of each caveat of
  the first session after the second):
  - *Measured almond data*: resolved. Compared with the public Austin
    benchmark at four frequencies (item 30). Say that the measured curves
    differ from both simulations by about 1 dB; quote the measurement
    uncertainty from Kelley et al. 2020 once read (T12). Not compared with
    Woo et al.'s original 1993 curves.
  - *Symmetry in the timings*: resolved as a statement, not as an advantage.
    General incidence is implemented for both shapes; for the cube it costs
    1.6-1.9 times the symmetric case (item 26). Bempp uses no symmetry and
    needs none.
  - *Independent cube reference*: improved from 0.05% to 0.02%. Graded Bempp
    with reference quadrature is 0.022% from the MFS field and its
    extrapolation 0.016% (item 29). Below that the reference is still the
    MFS itself (self-convergence 0.011% at J = 8, optical mismatch 4e-5).
  - *Exponents*: resolved. Edge, tip and vertex exponents are derived or
    computed (`DERIVATIONS.md` B-D), and the fitted or guessed sets were
    replaced with no loss. The derivations are written out but have not been
    read by anyone else; Timo should check B.
  - *Skeleton completeness*: now a proof for convex polyhedra
    (`DERIVATIONS.md` A), elementary but not found in the literature; for the
    almond (smooth, non-polyhedral) it remains the standard continuation
    argument plus numerical evidence, not a theorem.
  - *Novelty*: searched (T6). No precedent found for a 3D MFS with
    edge-singular sources for Maxwell or Helmholtz. Safe wording and the
    closest prior work are in `LITERATURE.md` section 1. Several key papers
    were seen as abstracts only (T12).
  - *Bempp settings*: partly resolved. Graded meshes and higher quadrature
    were tried and are reported; assembly is 98% of Bempp's time, so the
    linear solver is not the issue. Not tried: matrix compression,
    higher-order elements, other BEM codes (T16).
  - *New caveat*: on the cube the MFS does **not** beat a well-graded BEM
    mesh at 0.1-0.4% for a general incident wave (item 26), and is level with
    it at 0.15-0.2% in the symmetric case (item 27). Its advantage there is
    the convergence rate (0.011% in under 5 min), not the cost at engineering
    accuracy. On the sphere and the almond the advantage is large and not in
    doubt.
  - *New caveat*: single deep nulls of the almond RCS differ from the
    published reference simulation by up to 4.3 dB while the MFS is
    self-converged to 0.02 dB there (item 30). Report it; do not claim the
    MFS value is the right one without a third computation.
  - Everything is at most about 9 wavelengths in size, and the layouts are
    built by hand per shape.
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
  Bempp with more than about 5k dofs, where it is 3-4 times as fast
  (2592 dofs: 22 s against 75 s).
  Run commands there as `/home/sap/Projects/MFS/.venv/bin/python3 -W ignore
  script.py ...` from `~/Projects/MFS/validation`; background jobs there do
  not notify the session, so read their output file to check.
- **Session restarts.** The cloud worker was restarted several times in one
  day. Files in the repository survive; the scratch folder and running jobs
  do not. Put anything worth keeping in `validation/` and commit often.
  Long jobs: start them in the background and write results to a file.
- **Bempp timings** must exclude its one-off compilation: do a small warm-up
  solve first (the `fig_cost_*` scripts do).
- **Network** (second session): arxiv.org (PDFs and the export API),
  api.crossref.org and api.openalex.org (heavily rate-limited: one request
  every 10 s at best) work from the container; public GitHub repositories can
  be cloned. Publisher sites (sciencedirect, dl.acm.org, mdpi, and by
  implication ieeexplore) are blocked. `pip install` of packages that are
  not in CLAUDE.md's list was refused by the permission system (scikit-fem);
  write the few lines needed instead.
- **One job at a time, really** (item 31). Chain jobs in one background
  command; wait for a process with `tail --pid=PID -f /dev/null`.
- **Two containers, two speeds.** The second container of 16 November was
  about 20% slower than the first for both the MFS and Bempp. Never mix
  timings from different sessions in one table; `fig_cost_cube.py` data and
  all cube Bempp files were re-timed together.
- **Laptop worktree**: `~/Projects/MFS-ccr` on `flowx13` is a git worktree of
  the working branch (local branch `laptop-bem`), made so that the old clone
  `~/Projects/MFS` (detached at efac0ae with hand edits) stays untouched. It
  uses the old clone's `.venv`. Remove it with `git worktree remove
  ../MFS-ccr` from `~/Projects/MFS` when no longer wanted.
- **Scratch outputs** named `*_t_*.npy/json` and `almond_rcs_*.npz` are
  ignored by git on purpose.
- **Figures** follow one style (colour plus marker plus dash per method, so
  they print in greyscale); see `STYLE` in the `fig_*.py` scripts.
