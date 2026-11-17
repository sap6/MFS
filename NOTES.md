# Project notes

Last updated: 2026-11-17 (overnight session: Bempp and MFS both at their best, cube solver 25 times faster, first draft of the paper). All scripts are in
`validation/`. Companion files: `DERIVATIONS.md` (proofs and exponent tables),
`LITERATURE.md` (T6), `PAPER_OUTLINE.md` (T7). Problem throughout:
exterior scattering of a plane wave by a perfectly conducting (PEC) body,
tangential E_sc = -tangential E_inc on the surface, least-squares collocation.

## Start here (state at the end of 2026-11-16)

- Everything is committed and pushed on branch `ccr-116c520c-fts35f` (last
  commit of the session: "Notes: handover"). No pull request exists. No
  jobs are running anywhere.
- All experiment tasks on the list are done (T0-T11). What is left is
  reading, writing and optional extras: section 4, T12-T18. The next real
  step is **T13, the write-up**, which the author leads; `PAPER_OUTLINE.md`
  is the starting point and section 6 lists what the paper may and may not
  claim.
- Before changing any number in a table, read section 2, items 22-32: the
  cube table (item 27) was re-timed in one container and must be redone as a
  whole if any part of it is re-run.
- The one result that cuts against the paper's first draft of its story:
  on the cube an edge-graded Bempp mesh is as fast as the MFS at 0.15-0.2%
  and faster for a general incident wave at 0.1-0.4% (items 26-27).
- Open questions for the author (also at the end of `PAPER_OUTLINE.md`):
  keep or delete the third-party benchmark files in
  `validation/austin_almond/`; which BEM curve leads the cube figure; title.
- Nothing was ever pushed to `sap6/mfs-paper`; it is still the thesis text
  (commit 4f8731f).

### Fourth session, 2026-11-17, overnight (read this first)

The author asked: is Bempp compared with the MFS when both are at their
best; clear up open questions; start a first draft. Mid-session instruction:
"optimise our MFS so it does actually work significantly better than Bempp
where we show it".

- **Bempp had not been at its best.** It was assembling with its Numba
  backend. With `pip install pyopencl pocl-binary-distribution` it uses
  OpenCL and assembles about 20 times faster (items 39-40). Every Bempp time
  in items 13-36 is therefore obsolete; the *errors* are unchanged.
- **The cube MFS is 25 times faster** with identical fields (item 41):
  new line quadrature (`line_tree.py`), three mirror planes instead of two
  (`cube_mfs.py`), semi-normal equations (`lsq.py`). The almond solver is
  twice as fast.
- **All cost tables were re-timed in one container** with both codes in their
  best configuration (item 43; `retime_all.sh`). Use only those numbers.
- **On the cube Bempp was also given the mirror symmetry** that the MFS uses
  (item 55). On that equal footing, with the final MFS timings of item 56,
  Bempp is faster above 0.06%, the two are level from 0.06% to 0.04%, and
  the MFS is faster below 0.03% (1.4-3 times) and continues to an estimated
  0.003%. Do not write that the MFS is faster on the cube without saying
  which Bempp (with or without symmetry) and at what error.
- **The cube MFS numbers to use are those of item 56** (base grid 6 / 0.15
  from J = 6 on; `paper/tables/cube_mfs.tex`, `cube_oblique.tex`). The MFS
  times and errors in items 43, 47 and 55 are from the finer base grid and
  are superseded; the Bempp numbers there stand.
- **The optical theorem is an indicator, not a bound** (item 56): two fits
  with a too coarse interior grid were wrong by 2.4% and 0.03% with
  mismatches of 1e-5 and 3e-5.
- **The cube family with one source depth per octave stalls near 1.3e-4;
  two depths per octave converge** to 3.5e-5 with no stall (item 47). The
  reference for the cube is now the dense family at J = 10
  (`fig_cost_cube.npz`, key `dense_E10`). Against it the old J = 9 reference
  was wrong by 1.25e-4 and "0.011% at J = 8" was really 0.019%. All cube
  errors in items 8-44 are against the old reference: use item 47.
- **First draft**: `paper/` (LaTeX, compiles with `pdflatex`; `main.pdf` is
  committed). Tables are generated (`paper/make_tables.py`). Item 47.
- Do not use the old two-mirror path (`SymSources`, `fig_cost_cube.run_mfs`
  before this session) for timings; it still works and gives the same field.

### Third session of 2026-11-16

- New example: the **torus** (section 2, items 33-36; `pec_torus.py`,
  `torus_singularities.py`, `fig_cost_torus.py`, `fig_torus_sources.py`,
  `bem_torus.py`; `DERIVATIONS.md` E). Plain MFS (dipoles on an inner torus)
  reaches 7e-12 in 10 s; Bempp is at 0.8% after 334 s. Recommendation: put
  it in the paper as a short section between the sphere and the almond; the
  author decides.
- Second literature pass with full texts (item 37, `LITERATURE.md` section
  7): the 3D edge case is still open in the literature; new items to cite;
  two claims to soften (skeleton theorem: Gustafsson's mother bodies are the
  potential-theory analogue; torus: MAS was applied to a torus in the 1990s).
- An illustrated explainer for the author was published as a private
  artifact ("MFS Thesis to Paper"); its source is `explainer/index.html`.
- Still open: T12 (papers behind publisher walls, list in `LITERATURE.md`
  section 8), T13 write-up.

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
| Unit cube, k = 2pi | edge sources + edge and vertex tails, **one depth per octave** (`cube_mfs.balanced`) | 1.28%, 0.53%, 0.21%, 0.094%, 0.039%, 0.019%, 0.0125% for J = 3..9 (1.2 to 14 s), then stalls | MFS dense J = 10 |
| Unit cube, k = 2pi | the same with **two depths per octave** (`cube_mfs.dense`) | 1.39%, 0.51%, 0.20%, 0.081%, 0.032%, 0.012%, 0.0035% for J = 3..9 (2.7 to 45 s); optical mismatch down to 5.5e-6 at J = 10; no stall | its own J = 10 (estimated error 2.5e-5); graded Bempp with symmetry (28800 dofs, quadrature 8/10) is 0.010% away; cross-sections agree to 4e-5 |
| Unit cube, oblique wave (from theta 60, phi 25, polarisation angle 35 degrees) | the same, eight parity classes | one depth: 0.62% to 0.017% for J = 4..9 (2.2 to 27 s); two depths: 0.014% at J = 8 | MFS dense J = 10 for that wave; graded Bempp with symmetry, 28800 dofs, is 0.013% away |
| Torus R = 1, a = 0.5, k = 2pi (also a = 0.25, 0.8; any wave) | dipoles on the inner torus a_s = 0.5 a, solved mode by mode about the axis | field 5e-7 (8064 unknowns, 3 s), 2e-9 (13824, 7 s), 7e-12 (21120, 10 s) against the finest run; residual 1e-10, optical mismatch 9e-14 | finest MFS (40320 unknowns); Bempp 9600 dofs is 0.82% away and its 1/N extrapolation of the cross-section agrees to 7e-5 |
| NASA almond, 1.19 GHz, nose-on | confocal / shrunk source surface | unstable or 20%+ | - |
| NASA almond, 1.19 / 7 / 9.92 GHz, monostatic RCS, HH and VV | skeleton sources + tip tail, four parity classes | within 0.03 / 0.04 / 0.12 dB rms of Bempp; self-converged to 0.001 dB (1.19 GHz) and 0.01 dB rms (9.92 GHz) | Bempp 9558 / 13632 dofs |
| NASA almond, 1.19 GHz, nose-on | skeleton sources on mid-plane (0.028% in 2.1 s, 0.0056% in 3.2 s) | self-converged to 5e-5, optical mismatch 1e-4; Bempp approaches it like 1/N (0.67%, 0.29%, 0.20%) | Bempp 4212, 9558, 13632 dofs |
| NASA almond, 3.5 / 5.125 / 7 / 10.25 GHz, monostatic RCS, HH and VV | the same solver as the RCS row above, tip exponents from the cone | mean difference from the Austin reference simulation 0.04-0.14 dB (HH), 0.01-0.02 dB (VV); from the measurements 0.6-1.8 dB, the same as the reference simulation's own distance from them | Austin RCS benchmark III-A (Kelley et al. 2020) |

Bempp itself (OpenCL backend, one machine, item 43 and 47): sphere 0.37% with
9750 dofs (19 s); torus 0.82% (9600 dofs, 21 s); cube, quasi-uniform: 0.19%
(13302 dofs, 49 s); cube, graded to the edges: 0.106% (4608 dofs, 6.1 s),
0.029% (10368 dofs, quadrature 6/8, 56 s), 0.014% (14112 dofs, quadrature
8/10, 166 s); almond 0.30% (9558 dofs, 19 s). In every case Bempp converges
towards the MFS field.

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
    tau = 0.840, 1.20635. A vertex tail on the body diagonal (`VertexTail`,
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

### Third session, 2026-11-16 (items 33-38)

33. **Torus, plain MFS** (`pec_torus.py`). Electric dipoles (3 unknowns
    each) on rings of an inner torus; collocation on rings of the surface,
    two tangential components, rows weighted by sqrt(area). Both sets are
    tensor grids (ring, angle about the axis), so the least-squares problem
    splits exactly into one small problem per Fourier mode about the axis,
    for any incident wave (checked against the full matrix: 8e-15). Known
    exterior solution (three dipoles inside the tube): 4e-10 with 13824
    unknowns. Plane wave along the axis, R = 1, a = 0.5, a_s = 0.5 a: residual
    2.5e-2, 3.7e-4, 1.6e-6, 6.5e-9, 1.0e-10 for 1440, 3840, 8064, 13824, 21120
    unknowns; optical mismatch down to 9e-14; sigma_ext = 11.75257471.
    Oblique wave (from theta 60, phi 25, polarisation 35): residual 2.6e-10
    with 21120 unknowns, sigma_ext = 9.53615856. Not prototyped in 2D first:
    the 3D runs take seconds.
34. **Torus: the singular set** (`torus_singularities.py`,
    `fig_torus_sources.py`, `DERIVATIONS.md` E). The continued field is
    singular on the centre circle and on the focal ring rho = sqrt(R^2 - a^2);
    a concentric source torus needs a_s > R - sqrt(R^2 - a^2), i.e. a_s/a >
    0.127, 0.268, 0.5 for a/R = 0.25, 0.5, 0.8. Measured: distance 0.395-0.401
    for a/R = 0.8 with four waves (predicted 0.400), 0.49 from a centre moved
    by 0.1 (0.50); coefficient norm jumps from 3.6 to 7e4 between a_s/a = 0.5
    and 0.4. For a/R = 0.5 and 0.25 the fit is consistent but not clean. The
    reflection argument behind it is a sketch, not a proof.
35. **Torus: traps.** (a) An off-centre source circle with the default
    2 collocation rings per source ring gives garbage without warning
    (residual 1 to 1e4 at independent points): use 4 rings per source ring
    (`ntheta_col=4*ns`). (b) Fourier modes that the wave does not excite
    carry coefficients of order 1 (rounding noise in the null space); they
    do not harm the field but ruin any statistic of the coefficients unless
    excluded (`rhs_modes`). (c) The size of the coefficients on an off-centre
    circle is not a clean indicator of enclosure; the density-decay fit
    with a concentric circle is.
36. **Torus cost table** (`python3 fig_cost_torus.py`, third container, 4
    cores, error against the finest MFS run, times include residual and
    cross-section evaluation): MFS 1440, 3840, 8064, 13824, 21120 unknowns:
    0.28, 1.1, 3.1, 6.7, 10.1 s for 5.0e-3, 1.1e-4, 5.0e-7, 2.0e-9, 7.3e-12.
    The same problem with the full matrix (no symmetry): 1440, 2496, 3840
    unknowns: 1.7, 6.5, 25.5 s for 5.0e-3, 7.6e-4, 1.1e-4. Bempp (structured
    mesh, `bem_torus.py`) 864, 2400, 4704, 9600 dofs: 9, 33, 93, 334 s for
    10%, 3.3%, 1.7%, 0.82%; sigma_ext 11.255, 11.585, 11.668, 11.711, error
    like 1/dofs, extrapolated 11.7533 (MFS 11.7526). Oblique wave: Bempp
    4704 dofs is 2.0% from the MFS. The MFS uses the rotational symmetry of
    the body; Bempp uses none.
37. **Second literature pass** (`LITERATURE.md` sections 7-8). Full texts
    read this time. No 3D MFS/MAS with edge-singular sources found; the MAS
    literature rounds edges (Bogdanov, Karkashadze, Zaridze 1999). Gopal and
    Trefethen 2019 pose the 3D question explicitly. New to cite: Kangro 2012
    (convergence of interior source methods for 3D Maxwell), Hochman,
    Leviatan, White 2013 (JCP; rational-function / dipole fitting for 2D
    Laplace with corners, a forerunner of lightning), Gustafsson 1998 (mother
    bodies of convex polyhedra: the potential-theory analogue of the skeleton
    theorem), Karkashadze and Zaridze (MAS on a torus, auxiliary contour
    around the foci), Cai, Bagge, Martinsson 2026 (MFS on tori, acoustic),
    Yurkin and Kahnert 2013 (accuracy limits for a dielectric cube).
    arXiv:2310.01665 is Ginn's Oxford MMath dissertation (arXiv authors Ginn
    and Trefethen).
38. **Is the torus worth including?** Yes, as about one page: it is the one
    shape where the thesis's own layout (sources on an inner copy of the
    surface) works unchanged on a non-convex, ring-shaped body, it reaches
    1e-11, and it makes the paper's organising rule (sources must enclose the
    singularities of the continued field) visible with a formula that can be
    checked. It is not novel as an application (MAS and acoustic MFS have
    used tori) and it has no closed-form reference. Do not present the
    focal-ring statement as new before Millar 1983 has been read.

### Fourth session, 2026-11-17, overnight (items 39-55)

39. **Bempp was run with its slow backend.** bempp-cl assembles with OpenCL
    if `pyopencl` and a CPU driver are installed and falls back to Numba
    otherwise; the containers had neither. With
    `pip install pyopencl pocl-binary-distribution` (now in `CLAUDE.md`) the
    assembly is about 20 times faster with identical results
    (`bem_backend.py`: sphere, 4824 dofs: assembly 96 s -> 4.3 s, total 108 s
    -> 10 s; error 7.58e-3 both). Single precision changes neither the error
    nor, by more than 20%, the time; double is kept. All Bempp times of
    items 13-36 are obsolete, the errors are not.
40. **Three more things that held Bempp back** (all fixed, all in
    `bem_reference.efie_solve`, which every Bempp script now uses):
    (a) the right-hand side was a Numba closure defined inside the solver
    and so recompiled at every solve, about 3 s each time (`bem_overhead.py`:
    a 714-dof assembly takes 0.27 s once compiled). It is now one
    parameterised callable at module level. (b) PoCL takes its memory from
    one NUMA node as the kernel lists it (4 GB of 16) and refuses matrices
    above a quarter of that rounded up to a power of two: 1, 2 or 4 GB
    depending on the machine and on what else is running. `bem_env.py`
    (imported before bempp) tells it the real memory; the cap is then 4 GB =
    16384 dofs in double precision. (c) `bempp.linalg.lu` copies the matrix;
    the LU is now in place. Quadrature orders: default 4/6 limits graded
    cube meshes near 8e-4; 6/8 is the best trade; 8/10 for references.
    Not done: a boundary element code that uses the symmetry of the body,
    higher-order elements, compression.
41. **Cube MFS, 25 times faster, same fields** (changes verified against the
    stored fields to 1e-12 and then 1e-15 at every step):
    - `line_tree.py`: the line integrals. A binary tree of panels in the
      stretched variable with precomputed nodes and weights against all
      modes (moment rule with 20 nodes in t on the four coarsest levels, 10
      Gauss nodes in u below, deeper levels on the fly); a target accepts a
      panel when its length is below its distance. Electric and magnetic
      line share the kernel; the 9 tail lines share nodes and are summed
      before the modes; the 8 mirror images are summed inside the kernel so
      each matrix entry is written once. Accuracy 5e-10 against a finer
      adaptive rule (`test_line_tree.py`); looser settings (2e-6) change the
      field by 2e-7 and save little.
    - `cube_mfs.py`: three mirror planes (eight parity classes) instead of
      two; the standard wave needs two classes with half the unknowns each.
    - `lsq.py`: Householder QR without forming Q (LAPACK geqrf + unmqr with
      the optimal workspace); and, for well-conditioned matrices, the
      corrected semi-normal equations (Cholesky of the Gram matrix, two
      refinement steps, implicit column scaling, no copy of the matrix).
      cond(A) of a cube class: 3.8e4 (J = 6), 7.7e5 (J = 8). Falls back to
      QR by itself when the Cholesky factor looks ill conditioned.
    - The same factorisation in the almond solvers (2 times faster, RCS
      against Bempp unchanged: 0.023 / 0.029 dB rms at 1.19 GHz; all angles
      and both polarisations now take 43 s at 1.19 GHz).
    - Sphere: QR gives the same error as the SVD solver (2e-14 at N = 1200)
      and is 4 times faster.
    - Did not help: fewer collocation nodes. Two per panel instead of three
      gives 85% error at J = 5 (optical mismatch 0.35 is the only warning)
      and 0.2% at J = 6.
42. **The cube family has a floor near 5e-5** *(cause found and reference replaced: item 47)*. J = 10, 11, 12 differ from
    J = 9 by 5.0e-5, 7.5e-5, 1.0e-4 (from each other by 3-4e-5); variants at
    J = 10 (modes as at J = 9, four collocation nodes per panel, base grid
    10 with panels 0.08) differ from the family by 6e-5, 7e-5, 1e-5; the
    optical mismatch stops at 4-7e-5. Conditioning is not the cause
    (|c| about 1). So J = 9 is a reference good to about 5e-5, and the
    errors quoted for J = 7 and J = 8 carry that uncertainty. Suspect: the
    vertices (T15), not confirmed.
43. **Cost tables, final** *(cube errors here are against the old J = 9 reference; times stand; item 47 has the cube against the new reference)* (`validation/retime_all.sh`, one container, 4
    cores, 2026-11-17 02:25-02:46 UTC, one job at a time; Bempp with OpenCL,
    double precision, compiled right-hand side; MFS as in item 41; every
    table in `paper/tables/` is generated from these files by
    `paper/make_tables.py`). Times in seconds include assembly, solve and
    the field at the test points.
    - Sphere (against Mie): MFS 300, 600, 1200, 2400, 3600 unknowns: 0.03,
      0.08, 0.29, 1.7, 4.6 s for 2.7e-3, 2.7e-6, 2.0e-9, 1.7e-13, 2.3e-14.
      Bempp 1269, 2733, 4824, 9750 dofs: 0.77, 1.9, 4.4, 19 s for 3.2%,
      1.4%, 0.76%, 0.37%.
    - Torus (against the finest MFS): MFS 1440 to 21120 unknowns: 0.25,
      0.90, 2.4, 4.3, 6.3 s for 5.0e-3, 1.1e-4, 5.0e-7, 2.0e-9, 7.3e-12;
      full matrix 3840 unknowns: 19 s for 1.1e-4. Bempp 864, 2400, 4704,
      9600 dofs: 0.97, 1.5, 5.2, 21 s for 10%, 3.3%, 1.7%, 0.82%.
    - Almond, 1.19 GHz nose-on (against layout 8:5): MFS 795, 1377, 2103,
      2937 unknowns: 1.1, 1.5, 2.1, 3.2 s for 6.7%, 0.40%, 0.028%, 0.0056%.
      Bempp 2106, 4212, 9558 dofs: 1.0, 3.3, 19 s for 1.3%, 0.67%, 0.30%.
      At 1% the two are level; below, the MFS is far ahead.
    - Cube, standard wave (against MFS J = 9, 13 s): MFS J = 3..8: 1.3, 1.7,
      2.4, 3.9, 7.1, 9.4 s for 1.27%, 0.52%, 0.20%, 0.089%, 0.033%, 0.011%.
      No tails, J = 4, 6, 8: 0.9, 2.7, 6.9 s for 5.4%, 2.5%, 1.1%; plain MFS
      21-30%. Bempp quasi-uniform 2214 to 13302 dofs: 2.3, 6.8, 19, 49 s for
      0.68%, 0.39%, 0.26%, 0.19%. Bempp graded (g = 2), quadrature 4/6,
      1152 to 14112 dofs: 1.2, 2.8, 6.1, 14, 28, 56 s for 0.44%, 0.15%,
      0.105%, 0.092%, 0.084%, 0.077%; quadrature 6/8: 3.1, 7.1, 15, 30, 56,
      102 s for 0.44%, 0.13%, 0.060%, 0.038%, 0.030%, 0.026%; quadrature
      8/10, 10368 and 14112 dofs: 97 and 166 s for 0.022% and 0.018%.
    - Reading for the cube: above 0.2% level (Bempp 2.8 s for 0.15%, MFS
      2.4 s for 0.20%); at 0.1% the MFS is 1.6 times faster (3.9 against
      6.1 s); at 0.03-0.04% four to eight times (7.1 s against 30 s for
      0.038% and 56 s for 0.030%); Bempp's best, 0.018% in 166 s, the MFS
      passes in 9.4 s.
    - Cube, oblique wave (`cube_oblique_table.py`; against MFS J = 9 for
      that wave, 29 s, optical mismatch 9.9e-5): MFS, eight classes, J =
      4..8: 2.4, 3.4, 6.1, 12, 18 s for 0.61%, 0.27%, 0.10%, 0.039%, 0.011%.
      Bempp graded 4/6, 2592, 4608, 7200 dofs: 2.7, 6.4, 14 s for 0.21%,
      0.14%, 0.12%; 6/8, 4608, 7200, 10368 dofs: 15, 30, 57 s for 0.077%,
      0.049%, 0.039%. Same picture as the standard wave: Bempp slightly
      ahead above 0.2%, the MFS ahead below 0.1%.
    - Where the session started (item 27, slow Bempp, old MFS): cube MFS
      J = 6 and 8 took 95 and 284 s; graded Bempp 4608 dofs 161 s.
44. **What the cube comparison rests on (say it in the paper).** (a) The
    MFS uses the three mirror planes of the body (valid for any wave);
    Bempp uses none. A boundary element code that used them would need an
    eighth of the assembly, which is 40-90% of Bempp's time here. Without
    symmetry the dense MFS factorisation would cost 64 times as much as the
    eight class factorisations. (b) The MFS is not economical in unknowns
    on the cube: about 10^4 (all classes) for 0.09% against 4608 for the
    graded mesh. Its advantage is the rate (error / 2.3-2.9 per level for
    20-45% more unknowns) and the price of an unknown. (c) Errors below
    1e-4 are uncertain by 5e-5 (item 42).
45. **Sphere and wavenumber** (`sphere_scaling.py`, R = 0.5): sources needed
    for 1e-6 (interpolated): 224, 473, 797, 1628 for k = 2 pi, 4 pi, 6 pi,
    8 pi, i.e. 53, 28, 21, 24 unknowns per square wavelength of surface. At
    8 pi the error stalls near 3e-9 (coefficient growth) and the solver
    falls back to rank-revealing QR, which is slow (87 s at N = 2400).
46. **The almond null at 10.25 GHz, 35.5 degrees, is a matter of position**
    (`almond_null.py`, layout 6:3, HH). On a 0.02 degree grid the MFS null
    is at 35.66 degrees, -69.7 dBsm; the reference simulation, sampled every
    0.5 degrees, has -71.1 dBsm at 35.5 and its neighbours are lower on the
    left and higher on the right than the MFS, i.e. its null lies about
    0.1-0.15 degrees lower. Shifting the MFS curve by -0.10 degrees brings
    the rms difference over 33-38 degrees from 1.45 to 0.64 dB and the
    difference at 35.5 from 4.3 to 1.9 dB. A size (or frequency) change of
    0.01% moves the null by 0.02 degrees and changes the value at 35.5 by
    0.4 dB, so the 0.01% difference in length between our almond and the
    benchmark's explains a tenth of it; -0.05% gives rms 0.96 dB. So the
    two simulations agree on the null to about 0.1 degree; which position
    is right needs a third computation. The old caveat ("differ by 4.3 dB,
    not investigated") is replaced by this.
47. **Cube: the floor was the depth ratio.** (`cube_floor.py`,
    `cube_mfs.dense`, `fig_cost_cube.py`.) With two source depths per octave
    (ratio 2^(-1/2), same deepest depth 0.2 * 2^-J, four collocation nodes
    per panel) the stall disappears: at J = 9 the boundary residual falls
    from 7.6e-4 to 2.2e-4, the optical mismatch from 4.0e-5 to 1.3e-5, and
    two different collocation grids (4 nodes per panel with ratio 1/2; 3
    with ratio 2^(-1/2)) agree to 9.7e-7, where variants of the old family
    differed by 6e-5. Three nodes per panel with the dense depths is
    garbage (error above 90%, optical mismatch 0.5, |c| 7000-22000): the
    sampling rule of item 3 again. Dense family J = 3..10, standard wave:
    change to the next level 9.0e-3, 3.1e-3, 1.2e-3, 4.9e-4, 2.0e-4, 8.4e-5,
    3.5e-5 (ratio 2.4-2.9, no stall); optical mismatch 3.1e-3 down to
    5.5e-6; sigma_ext 2.927855, 2.927930, 2.927961 at J = 8, 9, 10 (limit
    about 2.92798). So the J = 10 reference is good to about 2.5e-5.
    Against it (times from `fig_cost_cube.py run`, same machine as item 43):
    - one depth per octave, J = 3..9: 1.2, 1.6, 2.4, 3.7, 7.1, 9.9, 13.6 s
      for 1.28%, 0.53%, 0.21%, 0.094%, 0.039%, 0.019%, 0.0125%; J = 10 is
      also 0.0125% (stalled);
    - two depths per octave, J = 3..9: 2.7, 3.7, 6.0, 10.7, 21.8, 31, 45 s
      for 1.39%, 0.51%, 0.20%, 0.081%, 0.032%, 0.012%, 0.0035%; J = 10: 64 s;
    - no tails (one depth), J = 4, 6, 8: 5.4%, 2.5%, 1.1% (so the tails are
      a factor 60 at J = 8, not 100);
    - Bempp graded, quadrature 4/6: 0.152%, 0.106%, 0.085% (2592, 4608,
      10368 dofs); 6/8: 0.060%, 0.037%, 0.029%, 0.025% (4608 to 14112);
      8/10: 0.031%, 0.0195%, 0.0143% (7200, 10368, 14112 dofs; 56, 97,
      166 s); extrapolated from these three 0.0098%, sigma_ext 2.927864.
      So the independent confirmation of the cube field is now 1.4e-4
      directly and 1.0e-4 by extrapolation, and the cross-sections agree to
      4e-5.
    - Reading: level above 0.2%; MFS (one depth) 1.6 times faster at 0.1%,
      4 times at 0.04%, 10 times at 0.02% (9.9 s against 97 s); Bempp's
      best (0.014%, 166 s) is passed in 14 s (one depth) or 31 s (two);
      two depths continue to 0.0035% in 45 s.
    - Oblique wave (`cube_oblique_table.py`, reference dense J = 10, optical
      mismatch 1.1e-5): one depth, J = 4..9: 2.2, 3.3, 6.0, 12.6, 18.2,
      27.1 s for 0.62%, 0.28%, 0.117%, 0.052%, 0.026%, 0.017%; two depths,
      J = 4..8: 5.3, 9.7, 19.7, 46, 105 s (the last two in class groups) for
      0.57%, 0.24%, 0.098%, 0.039%, 0.014%; Bempp graded 4/6: 2.7, 6.4,
      13.7 s for 0.21%, 0.144%, 0.122%; 6/8: 15, 30, 57 s for 0.082%,
      0.052%, 0.040%. Bempp slightly ahead above 0.2%, level at 0.14%, MFS
      2.4 times faster at 0.05% and alone below 0.04%.
    - Memory: eight classes of the dense family do not fit together from
      J = 8 (11.5 GB); `cube_mfs.run` splits the classes into groups
      (`max_gb`). One run was killed by the kernel for lack of memory before
      the size estimate was fixed.
48. **Figures made in this session** (all in `validation/`, print style):
    `fig_layouts.py` (cube: section through an edge, depths, stretched
    variable; almond: mid-plane rows, focal curve, tip tail);
    `fig_edge_residual.py` (residual on a face against the distance to an
    edge, J = 6, with and without tail: without, 1.6 at every distance below
    1e-3; with, 1.2e-2 for 1e-4..1e-3 and 0.13 for 1e-5..1e-4; below the last
    tail line, d_J/512, it rises again); `fig_sphere.py` (replaces the old
    six-panel PNG); `fig_cost_cube.py` now with both MFS families and graded
    Bempp at two quadrature settings.
49. **First draft of the paper**: `paper/` (see section 6). 40 pages in
    preprint format, 24 in the journal's layout (`main_3p.pdf`). `python3 make_tables.py && pdflatex main && bibtex
    main && pdflatex main && pdflatex main`. Every table body is generated
    from the data files; numbers in the running text were typed from the
    same output and must be re-checked if any table is re-run
    (`make_tables.py` prints them). `check_bib.py` compares every title in
    `refs.bib` with its Crossref record (all agree; four were corrected);
    the arXiv entries were checked against the arXiv API.
50. **Writing rules used for the draft** (from the author's thesis and the
    academic-humanizer skill, github.com/AIScientists-Dev/academic-humanizer):
    "we", British spelling, define before use, explain each formula in a
    sentence; no em-dashes; sentences under about 30 words; every claim
    carries its number, figure or citation; no "novel", "robust",
    "comprehensive", "paves the way"; verbs no stronger than the evidence
    ("show", "reach", "agree to"); hedges only where the result is uncertain
    (torus bound, skeleton statement for smooth bodies); limits stated in
    their own section.
51. **The almond's self-convergence is real** (`almond_floor.py`, 1.19 GHz
    nose-on; done because the cube's one-per-octave family had hidden a
    1.2e-4 error). The reference layout 8:5 changes by 8.3e-6 with row
    spacing 0.5 instead of 0.7, by 3.3e-6 with a finer collocation grid
    (s_col 0.25), and by 2.7e-6 with two rows per octave towards the focal
    curve (8:10, q = 2^-1/2). So 8:5 is good to about 1e-5 and "6:3 is
    5.6e-5 from it" is a true error. Not a valid member: 10:7 (under-sampled,
    coefficient norm 4e3, 4e-5 away; it falls back to QR and takes 54 s).
    The optical mismatch of 8:5 is 1.5e-4 in all variants, larger than the
    field error: it is an indicator, not a bound. Cases with more than about
    10000 unknowns per class ran out of memory.
52. **How robust is the cube reference?** (`cube_floor.py d9 ...`, dense
    family at J = 9, differences of the field from the family's own run.)
    More modes (n0 = 18 instead of 14): 6.3e-7. Tails of 14 lines instead of
    9: 6.3e-7. A second vertex exponent: 1.7e-7. Another collocation grid
    (3 nodes per panel, panels shrinking by 2^(-1/2)): 9.7e-7. Row-weight
    exponent beta = 0.3: 5.9e-5 (and further from J = 10: 9.2e-5; residual
    1.1e-3); beta = 0.7: 2.1e-5 (and closer to J = 10: 2.0e-5 instead of
    3.5e-5; optical mismatch 8.5e-6 instead of 1.3e-5). So only the level J
    and the weighting matter at the 1e-5 level, and beta = 0.7 may be a
    slightly better choice than 0.5 for the far field (not adopted: it would
    mean re-running every cube table; T24). A finer base grid (12, panels
    0.07) at J = 9 needs more than the 15 GB and was killed.
53. **The depth ratio sets a floor in 2D too** (`proto2d_tapered.py`,
    `proto2d_depth_ratio.py`; unit square, point sources on the corner
    bisectors). Free solves down to depth 1e-9, TE, far field against a
    converged solve: ratio 1/2 (30 depths): 8.1e-5; ratio 0.6 (40): 5.8e-6;
    ratio 2^(-1/2) (58): 1.4e-9; ratio 0.8 (89): 1e-11. With tails instead
    of deep depths the ratio-1/2 error is the same 8.4e-5 from a deepest
    depth of 4e-3 down to 2e-5, with one tail exponent per source type or
    six. So a fixed geometric ratio gives a fixed attainable accuracy, as
    the theory of lightning solvers says, and the 3D floor of 1.3e-4 for one
    depth per octave is the same effect. (The collocation must follow: two
    and three times as many nodes per panel for ratios 2^(-1/2) and 0.8.)
    Also found: the reference solve of `proto2d_square_tail.py` (70 depths,
    down to 1e-16) is wrong by 85% with numpy 2.5 / scipy 1.18 (far-field
    norm 49.75 instead of 35.95); 30, 40 and 50 depths agree to 3e-7. That
    script now uses 40.
54. **Tapered depths: a large gain in 2D, none in 3D** (T21, closed).
    Lightning-type spacing d_j = d0 exp(-sigma (sqrt(n) - sqrt(j))), dense
    near d0 and sparse near the corner. 2D (free solves, TE): 30 tapered
    depths give 9.5e-8 where 30 geometric ones (ratio 1/2) give 8.1e-5; 40
    tapered equal 58 geometric (2^(-1/2)) at 1e-9. 3D cube (`cube_mfs.tapered`,
    with tails, error against the reference): with as many depths as one
    per octave it is no better at J = 6 (9.4e-4) and worse at J = 8 and 9
    (4.4e-4 and 4.3e-4 against 1.9e-4 and 1.25e-4); 13 and 16 tapered
    depths at J = 9 give 1.27e-4 and 5.3e-5, where 19 uniform ones give
    3.5e-5 (34 s against 45 s). So with the tail in place the resolution
    near the deepest depth matters and the uniform ratios stay. Three
    collocation nodes per panel or four made no difference for the tapered
    layouts.
    Mixed layouts at J = 9 (two per octave only above, or only below, a
    depth; error against the reference): dense above 0.2 * 2^-6 (16
    depths, 35 s): 3.8e-5; above 2^-4 (14 depths, 31 s): 5.6e-5; dense only
    below 2^-4 (15 depths): 7.7e-5; only below 2^-6 (13 depths): 9.3e-5;
    uniform two per octave (19 depths, 45 s): 3.5e-5. So the spacing matters
    at all depths, a little more at the shallow ones, and a mixed layout
    saves about 20%: not worth a third family.
55. **Bempp with the cube's mirror symmetry: the comparison on an equal
    footing** (`bem_cube_sym.py`, T22 done). The mesh is built invariant
    under the eight reflections; a reflection is a signed permutation of the
    RWG coefficients (found numerically from the basis functions, so no
    orientation convention is assumed); only the rows of the octant's test
    functions are assembled (Bempp: `support_elements`), folded with the
    symmetrised basis, and each parity class is solved by a small LU. It
    reproduces the full Galerkin solution to 2.5e-5 (quadrature 4/6) and
    1.6e-6 (8/10): Bempp's quadrature is not reflection-invariant. It is 3.5
    to 8 times faster than the full solve and reaches 28800 dofs (41472
    failed, probably memory). Errors against the MFS reference:
    - standard wave, quadrature 4/6: 0.7, 1.3, 3.7 s for 0.152%, 0.107%,
      0.087% (2592, 4608, 10368 dofs); 6/8: 4.3, 5.8, 9.0, 13, 21, 43 s for
      0.060%, 0.037%, 0.029%, 0.026%, 0.023%, 0.020% (4608 to 28800);
      8/10: 9.0, 15, 22, 34, 69 s for 0.031%, 0.020%, 0.0145%, 0.012%,
      0.0100% (7200 to 28800). The convergence flattens beyond 14112 dofs
      (thin triangles at the edges, or the quadrature; not investigated).
    - **Reading, both with symmetry**: Bempp is faster above 0.05% (1.3 s
      against 3.7 s at 0.1%), level near 0.03% (5.8-9.0 s against 7.1 s),
      and the MFS is faster below 0.02% (9.9 s against 15 s at 0.019%; 14 s
      against 34 s at 0.012%); Bempp's best is 0.010% in 69 s, the MFS
      (two depths) reaches 0.0035% in 45 s.
    - oblique wave: 4/6: 0.8, 1.7 s for 0.21%, 0.144%; 6/8: 3.3, 5.0, 9.3 s
      for 0.082%, 0.052%, 0.040%; 8/10: 15, 24, 40, 91 s for 0.027%, 0.020%,
      0.016%, 0.013%. Against the MFS (6.0 s 0.117%, 13 s 0.052%, 18 s
      0.026%, 27 s 0.017%): Bempp faster down to 0.03%, level below.
    - The independent confirmation of the cube field is now 1.0e-4 directly
      (28800 dofs).
    - **What this means for the paper**: on the cube the cost advantage of
      the MFS over a graded low-order BEM is real only below about 2e-4 and
      is a factor of 1.5-2.5 there; against Bempp as distributed (no
      symmetry) it is 1.6 to 10 times from 1e-3 down. The draft states
      both. The cube's selling points are convergence where plain MFS
      fails, the rate, and reaching 3.5e-5; not speed at 0.1%.
      Item 44 (a) is superseded by this measurement.

56. **Last tuning of the cube family, and two cautions** (inline runs and
    `fig_cost_cube.py run`, `cube_oblique_table.py run`; adopted in
    `cube_mfs.balanced`).
    - *Fewer modes do not pay*: n0 - 2, n0 - 4 or a slower growth of the
      modes with depth all increase the error at equal J.
    - *The base grid was over-resolved from J = 7 on.* (n_base, hmax) =
      (6, 0.15) instead of (8, 0.1) gives the same fields for the dense
      family (difference 1.4e-6 to 2.2e-6 at J = 7..10) and is 25% faster.
      Adopted for J >= 6 in both families. All cube tables were re-run:
      - standard wave, one depth per octave, J = 3..9: 1.2, 1.6, 2.5, 3.7,
        5.2, 7.7, 10.5 s for 1.28%, 0.53%, 0.21%, 0.094%, 0.036%, 0.016%,
        0.015%;
      - two depths per octave, J = 3..10: 2.3, 3.7, 5.5, 10, 16, 24, 36,
        50 s for 1.39%, 0.51%, 0.20%, 0.081%, 0.032%, 0.012%, 0.0036%,
        reference (optical mismatch 5.5e-6, sigma_ext 2.927964; J = 8, 9:
        2.927854, 2.927930);
      - oblique wave, one depth, J = 4..9: 2.2, 3.4, 6.7, 9.8, 14, 22 s for
        0.62%, 0.28%, 0.117%, 0.053%, 0.027%, 0.019%; two depths, J = 6, 7,
        8: 20, 34, 71 s for 0.098%, 0.039%, 0.014%.
    - *Against Bempp with symmetry (item 55)*: Bempp faster above 0.06%
      (1.3 s against 3.7 s at 0.1%); level 0.06%-0.04% (4.3 s for 0.060%,
      5.8 s for 0.037% against 5.2 s for 0.036%); MFS faster below 0.03%:
      7.7 s for 0.016% against 15 s (0.020%) and 22 s (0.0145%); 24 s
      against 34 s at 0.012%; Bempp's best 0.010% in 69 s, MFS under 30 s.
      Against Bempp without symmetry: 1.6 times at 0.1%, 6 times at 0.04%
      (5.2 s against 30 s), 12 times at 0.02% (7.7 s against 97 s).
      Oblique wave: Bempp-sym faster down to 0.03% (5.0 s against 9.8 s at
      0.052%), level below (15 s against 14 s at 0.027%; 40 s for 0.016%
      against 22 s for 0.019%).
    - **Caution 1: the coarser base grid is not enough for the one-depth
      family beyond J = 9.** With (6, 0.15): J = 10, 11, 12 give 2.9e-4,
      4.6e-4, 9.0e-4 (rising); with (8, 0.1) they gave 1.3e-4 to 1.8e-4
      (item 47). J = 9 itself is 1.50e-4 against 1.25e-4. The family is
      only recommended down to 2e-4 (J <= 8), where the two grids agree.
      If the one-depth family is ever run deeper, use (8, 0.1). The dense
      family (four nodes per panel) is unaffected up to J = 10.
    - **Caution 2: the optical theorem missed two bad fits.** Base grid
      (4, 0.25) at J = 8: error 2.4%, optical mismatch 1.2e-5 (at J = 7:
      0.6% and 1.1e-3). One-depth J = 10 with (6, 0.15): error 2.9e-4,
      mismatch 3.0e-5. The draft said the mismatch "exposed every
      under-sampled fit"; corrected to "the fits under-sampled at the
      edges" plus the counter-example. Always check a new parameter set
      against the reference field, not only against the optical theorem.
    - **Caution 3: errors against a self-reference.** The dense family's
      level-to-level change falls by 2.3-2.9 per level (8.95e-3, 3.11e-3,
      1.19e-3, 4.86e-4, 2.04e-4, 8.42e-5, 3.61e-5). With that rate the
      true errors of J = 9 and J = 10 are about 6e-5 and 3e-5; the table
      entry 3.6e-5 for J = 9 is its difference from J = 10 and understates
      its error. The draft now says "estimated error 3e-5" for the finest
      run (abstract, introduction, highlights, conclusion) instead of
      "converges to 3.5e-5".
    - Edge appendix: the next terms of the wall currents have exponents
      larger by 2 (electric) and 4/3 (magnetic) (8/3, 11/3, 8/3, 11/3;
      `DERIVATIONS.md` B), not "by two or more" as the draft said.

## 3. Current best commands

```
python3 pec_sphere_convergence.py && python3 plot_pec_sphere.py
python3 pec_torus.py ; python3 torus_singularities.py                              # torus: self-checks; where the continued field is singular
python3 fig_cost_torus.py ; python3 fig_torus_sources.py                           # torus: cost table against Bempp; residual and |alpha| against the source radius (add `run`)
python3 bem_torus.py 12,20,28,40 0.5                                               # Bempp on the torus (optional: theta phi pol of the wave)
python3 fig_cost_cube.py ; python3 fig_cost_almond.py ; python3 fig_cost_sphere.py   # tables and figures from saved data (add `run` to recompute)
python3 fig_almond_rcs.py ; python3 fig_almond_benchmark.py                         # almond RCS against Bempp / against the Austin benchmark
python3 cube_mfs.py 4,5,6,7,8 [60 25 35]          # cube, the solver: standard wave, or any plane wave (from theta, phi; polarisation angle); eight parity classes
python3 cube_oblique_table.py                     # cube, oblique wave: both MFS families and graded Bempp against the dense reference (add `run`, 12 min)
python3 cube_floor.py J10 dense_p4 dense_colq     # what limits the one-depth-per-octave family (item 47)
python3 test_line_tree.py 5                       # accuracy of the fast line quadrature against the adaptive rule
python3 fig_layouts.py ; python3 fig_edge_residual.py 6 ; python3 fig_sphere.py   # layout sketches; residual near an edge; sphere convergence figure
python3 almond_null.py 6:3 fit                    # the HH null at 10.25 GHz: position, depth, best-fitting shift
python3 sphere_scaling.py                         # sphere: sources needed against the wavenumber
python3 bem_backend.py 0.14,0.1 ; python3 bem_overhead.py   # Bempp: Numba against OpenCL, single against double; fixed cost per assembly
./retime_all.sh                                   # every cost table on this machine, one job after another (20 min); then cd ../paper && python3 make_tables.py
python3 check_tail_exponents.py 5                 # edge and vertex tail variants at one depth
python3 cone_exponents.py ; python3 cone_exponents.py cube   # tip exponents of the almond cone; vertex exponents of the cube
python3 bem_optical.py 0.1,0.07,0.05,0.04         # Bempp, quasi-uniform cube meshes: field, cross-sections, time
python3 bem_cube_graded.py 8,12,16 2              # Bempp, meshes graded to the edges; add "8 10" for reference quadrature, --wave=60,25,35 for an oblique wave
python3 bem_cube_sym.py 8,12,16 2 [8 10] [--wave=60,25,35] [--check]   # the same using the three mirror planes (item 55); --check compares with the full solve
python3 pec_cube_lightning3d.py 4,5,6 tag=NAME kinds=em q=0.5 gamma=2.0 col_p=3 col_q=0.5 beta=0.5 sym=1 tail=9   # single-knob runs (no vertex tail, no balanced refinement)
python3 diag_cube.py 5 "tail_n=9"                 # residual map by zone, optical theorem
python3 pec_almond_skeleton.py 1.19 --bem         # nose-on case only (see its header)
python3 almond_rcs.py 1.19 6:3 10 --s=0.5         # monostatic RCS, HH and VV, compares with stored Bempp
python3 proto2d_square.py TE ; python3 proto2d_ellipse.py ; python3 proto2d_square_tail.py TE   # 2D prototypes, seconds
python3 proto2d_tapered.py TE                     # 2D: geometric against tapered depths; the ratio sets the attainable accuracy (items 53-54)
```

Costs: section 2, items 43 and 47 (everything before them is obsolete). Almond:
43 s per frequency at 1.19 GHz for all look angles and both polarisations
(6:3 layout).

Files: `mie.py` (reference); `cube_mfs.py` (the cube solver: `Cube`, `run`,
families `balanced` and `dense`); `line_tree.py` (fast line integrals and
point dipoles, mirror images summed inside); `lsq.py` (QR without Q,
semi-normal equations, `factor`); `pec_cube_lightning3d.py` + `line_kernel.py`
(source classes `Sources`, `TailLine`, `VertexTail`, collocation grids, and
the old two-mirror solver `SymSources`, kept as a check); `fig_cost_cube.py`
(cost table, both families); `cube_general.py` (old command line, calls
`cube_mfs`); `bem_env.py` (memory for Bempp's OpenCL driver); `diag_cube.py`
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

Done in the third session of 2026-11-16 (items 33-38), at the author's
request: **T19 Torus** (new example, solver, singular set, cost table
against Bempp, two figures) and a second literature pass with full texts.

Done in the overnight session of 2026-11-17 (items 39-50), at the author's
request: both codes at their best (Bempp: OpenCL and three smaller fixes; MFS:
cube solver 25 times faster), all cost tables re-timed (**T14**), part of
**T16** (Bempp's fast backend; no symmetric, high-order or compressed BEM),
**T18** in part (`cube_mfs.py` is the cube solver; old files kept), the missing
figures of T13, the almond null, the sphere against the wavenumber, the cube
floor (which turned out to be the depth ratio, not the vertices of T15), and a
first draft of the paper (`paper/`).

Open, in the order I would take them:

- [ ] **T13 Write-up, continued** (the author leads). The draft in `paper/`
      is complete in structure. To do: the author's pass over every section;
      decide the title, what goes to supplementary material and the length
      (section 6); acknowledgements and the AI declaration; then copy to
      `sap6/mfs-paper` (ask before pushing: it syncs to Overleaf).
- [ ] **T12 Read what could not be fetched** (the author offered to upload
      PDFs in chat): the list in `LITERATURE.md` section 8. Most urgent:
      Millar 1983 (torus singular set; the draft cites it as "not yet read"),
      Gustafsson 1998 (wording of the remark after the theorem), Kangro 2012,
      Ye et al. 2023, De Smedt and Van Bladel 1986 (compare the exponent
      tables), Kelley et al. 2020 (measurement uncertainty; a TODO in the
      almond section). Karkashadze and Zaridze's MAS review needs a proper
      bibliographic record (year, publisher).
- [ ] **T20 Timo to check** `DERIVATIONS.md` B (edge exponents) and the proof
      of the theorem in the draft (section 5 of the paper).
- [x] **T21 Why does one depth per octave stall?** Answered in items 53-54: a
      fixed ratio has a fixed attainable accuracy (2D model: 8e-5, 6e-6, 1e-9
      for ratios 0.5, 0.6, 0.71); tapered depths help in 2D but not in 3D.
- [ ] **T24 Row-weight exponent for the dense cube family (optional).**
      Tested at J = 5..8 against the reference: for two depths per octave
      beta = 0.7 lowers the error by 1.3-1.5 and beta = 0.9 by 1.5-1.7 at
      every level (0.81e-3 -> 0.60e-3 -> 0.53e-3 at J = 6; 1.20e-4 -> 0.82e-4
      -> 0.71e-4 at J = 8). For one depth per octave 0.5 is right: 0.7 is
      1.2 times better at J = 5, 6 but 2.5 times worse at J = 8, and 0.9 is
      worse throughout. Not adopted (a third of a level; it would mean
      re-running every cube table and the reference).
- [x] **T22 A boundary element code that uses the mirror symmetry**: done for
      the cube with Bempp itself (item 55).
- [ ] **T15 Vertex treatment beyond the diagonal (optional).** No longer
      needed for accuracy (the dense family converges without it).
- [ ] **T17 Non-convex polyhedra (optional, new work).** The skeleton
      theorem and the layouts cover convex bodies only.
- [ ] **T23 Tidy the old cube path (optional)**: `SymSources`, `line_kernel.py`
      and the single-knob `main` of `pec_cube_lightning3d.py` are only kept as
      a cross-check; `check_tail_exponents.py` and `diag_cube.py` still use
      them.

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
- Cube vertex exponents: nu = 0.45417 (known value reproduced), tau = 0.840. Vertex tail: 0.31 / 0.133 / 0.049% -> 0.20 / 0.086 / 0.029% (J = 5, 6, 7).
- Graded Bempp, grading scan at 2592 dofs: g = 1, 2, 2.5 (q 6/8), 3: 0.63%, 0.15%, 0.20%, 4.8%. Quadrature 4/6, 6/8, 8/10 at 4608 dofs: 0.104%, 0.060%, 0.057%.
- Graded Bempp reference (q 8/10), 4608 / 7200 / 10368 dofs: 0.058%, 0.032%, 0.022% from MFS J = 9; extrapolated 0.016%.
- Cube cost table re-timed in one container: MFS J = 3..8 22-284 s for 1.27%-0.011%; Bempp uniform 63-832 s for 0.68%-0.19%; graded 32-161 s for 0.44%-0.105%.
- Cube, oblique wave, J = 5, 6, 7: 86, 176, 370 s for 0.27%, 0.112%, 0.055% against graded Bempp; standard wave reproduced to 2e-15.
- Almond against the Austin benchmark, 3.5-10.25 GHz: 0.04-0.14 dB (HH), 0.01-0.02 dB (VV) from the reference simulation; about 1 dB from measurements, as the reference is.
- Almond 10.25 GHz, layouts 6:3 and 8:5: 0.0006 dB mean difference, 0.024 dB in the -66.7 dBsm null where the reference simulation differs by 4.3 dB.
- Two jobs at once: QR 20 times slower, timings of that batch discarded.

- Torus, known exterior solution (3 dipoles in the tube): 1.8e-3, 3.1e-5, 1.2e-7, 4.2e-10 for 1440 to 13824 unknowns; mode-by-mode solver equals the full matrix to 8e-15.
- Torus R = 1, a = 0.5, wave along the axis: residual 2.5e-2 to 1.0e-10 for 1440 to 21120 unknowns; sigma_ext 11.75257471; optical mismatch 9e-14 at the finest.
- Torus, fitted singularity distance, a/R = 0.8: 0.395-0.401 (axis), 0.39-0.40 (in plane), 0.391-0.400 (oblique), 0.398 (k = 4 pi); predicted R - c = 0.400. Centre moved to 1.1: 0.49 (0.50). a/R = 0.5, 0.25: no clean fit.
- Torus, |alpha| against a_s/a at 19200 unknowns: a/R = 0.8: 1.9, 3.6, 41, 7e4 at 0.55, 0.50, 0.45, 0.40; a/R = 0.5: 13, 49, 460 at 0.30, 0.25, 0.20; a/R = 0.25: 24, 770 at 0.15, 0.10.
- Torus, off-centre source circles with 2 collocation rings per source ring: residual 1 to 1e4 (silent); with 4: 6.5e-8.
- Torus cost: MFS 0.28-10.1 s for 5e-3 to 7e-12; full matrix 25.5 s for 1.1e-4; Bempp 9-334 s for 10% to 0.82%, extrapolated sigma_ext 11.7533 against 11.7526.
- Torus, oblique wave: MFS residual 2.6e-10 (21120 unknowns, 1.5 s solve); Bempp 4704 dofs 2.0% away.

- Bempp backends, sphere 4824 dofs: Numba double 108 s, Numba single 92 s, OpenCL double 10 s, OpenCL single 8.3 s; error 7.58e-3 in all four.
- Bempp, one 714-dof assembly: 1.9 s the first time, 0.27 s afterwards; right-hand side closure recompiled at every solve: about 3 s (fixed).
- PoCL memory cap: 1, 2 or 4 GB per buffer depending on machine and load (one NUMA node counted); with `bem_env.py` 4 GB.
- Cube line quadrature (`line_tree.py`) against the adaptive rule at twice the nodes: 3.4e-10 worst (J = 4), 4.5e-10 (J = 5); field unchanged to 2e-12.
- Cube, three mirrors instead of two: fields equal to 2e-12; J = 6 from 81 s to 7.5 s; then images inside the kernel and semi-normal equations: 3.8 s, field equal to 1e-15.
- Cube, cond of a scaled class matrix: 3.8e4 (J = 6), 7.7e5 (J = 8), 8.2e5 (J = 9); almond classes 5e4 to 3e5.
- Cube, two collocation nodes per panel instead of three (one depth per octave): 85% error at J = 5, 0.2% at J = 6.
- Cube, one depth per octave, J = 10, 11, 12: 0.0125%, 0.0145%, 0.018% against dense J = 10 (stalled); variants at J = 10 differ by 6e-5.
- Cube, two depths per octave at J = 9: residual 2.2e-4 (was 7.6e-4), optical 1.3e-5 (was 4.0e-5); two collocation grids agree to 9.7e-7; three nodes per panel: garbage (optical 0.5).
- Cube, dense family J = 3..10: change to next 9.0e-3 .. 3.5e-5 (ratio 2.4-2.9); sigma_ext 2.927961 at J = 10; graded Bempp extrapolated 2.927864.
- Cost tables re-timed on one machine (items 43, 47): sphere, torus, almond, cube standard and oblique.
- Cube boundary residual at J = 8 by distance to the nearest edge: equal shares of the weighted squared residual per decade (0.1-0.3); rms 0.43 at 1e-5..1e-4, 0.02 at 1e-3..1e-2.
- Residual near an edge, J = 6: without tail 1.6 below 1e-3; with tail 1.2e-2 (1e-4..1e-3), 0.13 (1e-5..1e-4).
- Almond null, 10.25 GHz HH: MFS minimum at 35.66 deg, -69.7 dBsm; size +-0.01%: 0.4 dB at 35.5 deg; shift by -0.10 deg: rms difference from the reference 1.45 -> 0.64 dB.
- Sphere against wavenumber, R = 0.5: N for 1e-6 is 224, 473, 797, 1628 at k = 2, 4, 6, 8 pi (21-28 unknowns per square wavelength from 4 pi).
- Sphere MFS, QR against SVD solver at N = 1200: 2.35e-14 against 2.14e-14; 3.4 s against 15 s.
- Bempp with mirror symmetry on the cube: equals the full solve to 2.5e-5 (4/6) and 1.6e-6 (8/10); 14112 dofs in 5.6 s (4/6), 22 s (8/10); 28800 dofs, 8/10: 69 s, 1.00e-4 from the MFS reference.
- 2D square, free solves to depth 1e-9, TE: ratio 0.5, 0.6, 0.71, 0.8: 8.1e-5, 5.8e-6, 1.4e-9, 1e-11; tapered 30 / 40 depths: 9.5e-8 / 1.0e-9.
- Cube with tapered depths, J = 9: 10, 13, 16 depths: 4.3e-4, 1.27e-4, 5.3e-5 (uniform: 1.25e-4 with 10, 3.5e-5 with 19).
- Cube, row-weight exponent 0.5 / 0.7 / 0.9 at J = 8: one depth per octave 1.9e-4 / 4.7e-4 / 2.4e-3; two depths 1.20e-4 / 0.82e-4 / 0.71e-4.
- Almond with the semi-normal equations: fields equal to 1e-15; RCS against Bempp 0.023 / 0.029 dB as before; 43 s per frequency.
- Almond against the Austin benchmark, re-run with the new solver: the four frequencies agree with the stored curves to 1e-13 dB; 66 s per frequency (361 angles, HH and VV).
- Almond 1.19 GHz nose-on, variants of layout 8:5: spacing 0.5: 8.3e-6; finer collocation: 3.3e-6; two rows per octave: 2.7e-6.
- Cube dense J = 9 variants: modes 6.3e-7, tail length 6.3e-7, second vertex exponent 1.7e-7, other collocation grid 9.7e-7; beta 0.3: 5.9e-5; beta 0.7: 2.1e-5.

- Cube base grid 6 / 0.15 instead of 8 / 0.1 from J = 7: dense family fields equal to 1.4e-6..2.2e-6 (J = 7..10), 25% faster; one depth J = 9: 1.50e-4 against 1.25e-4. Adopted.
- Cube, fewer modes (n0 - 2, n0 - 4) at J = 7, 8: error larger at equal J; not adopted.
- Cube base grid 4 / 0.25 at J = 7, 8: errors 6e-3, 2.4e-2 with optical mismatch 1.1e-3, 1.2e-5 (the mismatch does not flag the second).
- Cube one depth per octave with base 6 / 0.15, J = 10, 11, 12: 2.85e-4, 4.6e-4, 9.0e-4 (optical 3.0e-5, 3.1e-5, 1.2e-4); with 8 / 0.1 they were 1.3e-4..1.8e-4.
- Cube tables re-run with the adopted grid (item 56): one depth J = 7, 8, 9: 5.2, 7.7, 10.5 s for 3.6e-4, 1.6e-4, 1.5e-4; two depths J = 8, 9, 10: 24, 36, 50 s for 1.2e-4, 3.6e-5, reference.

## 6. Write-up session

**State on 2026-11-17: a complete first draft exists in `paper/`** (in this
repository, not in `sap6/mfs-paper`; nothing was pushed there). `main.pdf` is
committed. The author asked for it ("make a start on a first draft"), in the
voice of the thesis and with the academic-humanizer rules (item 50).

- Structure: introduction; formulation; sphere; torus; why an inner copy
  fails (theorem); cube; almond; discussion and limits; conclusion;
  appendices A-E (edge exponents, cone exponents, torus argument, line
  quadrature, Bempp settings). 13 figures and tables in the body.
- Open in the draft (marked `%% TODO`): title; acknowledgements; AI
  declaration; measurement uncertainty from Kelley et al.; Millar 1983;
  comparison of the exponent tables with De Smedt and Van Bladel.
- Length: 42 pages in the 12 pt preprint format (`main.pdf`), **26 pages in
  the journal's own layout** (`main_3p.pdf`, rebuilt after item 56). The author's target is 25, which the journal layout
  meets. If the preprint format must come down too, in the order I would
  try: (1) the three cost tables of sphere, torus and almond duplicate
  their figures; `validation/fig_cost_all.pdf` is one four-panel figure
  that can replace the four cost figures and those three tables;
  (2) Appendix E (Bempp settings) and the oblique-wave table to
  supplementary material; (3) the sphere-against-wavenumber paragraph to
  two sentences.
- The comparison with Bempp is now fair to both sides (items 39-41, 43, 47,
  55). For the cube the draft reports Bempp both without and with the mirror
  symmetry and compares on the equal footing first; the remaining
  qualification is the number of unknowns.
- Decisions still open for the authors: title; how much of the Bempp cost
  comparison stays in the main text; whether the one-depth-per-octave
  layout is presented as the default with the dense one as the high-accuracy
  option (as drafted) or the other way round.

The rest of this section is the state before the draft and still applies
where it is not superseded by items 39-50.


- Attach both `sap6/MFS` and `sap6/mfs-paper`. Read this file,
  `PAPER_OUTLINE.md` and the manuscript; do not re-run experiments.
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
  - *Torus (third session)*: say that MAS has been applied to tori before
    and that tori are standard MFS test bodies in acoustics; the content here
    is the Maxwell convergence study and the admissible source radius. The
    focal-ring location is argued and measured, not proved; no closed-form
    reference exists (Bempp confirms to 1e-4 by extrapolation only).
  - *Skeleton theorem*: cite Gustafsson 1998 as the analogue for Newtonian
    potentials; "to our knowledge not stated for scattering" is the most
    that can be said.
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

- **Read first (2026-11-17).** (1) Bempp needs `pyopencl` and
  `pocl-binary-distribution` (in `CLAUDE.md` now); check that a Bempp script
  prints or uses `opencl`, not `numba`. (2) The laptop is no longer needed
  for Bempp: with OpenCL the container is faster than the laptop was with
  Numba, and the largest mesh that fits (16384 dofs) takes under 3 min.
  (3) All timings of the paper come from one run of `retime_all.sh` plus
  `fig_cost_cube.py run` and `cube_oblique_table.py run` on one machine
  (2026-11-17, 02:25-03:25 UTC); if any of them is repeated on another
  machine, repeat all. (4) The worker was replaced once in this session
  (02:20 UTC); installed packages and uncommitted files survived, running
  jobs did not, and PoCL's memory cap changed with the machine (item 40).
  (5) Memory: the kernel killed one job (eight classes of the dense cube
  family at J = 8, 11.5 GB); `cube_mfs.run` now splits classes into groups.
  Keep single matrices below about 6 GB. (6) LaTeX is installed with the
  apt line in `CLAUDE.md` (2 min); `paper/main.pdf` builds in 20 s.
  (7) Network: api.crossref.org and export.arxiv.org work over HTTPS with
  10-60 s between requests (`paper/check_bib.py`, `paper/check_arxiv.py`).
- **Where to run.** The cloud container (4 cores, 15 GB) is the primary copy
  and is fast enough for everything except large Bempp meshes. Rough times
  (second container): cube J = 6: 95 s, J = 8: 284 s, J = 9: 419 s; cube, any
  wave, J = 6: 176 s; almond RCS, all angles and both polarisations: 3 min
  per frequency; Bempp 4608 dofs: 160 s, 13k dofs: 14 min. Setup on a fresh
  container: the three lines in `CLAUDE.md` (2-3 min).
- **The author's laptop** (`flowx13`, 16 cores, 15 GB, Linux) can be attached
  through Remote Control (folder `Projects`). Use it for Bempp with more
  than about 5k dofs, where it is 3-4 times as fast (2592 dofs: 22 s against
  75 s); its dense QR is slower than the container's, so keep the MFS here.
  Two checkouts exist there:
  - `~/Projects/MFS-ccr`: a git worktree of the working branch (local branch
    `laptop-bem`), made in the second session. **Use this one.** Before a
    run: `git pull origin ccr-116c520c-fts35f` in it. Run as
    `/home/sap/Projects/MFS/.venv/bin/python3 -W ignore script.py ...` from
    `~/Projects/MFS-ccr/validation`. To return results: commit there and
    `git push origin laptop-bem:ccr-116c520c-fts35f`, then pull in the
    container. Remove it with `git worktree remove ../MFS-ccr` (run in
    `~/Projects/MFS`) when no longer wanted.
  - `~/Projects/MFS`: the old clone, detached at efac0ae with hand-made
    edits to `pec_almond_skeleton.py` and untracked Bempp files. Not touched
    in either session; it only lends its `.venv` (numpy, scipy, numba,
    bempp-cl; gmsh's Python module is missing, so `bempp.shapes` does not
    work there, but `bem_cube_graded.py` builds its own meshes). Ask the
    author before changing it.
  Background jobs on the laptop do not notify the session: read their log
  file to check, and never assume they are finished.
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
- **Network, third session.** The publisher refusals are the publishers' own
  bot walls (Cloudflare challenge; Elsevier also reserves text mining), not
  the session's network policy: sciencedirect, epubs.siam.org, mdpi, dl.acm
  return 403 to any automated client. web.archive.org and
  perso.univ-rennes1.fr reset the connection. Open and useful: arxiv.org,
  lae.tsu.ge (the Tbilisi MAS papers as PDFs), core.ac.uk, link.springer.com
  (landing pages). The author offered to supply PDFs: ask for them in
  `~/Projects/MFS-ccr/papers/` on the laptop.
- **Worker restarts** happened again in the third session (three in a row):
  the scratch folder was lost, the repository survived. Downloaded papers
  were re-fetched; nothing else was lost because results were in
  `validation/`.
- **One job at a time, really** (item 31). Chain jobs in one background
  command; wait for a process with `tail --pid=PID -f /dev/null`.
- **Two containers, two speeds.** The second container of 16 November was
  about 20% slower than the first for both the MFS and Bempp. Never mix
  timings from different sessions in one table; `fig_cost_cube.py` data and
  all cube Bempp files were re-timed together.
- **Scratch outputs** named `*_t_*.npy/json/npz` and `almond_rcs_*.npz` are
  ignored by git on purpose (logs too). The scratchpad folder outside the
  repository is lost at every restart.
- **The manuscript** (`sap6/mfs-paper`) can be attached read-only with
  `add_repo`; it then sits at `/home/user/mfs-paper`. Rule 8 of `CLAUDE.md`
  stands: no push there without asking.
- **Figures** follow one style (colour plus marker plus dash per method, so
  they print in greyscale); see `STYLE` in the `fig_*.py` scripts.
