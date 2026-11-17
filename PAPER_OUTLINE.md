# Paper outline and critique (task T7), 2026-11-16

> **2026-11-17: a first draft now exists in `paper/`** (see `NOTES.md`,
> section 6 and items 39-50). This outline is kept as the plan it was written
> against. Three things in it are superseded: every cost figure (both codes
> were re-timed at their best), the cube accuracy ("0.011%, floor"; now
> 0.0035% with two source depths per octave, reference good to 2.5e-5,
> Bempp agreeing to 0.014%), and the almond null (a 0.1 degree shift).

For the author. Bullets and pointers only; the prose is yours. Every number
below is in `NOTES.md` with the script that prints it. Derivations are in
`DERIVATIONS.md`, references (with how each was checked) in `LITERATURE.md`.
The manuscript read for this was `sap6/mfs-paper` at commit 4f8731f (the
Overleaf import, i.e. still the thesis text). Nothing was pushed there.

## 1. What the paper can now claim, in one paragraph of bullets

- The MFS for 3D Maxwell (PEC, tangential-trace least squares) converges like
  exp(-c sqrt(N)) on the sphere down to 1e-13 against Mie, with the same
  stability picture as Barnett and Betcke's Helmholtz analysis.
- Placing the sources on an inner copy of the surface fails for the cube and
  the almond, and the reason is the same in both: the continued field has its
  singularities on the *skeleton* (bisector planes of the edges; the focal
  curve and mid-plane of the almond).
- For a convex PEC polyhedron this is a theorem with a ten-line proof
  (image principle): the scattered field is a layer of electric and magnetic
  current on the interior medial surface.
- Sources on the skeleton, clustered geometrically towards the singular set,
  plus a "singular tail" whose exponents are *computed* (edge: 2/3, 5/3, 4/3,
  7/3 from the Meixner expansion; cube vertex: 0.454 and 0.840 + 2; almond
  tip: 0.2996 + 2 from the elliptic cone), give convergence on both shapes:
  cube to 0.011% (self-convergence) and 0.02% against an independent
  graded-mesh BEM; almond to 0.006% (nose-on) and 0.01-0.14 dB against a
  published fine-mesh reference at four frequencies, inside the 1 dB scatter
  of the measurements.
- Cost, stated without spin: sphere and almond 30 to several hundred times
  faster than Bempp at far higher accuracy; cube 10-15 times faster than Bempp on a
  quasi-uniform mesh at 0.2%, about level with Bempp on an edge-graded mesh
  at 0.15-0.2% and ahead below 0.1%; for a general incident wave the cube
  costs 1.6-1.9 times the symmetric case and graded Bempp is then faster at
  0.1-0.4%. The cube's selling point is the rate, not the cost at 0.2%.

## 2. Proposed sections

1. **Introduction.** MFS / MAS / discrete sources in electromagnetics
   (Kupradze; Kaklamani and Anastassiu; Doicu, Eremin, Wriedt); the known
   rule that sources must enclose the singularities of the continued field
   (Barnett and Betcke 2008; Leviatan 1990; Fikioris 2006; Tsitsas et al.
   2018); lightning solvers as the 2D model for clustered sources (Gopal and
   Trefethen; Ginn and Trefethen for Helmholtz) and the open 3D question
   (quote Kong and Shen 2026, Brubeck and Trefethen 2022); recent MFS for
   Maxwell on smooth surfaces (Weed et al. 2026). Contribution list = the
   bullets of section 1. State that the paper grew out of the MSci thesis.
2. **Formulation.** Dyadic fundamental solution (electric and magnetic
   dipoles); PEC condition as two tangential components per collocation
   point; weighted least squares by QR; error measures: Mie, boundary
   residual, optical theorem (not automatic for collocation), a second
   solver. One paragraph on 2D Helmholtz as the model, not a section.
3. **Sphere.** Convergence in N and source radius R against Mie; coefficient
   norm and loss of digits; best R; sqrt(N) rate. Cost against Bempp.
   Figures: `pec_sphere_convergence`, `fig_cost_sphere`.
4. **Why inner surfaces fail.** Cube: stall at 21-30%. Almond: confocal /
   shrunk surface unstable or 22%. Explanation by continuation; the skeleton
   theorem (DERIVATIONS A) and its three consequences.
5. **Edge sources and the singular tail (cube).** Line dipoles on the
   bisector planes at depths d_j = d_0 2^-j; modes in a log-stretched
   variable; electric and magnetic; weights; the tail and its exponents
   (DERIVATIONS B, table); balanced refinement. Results: convergence table,
   optical theorem, graded-BEM reference, vertex tail, oblique incidence by
   parity classes. Figure: `fig_cost_cube`; a sketch of the layout is still needed.
6. **Skeleton sources and the tip tail (almond).** Rows on the mid-plane
   towards the focal curve; parity classes; tip cone exponents (DERIVATIONS
   C, table) and `cone_exponents.py`; results at 1.19 GHz against Bempp, and
   at 3.5, 5.125, 7, 10.25 GHz against the Austin benchmark (measured and
   simulated). Figures: `fig_almond_benchmark` (main), `fig_cost_almond`,
   `fig_almond_rcs` (could move to supplementary material).
7. **Discussion and limits.** Section 4 below, in the paper's own words.
8. **Conclusion.** Three sentences; no "robust and efficient for all shapes".
   Appendix: exponents of the cube vertex (DERIVATIONS D), parameters.

### Addition of the third session: a torus section (recommended, author decides)

Place: after the sphere, before "why inner surfaces fail" (new section 3b).
About one page and two figures (`fig_cost_torus`, `fig_torus_sources`).

- What it shows: the thesis's own layout (sources on an inner copy of the
  surface) converges on a non-convex, ring-shaped conductor: 1e-11 in 10 s,
  any incident wave, against 0.8% in 334 s for Bempp.
- What it adds to the argument: the inner copy may not be arbitrary. The
  continued field is singular on the centre circle and on the focal ring
  rho = sqrt(R^2 - a^2), so a_s > R - sqrt(R^2 - a^2) (DERIVATIONS E). The
  coefficient norm stays near 1 above that radius and grows by orders of
  magnitude below it: the picture the thesis drew for the sphere (error and
  |alpha| against R), now with a threshold that a formula predicts.
- It makes the sequence of shapes a single story: sphere (singular point at
  the centre, any radius works), torus (singular rings inside the tube, the
  radius is bounded), almond (singular curve within 6% of the rim, no inner
  copy works), cube (singular set reaches the surface at the edges, tails).
- Limits to state: no closed-form reference; Bempp confirms 0.8% directly
  and 7e-5 by extrapolating the cross-section; the mode-by-mode solver uses
  the rotational symmetry of the body; MAS has been applied to tori before
  (Karkashadze and Zaridze) and acoustic MFS uses tori as test bodies (Cai,
  Bagge, Martinsson 2026); the focal-ring statement is argued and measured,
  not proved (read Millar 1983 first).

New references to work into sections 1 and 4 (LITERATURE section 7): the
Gopal-Trefethen question about 3D (quote it in the introduction); Hochman,
Leviatan, White 2013; Kangro 2012; Bogdanov, Karkashadze, Zaridze on
rounding edges; Gustafsson 1998 next to the skeleton theorem.

## 3. Figures and tables that exist

| Item | Script | Status |
|---|---|---|
| Sphere convergence in N, R | `pec_sphere_convergence.py`, `plot_pec_sphere.py` | done |
| Sphere error against time | `fig_cost_sphere.py` | done |
| Cube error against time, five methods | `fig_cost_cube.py` | done, regenerated 16 Nov |
| Almond error against time | `fig_cost_almond.py` | done |
| Almond RCS against Bempp, 3 frequencies | `fig_almond_rcs.py` | done |
| Almond RCS against measurement and reference simulation, 4 frequencies | `fig_almond_benchmark.py` | done, new |
| Exponent tables (edge, tip cone, cube vertex) | `cone_exponents.py`, DERIVATIONS | done, new |
| Torus error against time; residual and coefficient norm against source radius | `fig_cost_torus.py`, `fig_torus_sources.py` | done, new (third session) |
| Source-layout sketches (cube cross-section at an edge; almond mid-plane) | `explainer/make_diagrams.py` draws them from the solvers' real coordinates (SVG, for the author's explainer) | drafts exist; redraw for print |
| Residual map near an edge with and without tail | `diag_cube.py` prints it | no figure yet |

## 4. What must be said plainly (limits)

- Below 0.02% the cube reference is the MFS itself. Independent evidence:
  graded Bempp (reference quadrature) at 4608, 7200, 10368 unknowns is
  0.058%, 0.032%, 0.022% from it, 0.016% when extrapolated, and its
  cross-section rises towards the MFS value; the optical theorem mismatch
  is 4e-5.
- The cube and nose-on almond timings use two mirror symmetries; Bempp uses
  none. General incidence is implemented for both shapes (cube: 86, 176,
  370 s for 0.27%, 0.11%, 0.055%).
- Bempp ran with dense assembly and LU. Assembly is 98% of its time,
  so a better linear solver would not change the comparison; matrix
  compression (FMM, H-matrices) or higher-order elements would, and were not
  tried. Graded meshes were tried and are reported.
- All three shapes are at most about ten wavelengths long (almond at
  10.25 GHz: 8.6). Nothing is claimed about larger problems; the dense QR
  scales like N^3.
- The source layouts are built by hand for each shape from its skeleton.
  There is no automatic construction for a general body.
- The polyhedron theorem is for convex bodies.
- Almond: agreement with the published reference simulation is 0.01-0.14 dB
  on average, but single deep nulls differ by up to 4.3 dB (HH, 10.25 GHz,
  35.5 degrees, a null at -66.7 dBsm). The MFS is self-converged to 0.02 dB
  there (NOTES item 30); which simulation is right in the null is open.

## 5. Critique of the current manuscript (referee style)

Major:

1. **The Maxwell results are not scattering problems.** Sections on the
   sphere, cube and almond impose E = p exp(ik d.x) on the surface and solve
   the *interior* problem, whose exact solution is the plane wave itself
   (the captions say "the exact modulus inside ... is 1"). Errors of 1e-9 to
   1e-12 there say nothing about scattering by these shapes. The titles
   "Plane Wave Scattering by a Cube / NASA Almond" are wrong as they stand.
2. **The one exterior case is under-determined and unvalidated.** "Exterior
   scattering on a perfect electric conductor" has no reference solution
   (Mie is available) and, per the notebook, 2M equations for 3N unknowns
   with M = N.
3. **Formulation.** The system enforces all three components of E (3M rows).
   For a PEC only the tangential components are prescribed; the normal
   component of the scattered field is not data. The new solvers use the
   tangential trace; the section must be rewritten around it.
4. **Convergence claims.** "Exponential in N" is not what the data show in
   3D (sqrt(N)); the coefficient-norm discussion should be tied to the
   position of the continuation's singularities, with Barnett and Betcke as
   the model.
5. **Solver paragraph.** The text credits GMRES and "exploiting the 3x3
   block structure"; the code solves a dense least-squares problem by QR.
   The complexity remark O(27 M N^2) is fine; the rest should go.
6. **Related work.** Twelve references, none on MAS, discrete sources,
   generalized multipoles, lightning solvers, or singular fields at edges
   and tips. `LITERATURE.md` section 5 has about fifty checked entries.
7. **Conclusion over-claims** ("robust and efficient approach", "versatility
   in handling diverse geometries") on the evidence of interior problems.

Minor:

- Thesis wording throughout ("this project", "foundational steps",
  "chapter"); abstract describes a project, not results.
- Helmholtz part is 264 lines for material that is known (cut to a
  paragraph or a short model-problem subsection).
- Journal name in `main.tex` is "Computational Physics" (should be Journal
  of Computational Physics); keywords fine.
- "Derivation of the Fundamental Solution" (110 lines) can be three lines
  and a citation; the final dyadic formula is correct (checked).
- The appendix mapping figures to notebooks should become a code and data
  availability statement pointing at `validation/`.
- Colour-only figures (3D surface plots with jet-like maps): most can go;
  the new figures use colour plus marker plus dash.
- Bibliography: in `kirsch` the second author is written "Hettlich F."
  (needs "Hettlich, F."); `colton1998inverse` lists Kress twice; three
  entries give the DOI as a URL (`doi = {https://doi.org/...}`, `url =
  {http://dx.doi.org/...}`) instead of a bare DOI.

## 6. Questions a referee is likely to ask, and where the answer is

| Question | Answer |
|---|---|
| Is the cube reference independent? | Section 4, first bullet; NOTES section 2, item 29 |
| Are the tail exponents tuned? | No: DERIVATIONS B and C; fitted sets were replaced by computed ones with no loss |
| Is the comparison with BEM fair? | Graded meshes, timing split, symmetry statement: NOTES items 26-28 |
| Does it work for other incidences? | Cube: `cube_general.py`; almond: all look angles in the benchmark figure |
| Is the almond compared with measurements? | Yes, public Austin benchmark data, four frequencies |
| What about the vertices of the cube? | Vertex tail with computed exponents (0.454, 0.840), 12 unknowns, error down by a third to 40%: NOTES item 25 |
| Novelty against lightning and MAS work? | LITERATURE section 1 |

## 7. Decisions for the authors

- Title: the working title ("stability and convergence of the MFS for 3D
  Maxwell scattering") undersells the edge and tip treatment; consider
  naming it.
- Whether the cube cost comparison leads with the quasi-uniform or the
  graded BEM curve. Recommendation: show both; Timo will ask for the graded
  one in any case.
- Whether to keep the 1.19 GHz Bempp comparison of the almond now that the
  benchmark comparison exists (recommendation: keep one cost figure, move
  the three-frequency Bempp figure to supplementary material).
- Third-party data: `validation/austin_almond/` holds unmodified benchmark
  files under their own licence (CC BY-SA 4.0 plus the suite's notice). If
  you prefer not to redistribute them, delete the folder and have
  `fig_almond_benchmark.py` point at a local download instead.
