# Referee-style reading of the draft, 18 November 2026

A fresh agent read the whole paper (text, tables, figures) as a JCP referee
would. Its findings, condensed; the status of each is tracked in
`NOTES.md` item 80.

## A. Wrong or unjustified

1. **Theorem 1 does not locate singularities.** The proof itself says the
   reflected field "is defined in all of Omega", so the continuation through
   each face is regular in the whole body. The only singularities are the
   edges and vertices, around which the continuation is multi-valued. The
   skeleton is one admissible branch cut; the same Stratton-Chu argument
   works for any walls that span the edges with one cell per face.
   Affected: abstract "which shows where the sources must go"; skeleton.tex
   "In this section we locate these singularities"; cube.tex "the sources
   belong on the twelve bisector planes"; conclusion.tex "for the cube it is
   the skeleton" and "The sources lie on the singular set" (contradicts
   "enclose" and the base points on the inner cube). Correct claim: sources
   must reach the edges; the medial surface is a sufficient, symmetric
   choice.
2. method.tex "corrected semi-normal equations [Bjorck 1987], which cost
   about half as much as QR": in Bjorck's method R comes from a QR
   factorisation. Half the cost means Cholesky of A^H A (`lsq.py`): normal
   equations with a correction step. Name it so and justify it: at
   condition 1e6, kappa^2 u is about 1e-4, above the claimed errors 2e-5 and
   6e-6. Give a QR cross-check at J >= 10.
3. smooth.tex "The sequence ... increases to its fixed point ... from any
   starting point inside the disc": it decreases for starting points
   between R - sqrt(R^2 - a^2) and a. Write "converges monotonically".
4. smooth.tex "stays of order one while a_s satisfies": Fig. 2 (right)
   shows a knee at the dotted line only for a/R = 0.8. For 0.25 and 0.5 the
   norm rises smoothly through the bound, and the residual minimum for
   a/R = 0.5 is at a_s/a = 0.25, below the bound 0.268. State the outcome
   for 0.25 and 0.5.
5. cube.tex "The vertex points are continued in the same way ... nu + 2
   and tau + 2": not the same way (the edge rule would give nu). The
   scaling argument is in the appendix "From exponents to weights" and is
   not referenced here; it assumes without argument that the electric-type
   term is carried by electric dipoles only.
6. cube.tex "This is the root-exponential convergence of the lightning
   solvers": the data say error ~ d_min^1.3 (2.5 per level; without tails
   d_min^0.6), and N_u ~ J^2 by construction. In lightning solvers the rate
   comes from tapered clustering, here from how unknowns are counted at a
   fixed ratio. State the d_min law and explain 1.3.
7. cube.tex "1.7e-5 at J = 10" (general wave) is measured against J = 11,
   whose own error is about 7e-6, 40% of the figure. cost.tex "continues to
   1.7e-5" inherits this.
8. cost.tex "The MFS does not have it" (near-field loss) is shown only for
   the sphere. On the cube Fig. 4 shows a residual of 1e-2 within 1e-3 of an
   edge. Show cube error against distance, or restrict the claim.

## B. Inconsistencies

- cube.tex "shows this (Fig. 5)", "at the levels 10 to 12 its error stays
  between": the one-depth curve in Fig. 5 ends at J = 9; levels 10-12 are in
  no figure or table.
- cube.tex "1.3 to 1.5 times ... up to J = 8": J = 8 is 41 s / 26 s = 1.6.
  "to within 20%" is 21% at J = 7.
- abstract "the reference uses 400,000"; almond.tex "397,674": write
  "about".
- cost.tex "between 6e-5 and 7e-12 in ten seconds or less": the sphere
  reaches 2.3e-14 in 5.0 s.
- Table 4 "Unknowns": per parity class for the MFS cube, whole mesh for
  Bempp; Fig. 5 uses N_u = 44160. One convention.
- "nose-on", "tip-on", "towards the tip"; Fig. 8 legend "ARCHIE-AIM" never
  named in the text.
- Appendix tables are numbered A.5 and C.6.

## Citations to verify

- kangro2012 (four-page proceedings) "proves" exp(-c sqrt N).
- desmedt1986 "0.45410 and 0.84108" for the cube corner; that paper treats
  cones.
- kelley2020: 397,674 unknowns, mesh differences 0.07/0.02/0.5 dB, the 1 and
  2 dB bands.
- Quotations of Gopal-Trefethen and Kong-Shen; aberra2002 and millar1983
  statements; ye2023 "the feature is smoothed".
- introduction "to smooth bodies or to bodies whose edges were rounded":
  refs.bib holds uncited tabatadze2025, eisler1989, shifman1997 (edges).

## C. Notation

- J: cube level, almond "n_mid:J", currents J, J_s, J_rho.
- The coordinate along the edge is t in 5.1 and Appendix B, s in 5.2; t is
  also the tangent vector and the superscript "total"; s is also the tangent
  index in (7).
- E is the scattered field in (2), the total field in (15)-(16) and
  Appendix A. H first appears in the proof, undefined.
- M boundary points, calligraphic M skeleton; a is torus radius, ellipsoid
  semi-axis and a_l(s).
- q defined and never used; "rim" and "coefficient norm" undefined.

## D. Clarity, over- and under-claims

- cube.tex heading "Sources on the skeleton": the base points are on an
  inner cube; say why that is allowed.
- cube.tex "w_i = A_i rho_i": say whether the boundary residual uses these
  weights; it is defined in 2.4 and never reported for the cube.
- "The two-dimensional analogue" gives no parameters. "Measured as an
  algebraic order" duplicates 7.3.
- skeleton.tex "With electric dipoles alone ... 5%": no level, no data.
- "over the range that we could test" three times. Under-claim: Fig. 9b
  shows the MFS without symmetry ahead of Bempp with symmetry.
- almond.tex "of the size of the difference between the two finest
  meshes": that comparison exists only at 7 GHz.

## E. Missing

- The standard wave travels along the torus axis and excites only azimuthal
  orders +-1; say so and report a general direction.
- Evaluation points of (8) for torus and almond; n_s, n_phi; M for torus and
  almond; the "smaller multiple" for almond boundary points; positions of
  the n_mid rows; waves and wavenumbers of the torus tests.
- Quadrature for sigma_sca in (9).
- Why the first magnetic exponent of the cube corner has only three digits.
- Whether "Time" in Table 2 covers both classes.
- Bempp meshes for sphere, torus and almond.

Consistent wherever they appear: 2e-14, 7e-12, 2e-5, 8e-5, 0.14 dB, 4193,
21-30%, the fit, tail share 5%, and all cross-references.

# Second round, 18 November, after the second timing run

## Number check (fresh agent, from the .npz files)

Every table entry and quoted error matches the data. Problems are
comparative claims that rest on log-log interpolation between sparse runs
("staircase" = fastest run with error at most the target):

1. General wave, "below that level the MFS is faster": no MFS run lies
   between J = 8 (1.62e-4, 37.8 s) and J = 9 (6.0e-5, 74.3 s); Bempp has
   1.44e-4 in 39.3 s and 1.08e-4 in 64.7 s. The crossing at 1.39e-4 exists
   only by interpolation.
2. Standard wave: the crossing at 1.46e-4 (24 s) holds under both
   conventions; below it the MFS leads by only 1.12 at 1.15e-4
   (interpolated) and by staircase Bempp wins between 1.40e-4 and 1.15e-4.
   The lead is robust only below 1.1e-4 (1.6 to 1.9 at 8e-5).
3. Three wavelengths, "level near 6e-4": at 6.0e-4 Bempp is still 2.1
   times faster; the two are level only at 5.2 to 5.7e-4 (Bempp's last two
   runs). "Four times at 1e-3" is 3.9 interpolated, 2.7 staircase.
4. Summary table: two Bempp rows are not the best runs at their cost.
5. Torus "2e-7" is fortuitous: the previous pair of meshes gives 1.2e-6.
   Say "about 1e-6".
Convention-dependent: "about ten at 1e-3" (9.9 / 7.1); two wavelengths
level at 3.0e-4 / 2.3e-4; one depth per octave level at 3.1e-4 / 2.6e-4;
"doubles with each wavelength" holds by interpolation only.
Minor: fit within 4.2% with the printed coefficients; "a third of the
time" is 0.44 to 0.30; Bempp's order on the finer meshes 0.99, 0.99, 0.77;
the J = 10 references at two and three wavelengths have about 15% of the
J = 8 error.

## Consistency read (fresh agent)

Leftovers of "the theorem locates the singularities": introduction
(contribution 2, roadmap), heading of Section 4, "must be placed". The
rule "enclose the singularities" against sources that reach them on cube
and almond; "surfaces that leave the focal curve outside" against rows
that approach it from inside; "the inner cube encloses the part of the
skeleton ... as the rule for smooth bodies requires". Captions say
"mirror symmetry" for the cube where the rotation is used too; "does not
change its error" against "differs by 3e-5". Seams: repeated sentence in
the cube comparison, the rate stated three times in one paragraph, a few
dangling "it" and "this". LaTeX: appendix table links point to Table 1;
Table 4 and eq. (4) are too wide; "Figures 10 and 7".
