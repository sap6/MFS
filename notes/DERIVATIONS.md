# Derivations behind the source layouts

Working notes for the author (not manuscript text). They write out the three
arguments that `NOTES.md` listed as unproved, with the numerical check for
each. Notation: time dependence exp(-i omega t), wavenumber k, perfectly
conducting body Omega, incident field E^i (entire), scattered field E^s,
total field E = E^i + E^s outside Omega, n x E = 0 on the surface.

## A. A convex polyhedron: the scattered field is a layer on the skeleton

**Claim.** Let Omega be a convex polyhedron with faces F. Split Omega into the
cells C_F = {x in Omega : F is the nearest face plane}; their common walls form
the interior medial surface M (for the cube: the twelve triangles that join
each edge to the centre, each lying in the bisector plane of its edge). Then
the scattered field outside Omega is exactly the field of an electric and a
magnetic surface current on M.

**Proof.**

1. *Reflection.* For a face F with unit normal n let x* be the mirror image
   of x in the plane of F and R = I - 2 n n^T. If (E, H) solves Maxwell's
   equations in a region U, then E'(x) = -R E(x*), H'(x) = R H(x*) solves them
   in the mirror image of U (E is a polar, H an axial vector). On the plane,
   E' has the tangential part -E_t and the normal part E_n; H' has H_t and
   -H_n.
2. *Continuation across a face.* On F the total field has E_t = 0 and
   therefore H_n = 0, so E' = E and H' = H on F: the field defined as (E, H)
   outside and (E', H') inside has continuous tangential traces across the
   open face and solves Maxwell's equations on both sides, hence across it.
   This is the image principle read backwards. Call the continuation E_F.
3. *It reaches the whole body.* E_F(x) = -R E(x*) is defined wherever x* lies
   outside Omega. For convex Omega the mirror image of Omega in the plane of F
   lies outside Omega, so E_F is a Maxwell solution in all of Omega; it is
   smooth up to the open face F and singular only at the edges of F.
4. *Patch.* Define V = E outside Omega and V = E_F in the cell C_F. V solves
   Maxwell's equations in the complement of the closure of M and is analytic
   across every open face. U = V - E^i equals E^s outside Omega and radiates.
5. *Representation.* Apply the Stratton-Chu formula to U on the two sides of
   M. The contributions of small tubes around the edges and small spheres
   around the vertices vanish, because the field grows only like
   (distance)^(-1/3) at an edge and (distance)^(nu - 1), nu = 0.454, at a
   vertex (section D), and the continuations inherit the same growth by
   reflection. What remains is

       E^s = field of the surface currents  J = n_M x [H],  K = -n_M x [E]  on M,

   where [.] is the jump between the continuations from the two neighbouring
   faces.  QED

**Consequences used in the code.**

- Sources on the skeleton are *complete* for the polyhedron, and nothing
  closer to the surface than the skeleton is needed. A shrunk copy of the
  surface cannot work: the edges are on the skeleton.
- Both electric and magnetic sources are needed (J and K are both non-zero
  in general): the numerical finding that electric-only edge sources stall.
- The densities are singular only towards the edges and vertices, which is
  where the sources are clustered.

**Status.** Steps 1-2 are the classical image principle; steps 3-5 are
elementary but I have not found this statement for polyhedra in the
literature (the continuation literature, Millar 1971-1986, Kyurkchan,
Sternin and Shatalov 1996, treats analytic boundaries and locates the
singularities of the continuation, which is the same idea for smooth bodies).
Numerical evidence: the cube solver, whose sources lie on M, converges to the
graded boundary-element solution (`NOTES.md`). Non-convex polyhedra need more
care (the mirror image of the body can re-enter the body) and are not covered.

## B. Exponents of the edge tail

Take one edge, with local coordinates s (along the edge), rho (distance to
the edge) and the wall of M that ends on it: the interior bisector half-plane.
Let m be the mirror in the bisector plane; it maps the exterior wedge (opening
3 pi / 2) to itself and swaps the two faces.

1. *Edge expansion (Meixner).* With phi measured from one face through the
   exterior, 0 < phi < 3 pi / 2, the field near the edge is a sum of

       TM-type:  E_s = a_n(s) rho^(2n/3) sin(2n phi / 3),  n = 1, 2, ...
       TE-type:  H_s = b_n(s) rho^(2n/3) cos(2n phi / 3),  n = 0, 1, 2, ...

   each multiplied by a power series in rho^2; the transverse components
   follow from E_s and H_s and are one power of rho lower (rho^(2n/3 - 1)).
2. *Continuation.* The reflection of step A.1 continues each term to the
   interior by the same formula, phi < 0 through one face and phi > 3 pi / 2
   through the other. With psi = phi - 3 pi / 4 (angle from the exterior
   bisector) the wall is psi = +pi seen from one side and psi = -pi from the
   other. The angular factors are cos or sin of 2 n psi / 3, so the jump
   across the wall is proportional to sin(2 n pi / 3) for the odd part and
   zero for the even part. In particular **terms with n a multiple of 3
   (integer powers) do not jump at all**: they are polynomial fields and need
   no sources.
3. *Parity.* A term is even under m (E(m x) = m E(x)) when E_s is even in
   psi or H_s is odd in psi: TM with n odd, TE with n odd. For an even field
   the tangential E is continuous across the wall and the tangential H jumps:
   **only electric current**. A term is odd under m for TM with n even and TE
   with n even: tangential H continuous, tangential E jumps: **only magnetic
   current**. Surface currents are tangential, so there is **no component
   normal to the bisector plane**.
4. *Densities.* Reading off the jumps (rho now the depth on the wall):

   | current on the wall | from | leading term | density | tail exponent e |
   |---|---|---|---|---|
   | electric, along the edge, J_s = [H_rho] | TM n = 1 (and TE n = 1 through d/ds) | rho^(2/3 - 1) | rho^(-1/3) | 2/3 |
   | electric, along the bisector, J_rho = [H_s] | TE n = 1 | rho^(2/3) | rho^(2/3) | 5/3 |
   | magnetic, along the edge, K_s = [E_rho] | TE n = 2, TM n = 2 | rho^(4/3 - 1) | rho^(1/3) | 4/3 |
   | magnetic, along the bisector, K_rho = [E_s] | TM n = 2 | rho^(4/3) | rho^(4/3) | 7/3 |

   Next terms: 8/3, 11/3, 8/3, 11/3 (from the rho^2 corrections and n = 4, 5).
5. *From density to weights.* The tail replaces the sheet below the deepest
   free depth d_J by lines at depths d_i = d_J 2^(-i). A line stands for the
   strip of the sheet around it, of width proportional to d_i, so a density
   rho^(e - 1) gives line weights (d_i / d_J)^e. That is the weight used in
   `TailLine`, and the exponents in the last column are the ones to use.

**Check** (`check_tail_exponents.py`, balanced family without vertex tail,
error against the J = 9 field): the derived set with one exponent per
component, (2/3 | 5/3 | 4/3 | 7/3), gives 0.314% at J = 5 and 0.049% at J = 7;
the set used before, which had a second, integer exponent per component
(1, 2, 1, 2), gives 0.327% and 0.049% with 8-10% more unknowns; adding the
derived second terms (8/3, 11/3) gives 0.316% at J = 5. The integer exponents
were not physical (item 2) and are dropped. Exponents one too low for the
transverse components gave garbage in an earlier session, in line with the
table.

**Check repeated on 18 Nov 2026 with the layout of the paper** (two depths
per octave, vertex tail, error against the J = 12 reference;
`check_tail_exponents.py` now uses `cube_mfs.py`): derived set 8.31e-4 at
J = 6 and 1.40e-4 at J = 8; with the second terms 7.99e-4 and 1.35e-4 (4%
and 3% lower); with the exponents of the two components along the bisector
lowered by one 3.45e-3 and 1.33e-3 (4 and 10 times larger); without the
vertex tail 1.33e-3 and 2.27e-4 (1.6 times larger).

## C. Exponents of the almond's tip

The tip is the vertex of an elliptic cone (half-angles 34.09 and 12.71
degrees). Near the vertex of a perfectly conducting cone the field is a sum
of two families (Van Bladel's book, ch. 5; De Smedt and Van Bladel 1986):

    electric type:  E = grad( r^nu  Y ),   E ~ r^(nu - 1)
    magnetic type:  H = grad( r^tau Z ),   H ~ r^(tau - 1),  E ~ r^tau

with Y (Z) a Dirichlet (Neumann) eigenfunction of the Laplace-Beltrami
operator on the part of the unit sphere outside the cone, eigenvalue
nu (nu + 1) (tau (tau + 1)). `cone_exponents.py` computes them (method of
particular solutions on the sphere, poles on the focal segment; exact to
8 digits for the circular cone, unchanged to 5 digits from 12 to 36 poles):

| parity of the potential (y, z) | nu (Dirichlet) | tau (Neumann) |
|---|---|---|
| even, even | **0.29959**, 1.46427, 2.02232 | 1.05604, 1.98580 |
| odd, even | 1.10372, 2.24730 | 0.96803, 1.96545 |
| even, odd | 1.04120, 2.10576 | 0.87905, 1.75454 |
| odd, odd | 2.01295, 3.04849 | 1.97538, 2.90094 |

Only one exponent is below 1: nu = 0.29959, so E ~ r^(-0.70) at the tip, and
only for fields that are even under both mirrors of the almond (the class
that drives current along the axis into the tip). Circular cones with the two
half-angles would give 0.377 and 0.222.

*Which potential belongs to which class.* A class (sy, sz) is defined by
E(R x) = s R E(x) for each mirror R. E = grad(r^nu Y) has this parity when Y
has the parities (sy, sz). H is an axial vector, H(R x) = -s R H(x), so the
magnetic potential Z has the parities (-sy, -sz).

*From field to dipole weights.* The tail sources lie on rows that run into
the tip with spacing proportional to r. Dipoles of strength r^e on such a row
produce a field r^(e - 3) at distance r from the tip. Hence electric dipoles
with e = nu + 2 for the electric family and magnetic dipoles with e = tau + 2
for the magnetic one (`tip_exponents` in `pec_almond_skeleton.py`).

**Check** (1.19 GHz, 19 look angles): the single computed exponent
e = 2.29959 for the class (+, +) gives the same boundary residual (0.24%) as
the fitted set {2.2, 2.4, 2.7, 3.0} used before, and the HH cross-section is
closer to the fine-layout result (0.0001 dB rms instead of 0.0003). More
exponents (second electric, first magnetic) change nothing at this accuracy.
The fitted set is no longer needed.

## D. Exponents at a vertex of the cube

Same two families as in C, on the sphere minus the spherical triangle cut out
by the cube (`python3 cone_exponents.py cube`):

    nu  (Dirichlet): 0.45417, 1.23087, 1.78427, 2.11761, 2.51450
    tau (Neumann):   0.840,  1.20635, 1.8061,  2 (regular), 2.44666

So E ~ r^(-0.546) (electric type) and H ~ r^(-0.160) (magnetic type) at a
vertex. 0.45417 is the known value for the exterior of a cube corner. They
bound the growth needed in step A.5.

**Vertex tail.** The free point dipoles on the body diagonal stop at the
deepest depth d_J, like the edge lines. Continuing them on the diagonal at
distances d_J 2^(-i) with weights (d_i / d_J)^e, e = nu + 2 = 2.454 for the
electric and e = tau + 2 = 2.840 for the magnetic dipoles (same counting as
in C), adds 12 unknowns in the symmetric solver (`VertexTail`).

**Check.** Error against the J = 9 field without / with the vertex tail:
J = 5: 0.31% / 0.20%; J = 6: 0.133% / 0.086%; J = 7: 0.049% / 0.029%; against
graded Bempp (10368 unknowns) at J = 7: 0.059% / 0.042%. A second exponent
per kind (nu_2, tau_2) changes nothing. The tail on the diagonal carries only
the radial law; the angular structure of the vertex singularity is left to
the edge lines, which is presumably why the gain is a constant factor (about
half a refinement level) rather than a change of rate.

## E. Torus: where the continued field is singular

Torus (rho - R)^2 + z^2 = a^2 about the z axis, a < R. Write
c = sqrt(R^2 - a^2): the ring rho = c, z = 0 is the focal ring of the toroidal
coordinates in which this torus is a coordinate surface (it lies inside the
tube, between the inner wall and the centre circle).

**Claim (argued, then checked numerically; not a theorem here).** The
scattered field of a plane wave continues analytically into the tube except
at the centre circle rho = R and at the focal ring rho = c (and possibly
between them); seen from the centre circle the singular set reaches out to
the distance

    delta = R - c = R - sqrt(R^2 - a^2).

A source torus concentric with the tube must therefore have minor radius
a_s > delta. As a fraction of a: 0.127 (a/R = 0.25), 0.268 (0.5), 0.5 (0.8),
tending to 1 for a -> R. A fat torus leaves little room.

**Why.** Three ways to see it, none a full proof for Maxwell:

1. *Reflection.* In the meridian half-plane write zeta = rho + i z. The
   Schwarz function of the meridian circle is S(zeta) = R + a^2 / (zeta - R).
   In two dimensions the continuation of the field is singular only where S
   is, at zeta = R. For a body of revolution the reduced equation has a
   singular coefficient on the axis rho = 0, and the reflection formula picks
   up a second singularity where the reflected point falls on the axis:
   zeta + S(zeta) = 0, i.e. zeta^2 = R^2 - a^2, zeta = c. (General theory:
   Millar, Arch. Ration. Mech. Anal. 81 (1983), the axially symmetric
   Helmholtz equation; not read, publisher access needed.)
2. *Images.* In a plane through the axis the tube is cut in two discs. The
   singularity at the centre of one disc is, for the other disc, an outside
   source; its image lies at distance a^2 / (2R) from that disc's centre,
   towards the axis; the image of the image at a^2 / (2R - a^2/(2R)), and so
   on. The sequence d -> a^2 / (2R - d) increases to its fixed point
   d = R - c. So the singularities between centre circle and focal ring
   accumulate at the focal ring.
3. *Static limit.* For Laplace's equation the field outside the torus is a
   series of toroidal harmonics P_{n-1/2}(cosh tau) exp(i n sigma), each
   singular only on the focal ring. With boundary data from an entire
   function the coefficients decay like exp(-2 n tau_0) (tau_0 the surface),
   so the series converges outside the torus tau = 2 tau_0, whose meridian
   circle has the diameter [R - a^2/R, R] and contains both c and R. The
   singular set of the Helmholtz and Laplace continuations has the same
   convex hull (cited by Hochman, Leviatan and White 2013 from the MAS
   literature).

**Check** (`python3 torus_singularities.py`, k = 2 pi unless stated). Sources
on a circle of radius b about (rho0, 0) in the meridian plane; if the
continuation is analytic outside the distance delta from (rho0, 0), the n-th
Fourier coefficient of the dipole density around that circle decays like
(delta / b)^n. Fitted delta for a/R = 0.8 (c = 0.6), source circle
(1.0, 0.60): wave along the axis 0.395-0.401; in the plane of the ring
0.39-0.40; oblique 0.391-0.400; k = 4 pi 0.398 (one clean window). Predicted
0.400; the first image alone would give 0.32, the point R - a^2/R 0.64.
With the centre moved outwards to rho0 = 1.1 the fitted distance becomes
0.49 (predicted |1.1 - c| = 0.50): the singularity is on the side of the
axis, not symmetric about the centre circle. For a/R = 0.5 and 0.25 the
decay is too fast to fit cleanly (0.15-0.17 against 0.134; 0.03-0.04 against
0.032): consistent, not decisive.

Second check (`fig_torus_sources.py`): with a fixed number of sources the
coefficient norm stays of order 1 down to the predicted radius and then
grows by orders of magnitude: a/R = 0.8: 1.9, 3.6, 41, 7e4 for a_s/a = 0.55,
0.50, 0.45, 0.40 (predicted 0.5); a/R = 0.5: 13, 49, 460, 2e4 for 0.30, 0.25,
0.20, 0.15 (0.268); a/R = 0.25: 24, 770, 2e4 for 0.15, 0.10, 0.07 (0.127).
This is the behaviour Barnett and Betcke prove for the disc (their Theorems
6 and 7).

**Status and literature.** The enclosure rule itself is old (MAS "main
singularities", Kyurkchan, Sternin and Shatalov; Barnett and Betcke). For a
torus with an elliptic meridian section Karkashadze and Zaridze place the
auxiliary contour around the foci of the ellipse, which for a circular
section is the centre circle only. In even dimensions the Schwarz potential
of the torus is singular only on the centre circle (Aberra 2002, as quoted
by Lundberg, arXiv:1009.5159); the extra ring in three dimensions comes from
the logarithmic branch of the Riemann function. I have not found the focal
ring stated for the 3D torus; Millar 1983 and Khavinson and Shapiro's
Schwarz-potential papers must be read before anything is claimed.
