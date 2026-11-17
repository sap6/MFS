# Literature check (task T6), 2026-11-16

How each item was checked is stated, because the levels differ:

- **[full]** full text read (arXiv PDF) for the points that matter here;
- **[abs]** abstract only (publisher sites cannot be reached from the cloud
  session; abstracts via Crossref / OpenAlex);
- **[meta]** authors, title, journal, volume, pages, year and DOI confirmed
  against the Crossref record, content not read in this session.

Hosts that work from the cloud session: arxiv.org, api.crossref.org,
api.openalex.org (rate-limited), github.com (public clones). Blocked:
sciencedirect.com, dl.acm.org, mdpi.com, ieeexplore (not tried, same class).

## 1. Is a 3D edge-singular ("lightning-type") MFS for Maxwell new?

Nothing found that does it. Evidence:

- arXiv searches (abstract fields) for lightning + 3D/polyhedra, fundamental
  solutions + edges + 3D, auxiliary sources + edges, MFS + cube/polyhedron +
  Helmholtz/scattering, MFS + Maxwell: the only hits are the papers below,
  none with sources clustered at the edges of a 3D body.
- Kong and Shen, arXiv:2610.26540 (October 2026) **[full, 3D section and
  conclusion]**: lightning method for the modified Helmholtz equation; its 3D
  part treats *isolated point* singularities only and says of the rest:
  "Although this approach can be in principle extended to deal with
  three-dimensional polyhedral edge and vertex singularities, how to
  efficiently solve the large scale least-squares problem required to resolve
  these singularities is particularly challenging."
- Brubeck and Trefethen, lightning Stokes solver **[full, remark]**:
  "Extension of these methods to three dimensions is undoubtedly possible, but
  not straightforward".
- Ginn and Trefethen, arXiv:2310.01665 **[full, method sections]**: 2D
  exterior Helmholtz; Hankel-function poles clustered exponentially along the
  interior bisectors of the corners plus a smooth part. This is the 2D
  version of our cube layout and must be cited as such.
- Weed, Wu, Huang, Cho, arXiv:2605.15527 **[full, sections 1-2]** (the paper
  the notes could not identify): MFS for Maxwell *transmission* problems in
  bi-periodic layered media, electric dipole sources offset from *smooth*
  interfaces, periodising proxy sources, errors near 1e-14. No edges or
  corners, no PEC body. It shows that MFS for 3D Maxwell on smooth surfaces is
  current work; it does not overlap with the cube or almond material.
- Ye, Hu, Xu, Yan, EABE 155 (2023) 401-406 **[meta + search snippet]**:
  regularised MFS for 3D PEC targets; sources *on* the surface with
  subtraction/add-back of the singularity, which avoids choosing an auxiliary
  surface. A different route (no exponential convergence, no edge treatment
  stated). To read before submission.
- Tabatadze, Iscan, Karacuha, J. Electromagn. Waves Appl. 39 (2025) 952-971
  **[abs]**: MAS for PEC square plates; edges are handled by wrapping the
  plates in an auxiliary dielectric interface, not by singular or clustered
  sources.
- Chen, Schaefer, Desbrun, "Lightning-fast method of fundamental solutions",
  ACM TOG 43(4) 2024 **[abs]**: a fast preconditioner for MFS/BEM systems;
  "lightning" refers to speed. Not a singularity method. Do not cite as a
  precedent.
- Levin, arXiv:2603.09926 **[full, skim]**: interior harmonic Dirichlet
  problems on the cube with singular *data*, via an image-series Green's
  function. Unrelated to edge singularities of scattering.
- Cai, Bagge, Martinsson, arXiv:2603.19113 **[full, sections 1-3]**: MFS
  scattering matrices for many bodies; in 2D, corners are handled by dyadic
  refinement with the source distance shrinking towards the corner, and the
  kinship with lightning solvers is noted; 3D examples are smooth.

Closest precedents for building the edge behaviour into a source model (all
2D):

- Eisler and Leviatan 1989 **[meta]**, fictitious filaments placed near
  edges; Shifman, Friedmann, Leviatan 1997 **[abs]**: "Previous approaches
  have required that some of the fictitious currents be located in the
  proximity of the edges... the accuracy... dependent on a judicious choice of
  source location"; they expand the near-edge field in infinite-wedge
  solutions and couple it to filaments elsewhere (a hybrid, not an MFS).
- Iatropoulos, Anastasiadou, Anastassiu 2020 **[abs]**: MAS sources moved
  towards the tip of curved wedges improve accuracy (2D, TM).
- Antunes and Valtchev 2010 **[meta]**; Marin 2010 **[meta]**: MFS enriched
  with corner functions, 2D Helmholtz. (The "J. Comput. Appl. Math. 470"
  item in the earlier notes could not be confirmed and should not be cited.)
- Barnett and Betcke 2010 **[meta]**: fundamental solutions plus
  Fourier-Bessel corner functions for polygons.

What is different here, stated carefully: (i) three dimensions and the
vector problem; (ii) the singular behaviour is carried by fixed-weight sums
of fundamental solutions ("tails") with exponents derived from the edge and
cone expansions (`DERIVATIONS.md`), so the method stays a pure MFS;
(iii) the same device at a conical tip. **Claim that is safe:** "to our
knowledge the first MFS for 3D Maxwell scattering that converges on bodies
with edges and tips". Not safe: "first lightning method in 3D" (Kong and
Shen have a 3D point-singularity version).

## 2. Where the sources must be: continuation of the field

- Barnett and Betcke 2008 **[meta; arXiv:0708.3533 available]**: MFS for
  Helmholtz converges exponentially and stays stable when the source curve
  encloses the singularities of the continued solution. The sphere study in
  the paper is the Maxwell analogue; cite as the model.
- Millar 1971, 1976, 1986 **[meta]**; Kyurkchan, Sternin, Shatalov 1996
  **[meta]**: location of the singularities of the continued scattered field
  for analytic boundaries.
- Leviatan 1990 **[meta]**; Zaridze et al. 2002 **[meta]**; Kaklamani and
  Anastassiu 2002 **[meta]**; Fikioris 2006 **[meta]**; Fikioris and Tsitsas
  2015 **[meta]**; Tsitsas, Zouros, Fikioris, Leviatan 2018 **[meta]**: the
  MAS literature on the same condition, and on divergent, oscillating
  auxiliary currents when it is violated (our failed shrunk-surface layouts).
- Doicu, Eremin, Wriedt 2000 (book) **[meta, chapter records]**: discrete
  sources on the axis or in the symmetry plane of the particle; the nearest
  relative of the almond layout. Hafner / Wriedt 1999 (GMT) **[meta]**.
- The statement and proof that the field of a convex PEC polyhedron is a
  layer on the interior medial surface (`DERIVATIONS.md`, A) were not found
  anywhere; it rests on the classical image principle. Present it as an
  observation with a short proof.

## 3. Singularity exponents

- Meixner 1972 **[meta]**: edge condition and expansion. Costabel and Dauge
  2000 **[meta]**: Maxwell singularities in polyhedra (edges, corners; the
  rigorous background for D).
- Van Bladel 1983; De Smedt and Van Bladel 1986; De Smedt 1988; Vafiadis and
  Sahalos 1984 **[meta]**; Van Bladel's book (1991, IEEE reissue 1996)
  **[meta]**: electric (Dirichlet) and magnetic (Neumann) tip exponents of
  cones, including elliptic cones and pyramidal corners. Boersma and Jansen
  1990 (Eindhoven report 90-WSK-01, no DOI) **[search snippet only]**.
  Our numbers (`cone_exponents.py`): almond tip nu = 0.29959; cube vertex
  nu = 0.45417, tau = 0.840. The value 0.4542 for the exterior of a cube
  corner is standard (Fichera 1975 **[meta]**); the tables of De Smedt and
  Van Bladel should be compared with ours once the paper is at hand.
- Betcke and Trefethen 2005 **[meta]**: the method of particular solutions
  used for the exponents.

## 4. Benchmarks and software

- Woo, Wang, Schuh, Sanders 1993 **[meta]**: almond definition.
- Kelley, Chamulak, Courtney, Yilmaz 2020 **[meta; data in the repo]**:
  measured almond RCS and the Austin benchmark reference data.
- Betcke and Scroggs 2021 (Bempp-cl) **[meta]**; Scroggs et al. 2017
  **[meta]**; Rao, Wilton, Glisson 1982 **[meta]**; Buffa and Hiptmair 2003
  **[meta]**; Bespalov, Heuer, Hiptmair 2010 **[meta]** (EFIE convergence on
  polyhedral surfaces; background for the slow rate of quasi-uniform meshes
  and the use of graded ones).
- MFS for Maxwell before this work: Young and Ruan 2005 (CMES 7, 223-232;
  DOI 10.3970/cmes.2005.007.223 is *not* in Crossref, check by hand);
  Liu's Dartmouth thesis 2016 (acoustic and electromagnetic periodic
  scattering, cited by Cai et al.; not seen); Liu and Barnett 2016 **[meta]**
  (acoustic, axisymmetric).

## 5. Verified reference data (Crossref record matched)

    Barnett, Betcke. J. Comput. Phys. 227 (2008) 7003-7026. 10.1016/j.jcp.2008.04.008
    Barnett, Betcke. SIAM J. Sci. Comput. 32 (2010) 1417-1441. 10.1137/090768667
    Gopal, Trefethen. SIAM J. Numer. Anal. 57 (2019) 2074-2094. 10.1137/19M125947X
    Gopal, Trefethen. PNAS 116 (2019) 10223-10225. 10.1073/pnas.1904139116
    Brubeck, Trefethen. SIAM J. Sci. Comput. 44 (2022) A1205-A1226. 10.1137/21M1408579
    Trefethen, Nakatsukasa, Weideman. Numer. Math. 147 (2021) 227-254. 10.1007/s00211-020-01168-2
    Herremans, Huybrechs, Trefethen. SIAM J. Numer. Anal. 61 (2023) 2580-2600. 10.1137/23M1551821
    Ginn, Trefethen. arXiv:2310.01665 (2023)
    Kong, Shen. arXiv:2610.26540 (2026)
    Weed, Wu, Huang, Cho. arXiv:2605.15527 (2026)
    Cai, Bagge, Martinsson. arXiv:2603.19113 (2026)
    Kupradze, Aleksidze. USSR Comput. Math. Math. Phys. 4(4) (1964) 82-126. 10.1016/0041-5553(64)90006-0
    Fairweather, Karageorghis. Adv. Comput. Math. 9 (1998) 69-95. 10.1023/A:1018981221740
    Fairweather, Karageorghis, Martin. Eng. Anal. Bound. Elem. 27 (2003) 759-769. 10.1016/S0955-7997(03)00017-1
    Liu, Barnett. J. Comput. Phys. 324 (2016) 226-245. 10.1016/j.jcp.2016.08.011
    Ye, Hu, Xu, Yan. Eng. Anal. Bound. Elem. 155 (2023) 401-406. 10.1016/j.enganabound.2023.05.050
    Kaklamani, Anastassiu. IEEE Antennas Propag. Mag. 44(3) (2002) 48-64. 10.1109/MAP.2002.1028734
    Fikioris. IEEE Trans. Antennas Propag. 54 (2006) 2022-2033. 10.1109/TAP.2006.877171
    Fikioris, Tsitsas. Comput. Math. Appl. 69 (2015) 636-649. 10.1016/j.camwa.2015.01.014
    Tsitsas, Zouros, Fikioris, Leviatan. IEEE Trans. Antennas Propag. 66 (2018) 5443-5452. 10.1109/TAP.2018.2855963
    Leviatan. IEEE Trans. Antennas Propag. 38 (1990) 1259-1263. 10.1109/8.56964
    Zaridze, Bit-Babik, Tavzarashvili, Economou, Uzunoglu. IEEE Trans. Antennas Propag. 50 (2002) 50-58. 10.1109/8.992561
    Anastassiu, Kaklamani, Economou, Breinbjerg. IEEE Trans. Antennas Propag. 50 (2002) 59-66. 10.1109/8.992562
    Eisler, Leviatan. IEE Proc. H 136(6) (1989) 431. 10.1049/ip-h-2.1989.0079
    Shifman, Friedmann, Leviatan. IEE Proc. Microw. Antennas Propag. 144(4) (1997) 235. 10.1049/ip-map:19971186
    Iatropoulos, Anastasiadou, Anastassiu. Appl. Sci. 10 (2020) 2309. 10.3390/app10072309
    Tabatadze, Iscan, Karacuha. J. Electromagn. Waves Appl. 39 (2025) 952-971. 10.1080/09205071.2025.2490919
    Antunes, Valtchev. J. Comput. Appl. Math. 234 (2010) 2646-2662. 10.1016/j.cam.2010.01.031
    Marin. Appl. Math. Model. 34 (2010) 1615-1633. 10.1016/j.apm.2009.09.009
    Doicu, Eremin, Wriedt. Acoustic and Electromagnetic Scattering Analysis Using Discrete Sources. Academic Press 2000. (chapters: 10.1016/B978-012219740-6/...)
    Wriedt (ed.). Generalized Multipole Techniques for Electromagnetic and Light Scattering. Elsevier 1999. 10.1016/B978-0-444-50282-7.X5012-X
    Bogdanov, Karkashadze, Zaridze. In the same volume, 143-172. 10.1016/B978-044450282-7/50019-1
    Kyurkchan, Sternin, Shatalov. Physics-Uspekhi 39 (1996) 1221-1242. 10.1070/PU1996v039n12ABEH000184
    Millar. Math. Proc. Camb. Phil. Soc. 69 (1971) 175-188. 10.1017/S0305004100046545
    Millar. SIAM J. Math. Anal. 7 (1976) 131-156. 10.1137/0507012
    Millar. IMA J. Appl. Math. 37 (1986) 155-171. 10.1093/imamat/37.2.155
    Meixner. IEEE Trans. Antennas Propag. 20 (1972) 442-446. 10.1109/TAP.1972.1140243
    Costabel, Dauge. Arch. Ration. Mech. Anal. 151 (2000) 221-276. 10.1007/s002050050197
    Van Bladel. Proc. IEEE 71 (1983) 901-902. 10.1109/PROC.1983.12686
    Vafiadis, Sahalos. Proc. IEEE 72 (1984) 1089-1091. 10.1109/PROC.1984.12976
    De Smedt, Van Bladel. IEEE Trans. Antennas Propag. 34 (1986) 865-870. 10.1109/TAP.1986.1143916
    De Smedt. IEEE Trans. Antennas Propag. 36 (1988) 152-155. 10.1109/8.1089
    Van Bladel. Singular Electromagnetic Fields and Sources. IEEE Press reissue 1996. 10.1109/9780470546420
    Fichera. Russian Math. Surveys 30(3) (1975) 107-127. 10.1070/RM1975v030n03ABEH001413
    Betcke, Trefethen. SIAM Review 47 (2005) 469-491. 10.1137/S0036144503437336
    Woo, Wang, Schuh, Sanders. IEEE Antennas Propag. Mag. 35(1) (1993) 84-89. 10.1109/74.210840
    Kelley, Chamulak, Courtney, Yilmaz. IEEE Antennas Propag. Mag. 62(1) (2020) 96-106. 10.1109/MAP.2019.2955702
    Betcke, Scroggs. J. Open Source Softw. 6(59) (2021) 2879. 10.21105/joss.02879
    Scroggs, Betcke, Burman, Smigaj, van 't Wout. Comput. Math. Appl. 74 (2017) 2897-2914. 10.1016/j.camwa.2017.07.049
    Rao, Wilton, Glisson. IEEE Trans. Antennas Propag. 30 (1982) 409-418. 10.1109/TAP.1982.1142818
    Buffa, Hiptmair. In: Lect. Notes Comput. Sci. Eng. 31 (2003) 83-124. 10.1007/978-3-642-55483-4_3
    Bespalov, Heuer, Hiptmair. SIAM J. Numer. Anal. 48 (2010) 1518-1529. 10.1137/090766620
    Mie. Ann. Phys. 330 (1908) 377-445. 10.1002/andp.19083300302
    Bohren, Huffman. Absorption and Scattering of Light by Small Particles. Wiley 1998 reissue. 10.1002/9783527618156
    Chen, Schaefer, Desbrun. ACM Trans. Graph. 43(4) (2024). 10.1145/3658199   (only to avoid confusion over the name)

Two DOIs recalled from memory were wrong and are corrected above (Eisler and
Leviatan; De Smedt and Van Bladel): check every DOI against Crossref when the
bibliography is typed.

## 6. Still to do by hand (needs publisher access)

Read: Ye et al. 2023; Tsitsas et al. 2018; Shifman et al. 1997; De Smedt and
Van Bladel 1986 (compare exponent tables); Kelley et al. 2020 (measurement
uncertainty, to quote next to the 0.4-1.5 dB differences); Liu's thesis
(MFS for Maxwell on axisymmetric bodies); Young and Ruan 2005.

## 7. Second pass, 2026-11-16 (third session): full texts

The author asked for a thorough reading rather than abstracts. What was read
in full (the passages that matter here, and the surrounding sections), and
what it changes. Quotations are copied from the PDFs.

### 7.1 Is the 3D edge case open? Yes, and it is posed as a question

- Gopal and Trefethen, SIAM J. Numer. Anal. 2019 **[full]**, closing
  section: "The extension to 3D is a bigger issue to be considered. Here,
  as in 2D, there are MFS papers going back many years. The question is
  whether a new more focussed connection with approximation theory near
  singularities can bring a new efficiency, as we have shown here for 2D
  problems, and we are open-minded on this question." This is the sentence
  the paper answers; quote it.
- Brubeck and Trefethen 2022 **[full]**: "Extension of these methods to three
  dimensions is undoubtedly possible, but not straightforward, and we have no
  view as to whether something practical can be developed here or not."
- Kong and Shen, arXiv:2610.26540 **[full: abstract, contributions, section 6,
  conclusion]**: the 3D part places modified-Bessel poles around *point*
  singularities only; conclusion as quoted in section 1 above (verified
  verbatim).
- Hochman, Leviatan, White, J. Comput. Phys. 238 (2013) 337-358
  (arXiv:1112.1643) **[full: abstract, conclusion]**: dipole sources (poles
  of a rational function) placed by vector fitting for 2D Laplace problems
  with sharp corners, errors near machine precision; a forerunner of the
  lightning solvers, in the target journal, and from the MAS community.
  Conclusion: the method "may be difficult to generalize to 3D Laplace
  problems". Must be cited. It also records that the convex hull of the
  singularities of the Helmholtz and Laplace continuations is the same.
- Cai, Bagge, Martinsson, arXiv:2603.19113 **[full: sections 1-3, 5.2]**:
  corners in 2D only, by dyadic refinement with the MFS distance shrinking
  towards the corner. Its 3D examples are smooth ellipsoids and **tori**
  (acoustic, sources at a fixed normal distance d = 0.18, errors to 4e-7),
  so a torus as an MFS test body is not new.
- Ginn and Trefethen, arXiv:2310.01665 **[full: front matter, sections 1-3]**:
  this is Henry Ginn's Oxford MMath dissertation (supervisor Trefethen; the
  arXiv record lists both as authors). Exterior 2D Helmholtz, Hankel poles
  along the interior bisectors. Cite as an arXiv preprint.
- Weed, Wu, Huang, Cho, arXiv:2605.15527 **[full: sections 1-2, 6]**:
  confirmed as described in section 1: smooth bi-periodic interfaces, dipoles
  offset along the normals, no edges.
- Kouroublakis, Tsitsas, Leviatan, arXiv:2607.28558 (July 2026) **[full:
  abstract, section III]**: MAS for metasurfaces in **2D**; at the edges of a
  finite sheet "sources are positioned closer to the edges and gradually
  placed farther away toward the interior" and extended beyond the edges.
  The current state of edge handling in MAS: ad hoc clustering, two
  dimensions, no singular functions.

### 7.2 What the MAS community does with edges in 3D: it rounds them

- Bogdanov, Karkashadze, Zaridze, "The method of auxiliary sources in
  electromagnetic scattering problems", in Wriedt (ed.) 1999 **[full, author
  copy at lae.tsu.ge]**. Section 3 defines the "main singularities" (the
  singularities of the continued field, interpreted as images of the primary
  sources) and section 4 gives the rules for the auxiliary surface: (2) it
  "should embrace the domain of location of the main singularities"; (3) the
  distance to the surface must be below the minimal radius of curvature,
  which "also indicates a way of investigating scattering problems upon
  bodies with non-smooth surfaces (edges), offering at edge points the
  finite rounding radius"; (4) the interior of the auxiliary surface must
  not resonate at the working frequency.
- Karkashadze and Zaridze, "The method of auxiliary sources in applied
  electrodynamics" (Tbilisi; a review of the mid-1990s, author copy at
  lae.tsu.ge) **[full]**: recommendation (d), "In the case of edges at the
  main surface ... It is more advantageous to smooth the edges under
  certain radius of curvature and to trace out its influence on the final
  results." Also: "for better efficiency of the method, we must locate the
  auxiliary sources in the scattered field's singularity region, i.e. in
  convexity centers of the scattered field or in basic contour angles (in
  the case of a surface with edges)" (stated for 2D contours, not carried
  out in 3D). **And a torus:** ring sources (their eq. 32) and the
  "excitation of the elliptic torus by magnetic dipole", where "the
  auxiliary contour l surrounds the foci of the given ellipse". So MAS on a
  torus exists (axisymmetric excitation); the rule used is the 2D one for
  the meridian section.
- Zaridze et al., IEEE TAP 50 (2002) 50-58 and the 2005 AP-S paper **[full]**:
  scattered-field singularities sit near caustic surfaces; "It is good to
  avoid sharp corners in MAS."
- Consequence for the paper: the sentence "the MAS literature treats edges
  by rounding them or by wrapping the body in an auxiliary surface
  (Tabatadze et al. 2025)" is accurate and can be supported by these
  quotations. Keeping the edge sharp and carrying its singularity by
  sources with derived power laws was not found anywhere.

### 7.3 The enclosure rule and its proofs

- Barnett and Betcke 2008 **[full: sections 2-3, arXiv:0708.3533]**. Unit
  disc, charge radius R, singularity radius rho: Theorem 3 (error
  C rho^(-N/2) if rho < R^2, C R^(-N) if rho > R^2), Theorem 6 (coefficients
  bounded if R < rho), Theorem 7 (if R > rho the coefficient norm grows at
  least like sqrt(N) (R/rho)^(N/2)). General analytic domains: Conjecture 9
  and 12, with Millar's result that the singularities of the Schwarz
  function are singularities of the continuation. The torus radius study
  (item 34 of the notes) is the 3D analogue of their Figure 3.
- Liu and Barnett 2016 **[full: sections 1, 3.2, arXiv:1506.05083]**: MFS on
  axisymmetric bodies mode by mode ("ring sources"), the same splitting as
  `pec_torus.py`; they state the enclosure condition as a conjecture and
  cite, for 3D Maxwell, **Kangro, AIP Conf. Proc. 1479 (2012) 2328-2331**
  ("Solution of three-dimensional electromagnetic scattering problems by
  interior source methods"): a convergence proof for analytic boundaries
  when the sources are the nodes of an exponentially convergent quadrature
  on the inner surface. Not read (AIP); must be cited next to the sphere
  study. Also Karageorghis and Fairweather 1998 (axisymmetric acoustic MFS).
- Kolezas, Fikioris, Roumeliotis, arXiv:2404.07914 **[full: abstract,
  introduction]**: MFS = MAS; convergence of the field while the source
  amplitudes diverge and oscillate when the auxiliary surface does not
  enclose the singularities (2D Laplace, Neumann).

### 7.4 The skeleton theorem has a potential-theory relative

- Gustafsson, "On mother bodies of convex polyhedra", SIAM J. Math. Anal. 29
  (1998) 1106-1117; Gustafsson and Sakai, Geom. Dedicata 76 (1999) 1-30
  **[statement read as quoted by Tran, arXiv:1710.06991, Theorem 1; the
  originals are behind SIAM / Springer]**: a convex polyhedron with uniform
  mass has a unique "mother body", a measure with the same exterior
  Newtonian potential whose support "is contained in a finite union of
  hyperplanes and reaches the boundary only at corners and edges (not at
  faces)". That support is the ridge of the distance function, i.e. the
  medial surface of `DERIVATIONS.md` A. Different problem (the potential of
  a filled body, not a field scattered by a conductor) and a different proof,
  but the same geometric statement. The skeleton observation must cite it,
  and should be presented as the scattering counterpart, not as a new kind
  of result.

### 7.5 The torus

- Not found: any statement of where the continued scattered field of a 3D
  torus is singular. Found: MAS with the auxiliary contour around the foci
  of an elliptic section (7.2); Aberra 2002 (the Schwarz potential of
  even-dimensional tori), described by Lundberg, arXiv:1009.5159 **[full:
  sections 2-4]**, where for the four-dimensional torus the singular part
  sits on the centre circle only; Millar 1983 (the axially symmetric
  Helmholtz equation; general theory, not read).
- The Khavinson-Shapiro conjecture (Lundberg, Conjecture 2.3) is the
  statement that continuations with entire data are singular only where the
  Schwarz potential is; it "holds in the plane and has been shown to hold
  generically in higher dimensions". Relevant background if the torus
  section makes a claim about the singular set.

### 7.6 Other

- Helsing and Perfekt, arXiv:1203.5997 **[abstract, introduction]**:
  integral-equation solver (recursive compressed inverse preconditioning)
  for the electrostatic cube, 5 to 10 digits. Shows what a specialised BIE
  can do on a cube in the static case; no Maxwell scattering.
- Yurkin and Kahnert, JQSRT 123 (2013) 176-183 **[abstract]**: accuracy
  limits of the discrete dipole approximation and the T-matrix method for a
  *dielectric* cube: DDA 1e-7 to 1e-3 after extrapolation, T-matrix 1e-4 to
  0.2. The only cube accuracy study found; different material, so not a
  reference for our numbers, but it sets the scale for "0.011%".
- Wriedt and Comberg 1998 compare DDA, FDTD and the extended boundary
  condition method on a dielectric cube (abstract only); the 3D MMP code
  has been applied to cubes (Scattport listing) without, as far as the
  abstracts show, any special edge treatment.

### 7.7 Verdict

- **Open problem?** Yes in the sense that matters: three papers by the
  people who started the lightning line of work say the 3D edge and vertex
  case is unsolved or doubtful (2019, 2022, 2026), and the MAS literature
  avoids it by rounding. Nothing found contradicts that.
- **Novel?** The combination is: (i) sources on the skeleton of a 3D
  conductor, as lines along the edges with geometrically clustered depths;
  (ii) the singular tail with exponents taken from the edge, vertex and tip
  expansions; (iii) the same at a conical tip of a smooth thin body;
  (iv) validation at the 1e-4 level against a graded BEM and a public
  benchmark. Not novel: the enclosure rule, exponential clustering in 2D,
  singular functions in 2D MFS, MFS for Maxwell on smooth bodies, tori as
  test bodies, the skeleton as the support of an equivalent source (in
  potential theory).
- **Limits of this search.** Web and arXiv searches plus full texts of what
  is openly hosted. Publisher-only papers were seen as abstracts at best
  (section 8). Russian and Georgian MAS and discrete-sources literature of
  the 1980s-90s is only partly online; a precedent there cannot be excluded.
  Safe wording stays "to our knowledge".

## 8. Papers to obtain (publisher access), in order of importance

    Millar. Arch. Ration. Mech. Anal. 81 (1983) 349-372. 10.1007/BF00250860            (torus singular set: is the focal ring covered?)
    Gustafsson. SIAM J. Math. Anal. 29 (1998) 1106-1117. 10.1137/S0036141097317918     (mother bodies; wording of the skeleton claim)
    Kangro. AIP Conf. Proc. 1479 (2012) 2328-2331. 10.1063/1.4756660                   (convergence of MFS for 3D Maxwell)
    Ye, Hu, Xu, Yan. Eng. Anal. Bound. Elem. 155 (2023) 401-406. 10.1016/j.enganabound.2023.05.050
    De Smedt, Van Bladel. IEEE Trans. Antennas Propag. 34 (1986) 865-870. 10.1109/TAP.1986.1143916   (compare exponent tables)
    Kelley, Chamulak, Courtney, Yilmaz. IEEE Antennas Propag. Mag. 62(1) (2020) 96-106. 10.1109/MAP.2019.2955702   (measurement uncertainty)
    Costabel, Dauge. Arch. Ration. Mech. Anal. 151 (2000) 221-276. 10.1007/s002050050197             (edge and corner singularities; check B and D against it)
    Shifman, Friedmann, Leviatan. IEE Proc. Microw. Antennas Propag. 144 (1997) 235. 10.1049/ip-map:19971186
    Tsitsas, Zouros, Fikioris, Leviatan. IEEE Trans. Antennas Propag. 66 (2018) 5443-5452. 10.1109/TAP.2018.2855963
    Yurkin, Kahnert. J. Quant. Spectrosc. Radiat. Transf. 123 (2013) 176-183. 10.1016/j.jqsrt.2012.10.001
    Young, Ruan. CMES 7 (2005) 223-232.                                                (first MFS for Maxwell; DOI not in Crossref)
    Aberra. Complex Var. Theory Appl. 47 (2002) 1-15. 10.1080/02781070290002895
    Hafner, Bomholt. The 3D Electromagnetic Wave Simulator. Wiley 1993.                (how MMP treats cubes)

Added and verified against Crossref in this pass:

    Hochman, Leviatan, White. J. Comput. Phys. 238 (2013) 337-358. 10.1016/j.jcp.2012.08.015
    Gustafsson, Sakai. Geom. Dedicata 76 (1999) 1-30. 10.1023/A:1005184009159
    Karageorghis, Fairweather. J. Acoust. Soc. Am. 104 (1998) 3212-3218.              (as cited by Liu and Barnett; not checked)
    Kolezas, Fikioris, Roumeliotis. arXiv:2404.07914 (2024)
    Kouroublakis, Tsitsas, Leviatan. arXiv:2607.28558 (2026)
    Lundberg. arXiv:1009.5159 (2010)
    Tran. arXiv:1710.06991 (2017)

## 9. Bibliography file of the draft (2026-11-17)

`paper/refs.bib` holds the 52 entries the draft cites. `paper/check_bib.py`
fetched the Crossref record of every DOI and compared titles; all agree
after these corrections to what had been typed from memory:

    Tabatadze, İşcan, Karaçuha 2025: "On electromagnetic scattering from perfectly conducting parallel square plates and their resonance characteristics using method of auxiliary sources"
    Millar 1983: "The analytic continuation of solutions of the generalized axially symmetric Helmholtz equation"   (so it is about exactly the torus question; still not read)
    Woo et al. 1993: "EM programmer's notebook - benchmark radar targets for the validation of computational electromagnetics programs"
    Kelley et al. 2020: "Rye Canyon radar cross-section measurements of benchmark almond targets [EM programmer's notebook]"
    Ye, Hu, Xu, Yan 2023: "A meshless regularized method of fundamental solution for electromagnetic scattering problems of three-dimensional perfect electric conductor targets" (given names Lv, ShengBo, GuangZhao, TingTing)

arXiv records (`paper/check_arxiv.py`): Ginn and Trefethen, "Lightning
Helmholtz solver" (2310.01665); Kong and Shen, "Lightning algorithm for
inhomogeneous modified Helmholtz equations on cornered domains and their
applications" (2610.26540); Weed, Wu, Huang, Cho, "Method of fundamental
solutions for Maxwell's equations in bi-periodic multilayered media"
(2605.15527); Cai, Bagge, Martinsson, "A stable and fast method for solving
multibody scattering problems via the method of fundamental solutions"
(2603.19113).

New in this session: Björck, "Stability analysis of the method of seminormal
equations for linear least squares problems", Linear Algebra Appl. 88-89
(1987) 31-48, 10.1016/0024-3795(87)90101-7 (Crossref record matched), for the
solver used on the cube and the almond.

Without a verifiable record, and marked so in the file: Young and Ruan 2005
(not in Crossref); Karkashadze and Zaridze's MAS review (author copy at
lae.tsu.ge; year and publisher unknown); the books of Doicu, Eremin and
Wriedt, Kirsch and Hettlich, Colton and Kress, Golub and Van Loan (no DOI
entered; editions as cited in the thesis).
