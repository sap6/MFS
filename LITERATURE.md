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
  nu = 0.45417, tau = 0.8400. The value 0.4542 for the exterior of a cube
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
