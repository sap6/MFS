# Project notes

Last updated: 2026-11-15. All scripts are in `validation/`. Problem throughout:
exterior scattering of a plane wave by a perfectly conducting (PEC) body,
tangential E_sc = -tangential E_inc on the surface, least-squares collocation.

## 0. Scope (from the author)

- The paper must include the **cube** and the **almond** as examples.
- Nothing else from the thesis is fixed. Frequency, body size, incidence
  direction, polarisation, formulation, source layout and which figures appear
  can all change if that gives a better or more standard result (for example,
  the almond at benchmark frequencies and nose-on incidence instead of the
  thesis's quarter-wavelength, broadside case).

## 1. Results so far

| Shape | Method | Best result | Reference |
|---|---|---|---|
| Unit sphere, k = 2pi | dipoles on inner sphere R = 0.5 | 2e-9 (N = 400, <1 s); 2e-14 (N = 1200) | Mie |
| Unit sphere, k = 4pi | same | 1.3e-12 (N = 1200) | Mie |
| Unit cube, k = 2pi | uniform dipoles on inner cube | stalls at 25% | Bempp |
| Unit cube, k = 2pi | lightning-type edge sources | 0.5% (J = 10); 0.8% at J = 12 (stall) | Bempp (~0.1%) |
| NASA almond, 1.19 GHz, nose-on | confocal / shrunk source surface | unstable or 20%+ | - |
| NASA almond, 1.19 GHz, nose-on | skeleton sources on mid-plane | self-converged to 5e-5; 0.3% from Bempp | Bempp 9558 dofs |

Bempp itself: sphere 0.4% with 9750 dofs (307 s); cube 0.12% at h = 0.05
(8514 dofs, 334 s); almond 98 s (4212 dofs) and 275 s (9558 dofs), still
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

## 3. Current best commands

```
python3 pec_sphere_convergence.py && python3 plot_pec_sphere.py
python3 bem_reference.py cube                      # writes bem_cube_h*.npz
python3 pec_cube_lightning3d.py 6,8,10 tag=NAME kinds=em q=0.5 gamma=2.0 col_p=3 col_q=0.5 beta=0.5 sym=1
python3 pec_almond_skeleton.py 1.19 --bem
python3 proto2d_square.py TE ; python3 proto2d_ellipse.py     # 2D prototypes, seconds
```

Cube cost at J = 10: 3126 unknowns, 35k rows, assembly ~90 s, solve ~85 s.
Almond: 979 sources -> 21 s solve; 1597 -> 50 s.

Files: `mie.py` (reference), `bem_reference.py` (Bempp EFIE),
`pec_cube_lightning3d.py` + `line_kernel.py` (current cube solver),
`pec_almond_skeleton.py` (current almond solver), `pec_almond.py` (geometry and
mesh only; its confocal solver is superseded), `pec_cube_lightning.py` (shared
helpers; its solver is superseded), `pec_cube_edge.py` (superseded).

## 4. Task list (one per session)

- [ ] **T0 Clean-up.** Delete one-off result files (`results_pec_cube_lightning3d_*`
      except the final run, `cube_lightning3d_field_*`, quick JSONs), remove
      superseded solvers, move shared helpers into one module.
- [ ] **T1 Cube stall.** Error rises from 0.58% (J = 10) to 0.8% (J = 12).
      Check in 2D what breaks at depth 5e-5 first; candidates: sampling at the
      deepest levels, rank cutoff `cond`, corners. Target 0.1%.
- [ ] **T2 Almond robustness.** The 2246-source run was under-sampled
      (rows/unknowns 1.6, residual 4.5%): refine collocation with the sources.
      Confirm MFS accuracy with a finer Bempp mesh.
- [ ] **T3 Almond, general incidence (optional).** Broadside incidence
      (d = -z, p = y, the thesis case) has one mirror symmetry and needs all six
      dipole components. Not required; do it only if the paper wants more than
      nose-on incidence.
- [ ] **T4 Almond at benchmark frequencies** (7 and 9.92 GHz) and monostatic
      RCS against published data (Woo et al. 1993).
- [ ] **T5 Error-versus-time plots,** MFS against Bempp, for sphere, cube and
      almond. One script per figure.
- [ ] **T6 Literature check.** Is a 3D lightning-type method for Maxwell new?
      Must cite and compare with: method of auxiliary sources (Fikioris et
      al.), discrete sources method (Eremin, Doicu, Wriedt; sources on the
      axis/skeleton are close to our almond layout), generalized multipole
      technique, Gopal & Trefethen lightning solvers, 2023 regularized MFS for
      PEC, arXiv 2605.15527, Barnett & Betcke 2008. Read full texts.
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
