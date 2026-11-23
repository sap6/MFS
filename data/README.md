# Data

The results behind the figures and tables of the paper, written by the
`--run` option of the scripts in `../figures/`. Errors are relative 2-norms of
the scattered field at the test points; times are wall-clock seconds on two
cores (Intel Xeon, 2.3 GHz, 8 GB RAM). Some files hold extra keys or records from the
original runs (for example the k = 4π runs in `sphere_convergence.json`); the
scripts use only those listed here. To look inside one:

    import numpy as np
    d = np.load("data/cube_convergence.npz"); print(list(d))

| File | Contents |
|---|---|
| `sphere_convergence.json` | sphere: for each source radius `R` and number of sources `N`, the error `mid` on the sphere of radius 1.5, the boundary residual `bc` and the coefficient norm `alpha_norm` (Figure 1) |
| `sphere_nearfield.npz` | sphere: error against the distance from the surface, MFS and Bempp (Figure 12) |
| `torus_sources.npz` | torus: keys `a0.25`, `a0.5`, `a0.8` (a/R), rows (a_s/a, boundary residual, coefficient norm, optical mismatch) (Figure 2) |
| `edge_residual.npz` | cube, J = 6: tangential residual near an edge, with and without tails (Figure 5) |
| `cube_convergence.npz` | cube, standard wave: for each layout family (`dense`: two depths per octave, `tail`: one, `notail`, `plain`) rows (J, time, unknowns of both classes, equations, optical mismatch) and the field `<family>_E<J>` at the 400 test points `pts`; the reference `ref_E` (J = 12) (Figure 6, Table 2) |
| `cube_bempp_uniform.npz` | cube: Bempp on quasi-uniform meshes, no symmetry: unknowns, time, field, element size (Figures 6, 11d) |
| `cube_cases_<case>.npz` | cube, cases `k1` (standard wave), `k1o` (general wave), `k2`, `k3` (two and three wavelengths): `mfs` rows (depths per octave: 2 for the layout of the paper, 1 for one depth per octave; J, unknowns per class, unknowns of the classes solved, equations, time, error, optical mismatch), `bem` (all timed Bempp runs) and `env` (their lower envelope), Bempp rows (grading g, regular and singular quadrature order, intervals per edge, unknowns, time, error) (Figures 6, 7, 11d, Tables 5, C.1) |
| `cube_cases_<case>_ref.npz` | the reference field `E` of each case at the points `pts`, and `info` (J, time, unknowns per class, equations, optical mismatch, extinction cross-section) |
| `cube_scan_<case>.npz` | Bempp: the scan over grading, quadrature and mesh size from which the timed runs were chosen |
| `tetrahedron_convergence.npz` | regular tetrahedron, standard wave at k = 2π: MFS layout families (`dense`, `one`, `notail`, `plain`) rows (J, time, unknowns per class, equations, optical mismatch, extinction cross-section) and fields `<family>_E<J>`, reference `ref_E` (J = 10), and graded symmetric Bempp runs `bem` (Figure 8) |
| `scalene_tetrahedron.npz` | asymmetric convex tetrahedron, zero symmetry at k = 2π: two-depth MFS rows `mfs` (J, time, unknowns, equations, optical mismatch, error vs J = 7), reference `ref_E` (J = 7), and non-symmetric graded Bempp runs `bem` up to 8 GB (Section 5.6) |
| `iterative_check.npz` | matrix-free coarse-bordered block-QR LSMR check on the cube (`cube`: J = 4, 6, 8 at k = 2π and J = 6 at k = 4π) and zero-symmetry asymmetric tetrahedron (`scalene`: J = 4, 5 at k = 2π): preconditioner vs dense-A memory, iteration counts, and relative field differences against direct solves (Section 8) |
| `cost_sphere.npz` | `mfs`, `bem`: rows (unknowns, time, error against Mie) (Figure 11a, Table 4) |
| `cost_torus.npz` | `mfs`: rows (unknowns, time, error, boundary residual, optical mismatch, coefficient norm); `bem`: (unknowns, time, error); `ref`: the finest run (unknowns, residual, optical mismatch, extinction cross-section) (Figure 11b, Table 4) |
| `cost_almond.npz` | `mfs`: rows (sources, unknowns, time, optical mismatch) and the fields `mfs_E<i>` of the layouts 3:0 to 8:5 at 1.19 GHz; `bem`: (unknowns, time) and the fields `bem_E<i>` (Figure 11c, Table 4) |
| `cost_almond35.npz` | nose-on almond at 3.5 GHz (`mfs` layouts 3:0 to 8:5 and symmetric `bem` meshes, Section 7.1) |
| `almond_rcs_<f>GHz.npz` | almond, `f` in `3.5`, `5.125`, `7.0`, `10.25`: monostatic cross-section `sigma_HH`, `sigma_VV` [m²] at 361 look angles `phi_deg` (Figure 10, Table 3) |
| `austin_almond/` | measurements and reference simulation of the Austin RCS Benchmark Suite, problem III-A, unmodified, with their licences |
