# Data

The results behind the figures and tables of the paper, written by the
`--run` option of the scripts in `../figures/`. Errors are relative 2-norms of
the scattered field at the test points; times are wall-clock seconds on four
cores (Intel Xeon, 2.1 GHz).

| File | Contents |
|---|---|
| `sphere_convergence.json` | sphere: error, boundary residual and coefficient norm for each source radius R and number of sources N (Figure 1) |
| `sphere_nearfield.npz` | sphere: error against the distance from the surface, MFS and Bempp (Figure 11) |
| `torus_sources.npz` | torus: rows (a_s/a, boundary residual, coefficient norm, optical mismatch) for a/R = 0.25, 0.5, 0.8 (Figure 2) |
| `edge_residual.npz` | cube, J = 6: tangential residual near an edge, with and without tails (Figure 5) |
| `cube_convergence.npz` | cube, standard wave: for each layout family (`dense`: two depths per octave, `tail`: one, `notail`, `plain`) rows (J, time, unknowns of both classes, equations, optical mismatch) and the field `<family>_E<J>` at the 400 test points `pts`; the reference `ref_E` (J = 12) (Figure 6, Table 2) |
| `cube_bempp_uniform.npz` | cube: Bempp on quasi-uniform meshes, no symmetry: unknowns, time, field (Figure 10d) |
| `cube_cases_<case>.npz` | cube, cases `k1` (standard wave), `k1o` (general wave), `k2`, `k3` (two and three wavelengths): `mfs` rows (1, J, unknowns per class, all unknowns, equations, time, error, optical mismatch), `bem` (all timed Bempp runs) and `env` (their lower envelope), Bempp rows (grading g, regular and singular quadrature order, intervals per edge, unknowns, time, error) (Figures 7, 10d, Tables 5, C.1) |
| `cube_cases_<case>_ref.npz` | the reference field of each case |
| `cube_scan_<case>.npz` | Bempp: the scan over grading, quadrature and mesh size from which the timed runs were chosen |
| `cost_sphere.npz`, `cost_torus.npz`, `cost_almond.npz` | MFS and Bempp: unknowns, time and error (for the almond the fields, against the finest layout) (Figure 10, Table 4) |
| `almond_rcs_<f>GHz.npz` | almond: monostatic cross-section `sigma_HH`, `sigma_VV` [m²] at 361 look angles `phi_deg` (Figure 9, Table 3) |
| `austin_almond/` | measurements and reference simulation of the Austin RCS Benchmark Suite, problem III-A, unmodified, with their licences |
