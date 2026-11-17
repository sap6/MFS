# Code behind the paper

Everything is run from this folder (`cd validation`). The history of what
was tried, and why the parameters are what they are, is in `../NOTES.md`.

## The four solvers

| Shape | Solver | Lines | What is in it |
|---|---|---|---|
| Sphere | `pec_sphere_convergence.py` | 140 | dipoles on an inner sphere, least squares, Mie reference (`mie.py`) |
| Torus | `pec_torus.py` | 270 | dipoles on an inner torus, solved mode by mode about the axis |
| Almond | `pec_almond_skeleton.py`, `almond_rcs.py` | 500 | dipole rows on the mid-plane, tail at the tip, two mirror planes |
| Cube | `cube_mfs.py` (layout and solver), `pec_cube_lightning3d.py` (sources and kernels), `line_tree.py` (line integrals) | 1400 | edge lines on the bisector planes, tails, three mirror planes |

`lsq.py` holds the least-squares solvers (QR and corrected semi-normal
equations). The cube is the only solver that is long: about half of it is
the quadrature for the line sources and the bookkeeping of the eight parity
classes. `pec_cube_lightning3d.py` still contains the first, slower version
of the solver (`SymSources`, `line_kernel.py`), kept as a cross-check.

Self-check of a fresh installation: `python3 mie.py` (errors near 1e-14).

## Figures and tables of the paper

| In the paper | Script | Data |
|---|---|---|
| Fig. 1 sphere | `pec_sphere_convergence.py`, `fig_sphere.py short` | `results_pec_sphere.json` |
| Fig. 2 torus source radius | `fig_torus_sources.py` | `fig_torus_sources.npz` |
| Fig. 3, 7 layouts | `fig_layouts.py` | none |
| Fig. 4 residual near an edge | `fig_edge_residual.py` | `fig_edge_residual.npz` |
| Table 2, Fig. 5 cube convergence | `fig_cost_cube.py run`, `fig_cube_convergence.py` | `fig_cost_cube.npz` |
| Fig. 6 cube, general wave and two wavelengths | `cube_cases.py`, `fig_cube_cases.py` | `cube_cases_*.npz` |
| Table 3, Fig. 8 almond benchmark | `fig_almond_benchmark.py` | `fig_almond_benchmark_*.npz`, `austin_almond/` |
| Fig. 9, Table 4 cost comparison; Fig. 10 near the surface (`sphere_nearfield.py`) | `fig_cost_sphere.py`, `fig_cost_torus.py`, `fig_cost_almond.py`, `cube_duel.py`, `fig_cost_all.py` | `fig_cost_*.npz` |
| Appendix A exponents | `cone_exponents.py` | none |
| Table C.6 best Bempp runs | `bem_grading_scan.py`, `cube_duel.py` | `bem_grading_scan_*.npz` |

All table bodies are written by `../paper/make_tables.py`.

**Timings.** The errors do not depend on the machine; the seconds do.
`sh time_paper.sh` measures every time of the paper again in one go (about
two hours; run nothing else meanwhile) and the figure scripts and
`make_tables.py` then pick the new times up.

## Boundary element side (Bempp-cl)

`bem_env.py` (OpenCL set-up), `bem_reference.py` (EFIE solve for any mesh),
`bem_optical.py` (cube, quasi-uniform meshes), `bem_cube_graded.py` (cube,
meshes graded to the edges), `bem_cube_sym.py` (the same with the three
mirror planes), `bem_torus.py`, `bem_almond_rcs.py`.

## Everything else

Experiments that settled a question recorded in `../NOTES.md` and are not
needed to reproduce the paper: `proto2d_*.py` (two-dimensional prototypes),
`cube_floor.py`, `cube_k3_tune.py`, `cube_wavenumber.py`,
`cube_oblique_table.py`, `cube_general.py`, `cube_speed.py`,
`cube_residual_map.py`, `check_tail_exponents.py`, `diag_cube.py`,
`torus_singularities.py`, `sphere_scaling.py`, `almond_floor.py`,
`almond_null.py`, `bem_backend.py`, `bem_overhead.py`, `test_line_tree.py`,
`pec_almond.py` (the first almond solver, which failed).
