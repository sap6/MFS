#!/bin/sh
# Re-times every cost table on this machine, one job after another (never two
# at once). Bempp must have its OpenCL backend: pip install pyopencl pocl-binary-distribution.
# Logs: retime_t_*.txt (ignored by git); data: the fig_cost_*.npz and bem_*.npz files.
set -x
python3 -W ignore fig_cost_sphere.py run > retime_t_sphere.txt 2>&1
python3 -W ignore fig_cost_torus.py run > retime_t_torus.txt 2>&1
python3 -W ignore fig_cost_almond.py run > retime_t_almond.txt 2>&1
python3 -W ignore fig_cost_cube.py run > retime_t_cube.txt 2>&1
python3 -W ignore bem_optical.py 0.1,0.07,0.05,0.04 > retime_t_cube_bem_uniform.txt 2>&1
python3 -W ignore bem_cube_graded.py 8,12,16,20,24,28 2 > retime_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 8,12,16,20,24,28 2 6 8 >> retime_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 12,16,20,24,28 2 8 10 >> retime_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 12,16,20 2 --wave=60,25,35 > retime_t_cube_bem_oblique.txt 2>&1
python3 -W ignore bem_cube_graded.py 16,20,24 2 6 8 --wave=60,25,35 >> retime_t_cube_bem_oblique.txt 2>&1
python3 -W ignore cube_oblique_table.py run > retime_t_cube_oblique.txt 2>&1
# Bempp with the mirror symmetry of the cube (NOTES.md item 55): standard wave, then the oblique wave
python3 -W ignore bem_cube_sym.py 8,12,16,20,24,28 2 > retime_t_cube_bem_sym.txt 2>&1
python3 -W ignore bem_cube_sym.py 12,16,20,24,28,32,40 2 6 8 >> retime_t_cube_bem_sym.txt 2>&1
python3 -W ignore bem_cube_sym.py 8,16,20,24,28,32,40 2 8 10 >> retime_t_cube_bem_sym.txt 2>&1
python3 -W ignore bem_cube_sym.py 12,16,20,24 2 --wave=60,25,35 > retime_t_cube_bem_sym_oblique.txt 2>&1
python3 -W ignore bem_cube_sym.py 16,20,24,28,32 2 6 8 --wave=60,25,35 >> retime_t_cube_bem_sym_oblique.txt 2>&1
python3 -W ignore bem_cube_sym.py 20,24,28,32,40 2 8 10 --wave=60,25,35 >> retime_t_cube_bem_sym_oblique.txt 2>&1
# the cube at two wavelengths (item 57): MFS table, Bempp against it, and the oblique wave
python3 -W ignore cube_wavenumber.py 2 > cube_wavenumber_t_k2.txt 2>&1
python3 -W ignore cube_wavenumber.py 2 bem > cube_wavenumber_t_k2_bem.txt 2>&1
python3 -W ignore cube_wavenumber.py 2 oblique > cube_wavenumber_t_k2_obl.txt 2>&1
(cd ../paper && python3 make_tables.py > make_tables_t_log.txt 2>&1)
echo done
