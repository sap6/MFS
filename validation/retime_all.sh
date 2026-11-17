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
python3 -W ignore cube_mfs.py 4,5,6,7,8 60 25 35 > retime_t_cube_oblique.txt 2>&1
python3 -W ignore bem_cube_graded.py 12,16,20 2 --wave=60,25,35 > retime_t_cube_bem_oblique.txt 2>&1
python3 -W ignore bem_cube_graded.py 16,20,24 2 6 8 --wave=60,25,35 >> retime_t_cube_bem_oblique.txt 2>&1
echo done
