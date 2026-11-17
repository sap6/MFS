#!/bin/sh
# Everything of retime_all.sh except the MFS table of the cube for the standard wave and the Bempp
# scan, which cube_duel.py times alternately (NOTES.md item 60). Run it in the same session as
# cube_duel.py, so that all timings of the paper come from one machine.
set -x
python3 -W ignore fig_cost_sphere.py run > retime_t_sphere.txt 2>&1
python3 -W ignore fig_cost_torus.py run > retime_t_torus.txt 2>&1
python3 -W ignore fig_cost_almond.py run > retime_t_almond.txt 2>&1
python3 -W ignore bem_optical.py 0.1,0.07,0.05,0.04 > retime_t_cube_bem_uniform.txt 2>&1
python3 -W ignore bem_cube_graded.py 8,12,16,20,24,28 2 > retime_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 8,12,16,20,24,28 2 6 8 >> retime_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 12,16,20,24,28 2 8 10 >> retime_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_sym.py 8,12,16,20,24,28 2 > retime_t_cube_bem_sym.txt 2>&1
python3 -W ignore bem_cube_sym.py 12,16,20,24,28,32,40 2 6 8 >> retime_t_cube_bem_sym.txt 2>&1
python3 -W ignore bem_cube_sym.py 8,16,20,24,28,32,40 2 8 10 >> retime_t_cube_bem_sym.txt 2>&1
python3 -W ignore bem_cube_graded.py 12,16,20 2 --wave=60,25,35 > retime_t_cube_bem_oblique.txt 2>&1
python3 -W ignore bem_cube_graded.py 16,20,24 2 6 8 --wave=60,25,35 >> retime_t_cube_bem_oblique.txt 2>&1
python3 -W ignore cube_oblique_table.py run > retime_t_cube_oblique.txt 2>&1
python3 -W ignore bem_cube_sym.py 12,16,20,24 2 --wave=60,25,35 > retime_t_cube_bem_sym_oblique.txt 2>&1
python3 -W ignore bem_cube_sym.py 16,20,24,28,32 2 6 8 --wave=60,25,35 >> retime_t_cube_bem_sym_oblique.txt 2>&1
python3 -W ignore bem_cube_sym.py 20,24,28,32,40 2 8 10 --wave=60,25,35 >> retime_t_cube_bem_sym_oblique.txt 2>&1
python3 -W ignore cube_wavenumber.py 2 > cube_wavenumber_t_k2.txt 2>&1
python3 -W ignore cube_wavenumber.py 2 bem > cube_wavenumber_t_k2_bem.txt 2>&1
python3 -W ignore cube_wavenumber.py 2 oblique > cube_wavenumber_t_k2_obl.txt 2>&1
echo done
