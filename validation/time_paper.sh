#!/bin/sh
# Every timing of the paper, on one machine, in one go (NOTES.md item 64). The errors do not
# depend on the machine and are not touched; only the seconds are measured again. Runs that are
# compared with each other are timed alternately and keep their smallest time. Run nothing else
# meanwhile (not even a LaTeX build). About three hours on four cores. Afterwards:
#     python3 fig_cost_all.py; python3 fig_cube_cases.py; cd ../paper; python3 make_tables.py; latexmk -pdf main.tex
set -x
uname -a > time_paper_t_machine.txt; lscpu | grep -i "model name\|^CPU(s)" >> time_paper_t_machine.txt; date -u >> time_paper_t_machine.txt
python3 -W ignore cube_duel.py 2 > time_paper_t_cube.txt 2>&1                 # cube, standard wave: MFS and best Bempp runs
python3 -W ignore cube_cases.py k1o 2 > time_paper_t_k1o.txt 2>&1             # general wave
python3 -W ignore cube_cases.py k2 2 > time_paper_t_k2.txt 2>&1               # two wavelengths
python3 -W ignore cube_cases.py k3 2 > time_paper_t_k3.txt 2>&1               # three wavelengths
python3 -W ignore fig_cost_sphere.py run > time_paper_t_sphere.txt 2>&1
python3 -W ignore fig_cost_torus.py run > time_paper_t_torus.txt 2>&1
python3 -W ignore fig_cost_almond.py run > time_paper_t_almond.txt 2>&1
python3 -W ignore bem_optical.py 0.1,0.07,0.05,0.04 > time_paper_t_cube_bem_uniform.txt 2>&1     # Bempp without symmetry, for panel (d) of the cost figure
python3 -W ignore bem_cube_graded.py 8,12,16,20,24,28 2 > time_paper_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 8,12,16,20,24,28 2 6 8 >> time_paper_t_cube_bem_graded.txt 2>&1
python3 -W ignore bem_cube_graded.py 12,16,20,24,28 2 8 10 >> time_paper_t_cube_bem_graded.txt 2>&1
date -u >> time_paper_t_machine.txt
echo done
