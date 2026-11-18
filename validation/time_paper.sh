#!/bin/sh
# Every timing of the paper, on one machine, in one go (NOTES.md items 66, 68). The errors do not
# depend on the machine and are not touched; only the seconds are measured again. Runs that are
# compared with each other are timed alternately and keep their smallest time. Run nothing else
# meanwhile (not even a LaTeX build). About six hours on four cores (three and a half without the finest Bempp mesh). Afterwards:
#     python3 fig_cost_all.py; python3 fig_cube_cases.py; python3 fig_cube_convergence.py
#     cd ../paper; python3 make_tables.py; latexmk -pdf main.tex
# Needed once beforehand (fields, not timings): python3 fig_cost_cube.py run;
#     python3 cube_cases.py k1o ref; python3 cube_cases.py k2 ref; python3 cube_cases.py k3 ref
set -x
(uname -a; lscpu | grep -i "model name\|^CPU(s)"; date -u) > time_paper_t_machine.txt
python3 -W ignore cube_duel.py 2 > time_paper_t_cube.txt 2>&1                 # cube, standard wave: MFS and best Bempp runs
python3 -W ignore cube_cases.py k1o 2 > time_paper_t_k1o.txt 2>&1             # general wave
python3 -W ignore cube_cases.py k2 2 > time_paper_t_k2.txt 2>&1               # two wavelengths
python3 -W ignore cube_cases.py k3 2 > time_paper_t_k3.txt 2>&1               # three wavelengths
for c in k1 k1o k2 k3; do python3 -W ignore bem_retime_extra.py $c 2 > time_paper_t_extra_$c.txt 2>&1; done   # Bempp runs the scan had excluded on old times
for c in k1 k1o k2 k3; do python3 -W ignore bem_finest.py $c 2 > time_paper_t_finest_$c.txt 2>&1; done   # Bempp on its finest mesh (41472 unknowns); 2.7 hours
python3 -W ignore fig_cost_sphere.py run > time_paper_t_sphere.txt 2>&1
python3 -W ignore fig_cost_torus.py run > time_paper_t_torus.txt 2>&1
python3 -W ignore fig_cost_almond.py run > time_paper_t_almond.txt 2>&1
python3 -W ignore almond_cost_freq.py 3.5 > time_paper_t_almond35.txt 2>&1      # almond at 3.5 GHz (one sentence in Section 7.1)
python3 -W ignore bem_optical.py 0.1,0.07,0.05,0.04 > time_paper_t_cube_bem_uniform.txt 2>&1     # Bempp as distributed, quasi-uniform meshes
date -u >> time_paper_t_machine.txt
echo done
