#!/bin/sh
# Self-test of the code behind the paper (about five minutes on four cores). Each line prints its own verdict;
# the numbers to expect are in the comments.
set -e
python3 mie.py | tail -3                                  # Mie series against itself: errors near 1e-14
python3 pec_torus.py | tail -9                            # torus: mode by mode against the full matrix (1e-14), known solution, plane wave
python3 test_cube_rotation.py 4 | tail -2                 # cube: four class matrices against eight (1e-15); stored field (1e-15)
python3 test_line_tree.py 5 | tail -1                     # fast line quadrature against the adaptive rule (about 1e-7)
python3 cone_exponents.py | tail -4                       # tip exponents of the almond cone
python3 -W ignore test_bem_symmetry.py 2>&1 | grep -v Warn | tail -7   # Bempp with symmetry against full solves
echo "all checks ran"
