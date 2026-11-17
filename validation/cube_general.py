"""Unit cube, any incident plane wave. Kept for the old command line; the
solver is cube_mfs.py (eight parity classes of the three mirror planes).

    python3 cube_general.py J [theta_deg phi_deg psi_deg]
"""
import runpy

from cube_mfs import CLASSES, COL, REFL, Cube, balanced, class_part, classes_of, wave  # noqa: F401

CubeSolver = Cube

if __name__ == "__main__":
    runpy.run_module("cube_mfs", run_name="__main__")
