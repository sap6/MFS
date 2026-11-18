"""
The Method of Fundamental Solutions for Maxwell's equations: scattering of a
plane wave by perfect conductors (a sphere, a torus, a cube and the NASA
almond), as in the paper by Bandyopadhyay and Betcke.

    from mfs import sphere, torus, cube, almond, mie
    y, alpha = sphere.solve(400)            # 400 dipoles on the sphere of radius 0.5

The Bempp comparison (mfs.bem) needs bempp-cl and is imported on its own.
"""
from . import almond, cone, core, cube, mie, sphere, torus  # noqa: F401
from .core import plane_wave, wave  # noqa: F401
