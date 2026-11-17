"""Fixed cost of one Bempp assembly with the OpenCL backend: the same small
EFIE matrix assembled four times in one process, then two larger ones.

    python3 bem_overhead.py
"""
import time

import bempp_cl.api as bempp
import numpy as np

import bem_reference as br

K = 2 * np.pi
if __name__ == "__main__":
    print("backend:", bempp.DEFAULT_DEVICE_INTERFACE, bempp.DEFAULT_PRECISION, flush=True)
    for h, rep in ((0.3, 4), (0.14, 3), (0.1, 2)):
        grid = bempp.shapes.sphere(h=h)
        for i in range(rep):
            t0 = time.time()
            div = bempp.function_space(grid, "RWG", 0)
            curl = bempp.function_space(grid, "SNC", 0)
            t1 = time.time()
            efie = bempp.operators.boundary.maxwell.electric_field(div, div, curl, K)
            A = efie.weak_form()
            t2 = time.time()
            print(f"h={h} dofs={div.global_dof_count}: spaces {t1 - t0:.2f}s, assembly {t2 - t1:.2f}s", flush=True)
