"""Import before bempp_cl. Tells the OpenCL CPU driver (PoCL, through hwloc)
how much memory and how many cores the machine really has.

On the cloud machines PoCL sizes itself from one NUMA node as the kernel
lists it (4 GB of the 16 GB), and then refuses any matrix above a quarter of
that, rounded up to a power of two: 1 to 4 GB depending on the machine, which
caps Bempp at 8000 to 16000 unknowns in double precision. With the real
memory the cap is 4 GB (16384 unknowns). Has no effect if HWLOC_SYNTHETIC is
already set, or on the Numba backend.
"""
import os

if "HWLOC_SYNTHETIC" not in os.environ:
    gb = int(0.9 * os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2 ** 30)
    os.environ["HWLOC_SYNTHETIC"] = f"node:1(memory={gb}GB) core:{os.cpu_count()} pu:1"
