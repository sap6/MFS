"""
Table A.1 and the vertex exponents of Appendix A: the exponents of the field at the
tip of the almond (an elliptic cone with half-angles 34.09 and 12.71 degrees), by the
parities of the potential under y → -y and z → -z, and at a vertex of the cube.
Computed directly (about five minutes); the paper's values are in the comments.

    python figures/cone_exponents.py
"""
from mfs import cone

# Almond tip: nu (electric type, Dirichlet) and tau (magnetic type, Neumann), the two smallest of each class
tip = cone.EllipticCone()
print("parity (y, z)   nu (electric type)     tau (magnetic type)")
for sy, sz in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
    nu, tau = tip.exponents(sy, sz)[:2], tip.exponents(sy, sz, neumann=True)[:2]
    print(f"   ({sy:+d}, {sz:+d})      " + ", ".join(f"{v:.5f}" for v in nu) + "       " + ", ".join(f"{v:.5f}" for v in tau))
# Paper: (+, +) 0.29959, 1.46427 | 1.05604, 1.98580;  (-, +) 1.10372, 2.24730 | 0.96803, 1.96545
#        (+, -) 1.04120, 2.10576 | 0.87905, 1.75454;  (-, -) 2.01295, 3.04849 | 1.97538, 2.90094

# Cube vertex: the exterior of a solid octant
vertex = cone.CubeVertex()
print("cube vertex: nu =", ", ".join(f"{v:.5f}" for v in vertex.exponents()[:2]),
      "  tau =", ", ".join(f"{v:.5f}" for v in vertex.exponents(neumann=True)[:2]))
# Paper: nu = 0.45417, 1.23087; tau = 0.840, 1.20635 (De Smedt and Van Bladel 1986: 0.45410, 0.84108)
