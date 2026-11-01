"""Neighborhood of medial spheres: regular triangulation (power diagram dual).

    python examples/04_regular_triangulation.py
"""
import numpy as np

import libmat
from _common import data, parse_args, polyscope_or_none

args = parse_args(__doc__)

V, F = libmat.read_surface_mesh(data("mesh_69058.obj"))
C, R, _ = libmat.shrink_spheres(V, F, num_spheres=1000)
T, E, valid = libmat.regular_triangulation(C, R)
print(f"{len(C)} spheres, {(~valid).sum()} hidden")
print(f"regular triangulation: {len(E)} edges, {len(T)} tets")
degree = np.bincount(E.ravel(), minlength=len(C))
print(f"neighbors per sphere: mean {degree.mean():.1f}, max {degree.max()}")

ps = polyscope_or_none(args)
if ps:
    ps.register_surface_mesh("surface", V, F, transparency=0.2)
    ps.register_point_cloud("centers", C, radius=0.003)
    cn = ps.register_curve_network("power diagram neighbors", C, E, radius=0.0007)
    cn.add_scalar_quantity("degree", degree, defined_on="nodes")
    ps.show()
