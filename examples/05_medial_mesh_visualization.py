"""Visualize a medial mesh (.ma) and the shape it represents.

The shape is the union of medial spheres, cones and slabs. Its boundary is
made of the triangles tangent to the three spheres of every medial face
(plus pieces of cones, not shown here).

    python examples/05_medial_mesh_visualization.py [--ma data/bug.ma]
"""
import numpy as np

import libmat
from _common import data, parse_args, polyscope_or_none

args = parse_args(__doc__, ma=dict(default=data("bug.ma")),
                  surface=dict(default=data("bug.obj")))

V, R, E, F = libmat.read_ma(args.ma)
print(f"medial mesh: {len(V)} spheres, {len(E)} edges, {len(F)} faces")
print(f"Euler characteristic V - E + F = {len(V) - len(E) + len(F)}")

T0, T1, valid = libmat.slab_tangent_triangles(V, R, F)
print(f"{valid.sum()}/{len(F)} slabs have tangent triangles")
# outward oriented triangle soup
tris = np.concatenate([T0[valid], T1[valid][:, [0, 2, 1]]])
env_V = tris.reshape(-1, 3)
env_F = np.arange(len(env_V)).reshape(-1, 3)

ps = polyscope_or_none(args)
if ps:
    mm = ps.register_surface_mesh("medial mesh", V, F)
    mm.set_back_face_policy("identical")
    mm.add_scalar_quantity("radius", R, enabled=True)
    ps.register_curve_network("medial edges", V, E, radius=0.0005, enabled=False)
    ps.register_surface_mesh("slab envelope", env_V, env_F, enabled=False,
                             transparency=0.5)
    SV, SF = libmat.read_surface_mesh(args.surface)
    ps.register_surface_mesh("input surface", SV, SF, enabled=False)
    ps.show()
