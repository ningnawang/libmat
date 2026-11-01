"""Sample medial spheres with the ball-shrinking algorithm.

    python examples/03_shrink_spheres.py [--num-spheres 2000]
"""
import numpy as np

import libmat
from _common import data, output, parse_args, polyscope_or_none

args = parse_args(__doc__, num_spheres=dict(type=int, default=2000),
                  mesh=dict(default=data("lblock.obj")))

V, F = libmat.read_surface_mesh(args.mesh)
C, R, types = libmat.shrink_spheres(V, F, num_spheres=args.num_spheres,
                                    preserve_features=True)
print(f"{len(C)} medial spheres, radius in [{R.min():.4f}, {R.max():.4f}]")
for t, n in zip(*np.unique(types, return_counts=True)):
    print(f"  {libmat.SPHERE_TYPES.get(t, t)}: {n}")

libmat.write_ma(output("shrink_spheres.ma"), C, R)
print("saved", output("shrink_spheres.ma"))

ps = polyscope_or_none(args)
if ps:
    ps.register_surface_mesh("surface", V, F, transparency=0.3)
    pc = ps.register_point_cloud("medial spheres", C)
    pc.add_scalar_quantity("radius", R, enabled=True)
    pc.set_point_radius_quantity("radius", autoscale=False)
    ps.register_point_cloud("sphere centers", C, radius=0.002, enabled=False)
    ps.show()
