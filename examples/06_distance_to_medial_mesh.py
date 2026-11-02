"""Measure how well a medial mesh approximates a surface.

Computes the signed distance from every surface vertex to the volume of the
medial mesh (union of spheres, cones and slabs), on the CPU and, if libmat
was built with CUDA, on the GPU.

    python examples/06_distance_to_medial_mesh.py
"""
import time

import numpy as np

import libmat
from _common import data, parse_args, polyscope_or_none

args = parse_args(__doc__)

SV, SF = libmat.read_surface_mesh(data("bug.obj"))
V, R, E, F = libmat.read_ma(data("bug.ma"))
diag = np.linalg.norm(SV.max(0) - SV.min(0))

devices = ["cpu"] + (["cuda"] if libmat.has_cuda() else [])
for device in devices:
    t = time.time()
    D, prim_type, prim_id = libmat.distance_to_medial_mesh(SV, V, R, E, F,
                                                           device=device)
    print(f"[{device}] {len(SV)} points x {len(V) + len(E) + len(F)} "
          f"primitives in {time.time() - t:.2f}s")
    print(f"[{device}] one-sided Hausdorff distance surface -> MAT: "
          f"{np.abs(D).max():.4f} ({100 * np.abs(D).max() / diag:.3f}% of bbox)")
    print(f"[{device}] closest primitives: "
          f"{np.bincount(prim_type, minlength=3)} (spheres, cones, slabs)")

ps = polyscope_or_none(args)
if ps:
    sm = ps.register_surface_mesh("surface", SV, SF)
    sm.add_scalar_quantity("distance to MAT", np.abs(D), enabled=True,
                           cmap="reds")
    sm.add_scalar_quantity("closest primitive type", prim_type)
    mm = ps.register_surface_mesh("medial mesh", V, F, enabled=False)
    mm.set_back_face_policy("identical")
    ps.show()
