"""Full medial mesh pipeline on the GPU (requires libmat built with CUDA).

tet mesh -> boundary surface -> ball-shrinking spheres -> restricted power
diagram (CUDA) -> medial mesh -> thinning -> .ma file.

    python examples/08_medial_mesh_from_spheres_cuda.py
"""
import sys

import numpy as np

import libmat
from _common import data, output, parse_args, polyscope_or_none

args = parse_args(__doc__, num_spheres=dict(type=int, default=1500),
                  threshold=dict(type=float, default=0.3))

if not libmat.has_cuda():
    print("libmat was built without CUDA; rebuild with LIBMAT_WITH_CUDA=ON:")
    print("  LIBMAT_WITH_CUDA=ON pip install .")
    sys.exit(0)

V, T = libmat.read_tet_mesh(data("lblock.tet"))
SV, SF, _ = libmat.extract_surface(V, T)
C, R, types = libmat.shrink_spheres(SV, SF, num_spheres=args.num_spheres)
print(f"{len(C)} medial spheres")

E, F, importance = libmat.medial_mesh_from_spheres(V, T, C, R)
print(f"medial mesh: {len(E)} edges, {len(F)} faces")
E, F = libmat.thin_medial_mesh(C, R, E, F, importance, args.threshold)
print(f"after thinning: {len(E)} edges, {len(F)} faces")

libmat.write_ma(output("lblock_mat.ma"), C, R, E, F)
print("saved", output("lblock_mat.ma"))

ps = polyscope_or_none(args)
if ps:
    ps.register_surface_mesh("surface", SV, SF, transparency=0.3)
    mm = ps.register_surface_mesh("medial mesh", C, F)
    mm.set_back_face_policy("identical")
    mm.add_scalar_quantity("radius", R, enabled=True)
    ps.show()
