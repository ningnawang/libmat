"""Topology-preserving thinning of a medial mesh.

A small medial mesh with a medial tet (4 spheres whose 4 triangles are all
medial faces) attached to a sheet. Thinning removes the tet and the faces of
low importance, without changing the topology.

    python examples/07_thin_medial_mesh.py
"""
import numpy as np

import libmat
from _common import parse_args, polyscope_or_none

args = parse_args(__doc__)

# a 4x4 grid sheet in the plane z = 0 ...
n = 4
xs, ys = np.meshgrid(np.arange(n), np.arange(n))
C = np.stack([xs.ravel(), ys.ravel(), np.zeros(n * n)], axis=1).astype(float)
F = []
for i in range(n - 1):
    for j in range(n - 1):
        a, b, c, d = i * n + j, i * n + j + 1, (i + 1) * n + j, (i + 1) * n + j + 1
        F += [[a, b, d], [a, d, c]]
# ... plus a medial tet on top of the triangle (5, 6, 10)
C = np.vstack([C, [1.3, 1.6, 1.0]])
apex = len(C) - 1
F += [[5, 6, apex], [6, 10, apex], [5, 10, apex]]
F = np.array(F)
R = np.full(len(C), 0.2)

# importance: 1 for the sheet, small for the tet faces above it
importance = np.ones(len(F))
importance[-3:] = [0.1, 0.2, 0.3]
E_out, F_out = libmat.thin_medial_mesh(C, R, None, F, importance, threshold=0.5)



def euler(E, F):
    """V - E + F - T, counting medial tets as (filled) 3-cells."""
    faces = {tuple(sorted(f)) for f in F}
    verts = np.unique(np.concatenate([E.ravel(), F.ravel()]))
    tets = {tuple(sorted(f + (d,))) for f in faces for d in verts
            if d not in f and all(tuple(sorted(g)) in faces
                                  for g in [(f[0], f[1], d), (f[0], f[2], d), (f[1], f[2], d)])}
    return len(verts) - len(E) + len(F) - len(tets), len(tets)


E_in = np.unique(np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [0, 2]]]), axis=1), axis=0)
for name, (E_, F_) in [("before", (E_in, F)), ("after", (E_out, F_out))]:
    chi, num_tets = euler(E_, F_)
    print(f"{name:6s}: {len(F_)} faces, {len(E_)} edges, {num_tets} medial tets, "
          f"Euler characteristic {chi}")

ps = polyscope_or_none(args)
if ps:
    ps.register_surface_mesh("before", C, F, transparency=0.4).set_back_face_policy("identical")
    after = ps.register_surface_mesh("after", C, F_out)
    after.set_back_face_policy("identical")
    after.translate((n + 1.0, 0, 0))
    ps.show()
