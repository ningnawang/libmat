"""Detect sharp edges, concave edges and corners of a CAD-like model.

    python examples/02_detect_sharp_features.py
"""
import numpy as np

import libmat
from _common import data, parse_args, polyscope_or_none

args = parse_args(__doc__)

V, F = libmat.read_surface_mesh(data("lblock.obj"))
E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(
    V, F, thres_concave=0.2, thres_convex=40.0)
print(f"{np.sum(E_type == 1)} sharp edges, {np.sum(E_type == 2)} concave edges")
print(f"{len(np.unique(E_line))} feature lines")
print(f"{np.sum(corner_type == 1)} sharp corners, "
      f"{np.sum(corner_type == 2)} concave corners")

ps = polyscope_or_none(args)
if ps:
    ps.register_surface_mesh("surface", V, F, color=(0.9, 0.9, 0.9))
    sharp = ps.register_curve_network("sharp edges", V, E[E_type == 1],
                                      color=(0.1, 0.3, 0.9), radius=0.004)
    concave = ps.register_curve_network("concave edges", V, E[E_type == 2],
                                        color=(0.9, 0.2, 0.1), radius=0.006)
    lines = ps.register_curve_network("feature lines", V, E, enabled=False)
    lines.add_scalar_quantity("line id", E_line, defined_on="edges",
                              enabled=True, cmap="rainbow")
    ps.register_point_cloud("corners", V[corners], radius=0.01)
    ps.show()
