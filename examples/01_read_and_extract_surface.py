"""Read a tet mesh and extract its boundary surface.

    python examples/01_read_and_extract_surface.py
"""
import libmat
from _common import data, parse_args, polyscope_or_none

args = parse_args(__doc__)

V, T = libmat.read_tet_mesh(data("lblock.tet"))
SV, SF, sf2tet = libmat.extract_surface(V, T)
print(f"tet mesh: {len(V)} vertices, {len(T)} tets")
print(f"surface:  {len(SV)} vertices, {len(SF)} triangles")

ps = polyscope_or_none(args)
if ps:
    tm = ps.register_volume_mesh("tet mesh", V, tets=T)
    tm.set_cull_whole_elements(True)
    ps.register_surface_mesh("surface", SV, SF, enabled=False)
    ps.show()
