"""LibMAT: a library for computing the 3D medial axis transform.

Every function lives in its own module and is re-exported here, e.g.::

    import libmat
    V, F = libmat.read_surface_mesh("data/mesh_69058.obj")
    C, R, types = libmat.shrink_spheres(V, F, num_spheres=500)
"""
from ._libmat import __version__, has_cuda, set_verbose, get_verbose

from .detect_sharp_features import detect_sharp_features
from .distance_to_medial_mesh import distance_to_medial_mesh
from .extract_surface import extract_surface
from .medial_mesh_from_spheres import medial_mesh_from_spheres
from .normalize_vertices import normalize_vertices
from .read_ma import read_ma
from .read_surface_mesh import read_surface_mesh
from .read_tet_mesh import read_tet_mesh
from .regular_triangulation import regular_triangulation
from .shrink_spheres import shrink_spheres, SPHERE_TYPES
from .slab_tangent_triangles import slab_tangent_triangles
from .thin_medial_mesh import thin_medial_mesh
from .write_ma import write_ma

__all__ = [
    "__version__",
    "has_cuda",
    "set_verbose",
    "get_verbose",
    "detect_sharp_features",
    "distance_to_medial_mesh",
    "extract_surface",
    "medial_mesh_from_spheres",
    "normalize_vertices",
    "read_ma",
    "read_surface_mesh",
    "read_tet_mesh",
    "regular_triangulation",
    "shrink_spheres",
    "SPHERE_TYPES",
    "slab_tangent_triangles",
    "thin_medial_mesh",
    "write_ma",
]
