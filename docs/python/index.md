---
title: Python API
---

# Python API

All functions take and return numpy arrays, and are available directly from the `libmat` module:

```python
import libmat
V, F = libmat.read_surface_mesh("data/mesh_69058.obj")
C, R, types = libmat.shrink_spheres(V, F)
```

Every function has a C++ counterpart with the same name in namespace `libmat`, declared in `<libmat/<name>.h>`.

## Input / output

| function | description |
| --- | --- |
| [`read_tet_mesh`](read_tet_mesh.md) | Read a tet mesh (.msh, .tet, .vtk) |
| [`read_surface_mesh`](read_surface_mesh.md) | Read a triangle mesh (.obj, .off, .ply, .stl, ...) |
| [`read_ma`](read_ma.md) | Read a medial mesh (.ma) |
| [`write_ma`](write_ma.md) | Write a medial mesh (.ma) |

## Mesh processing

| function | description |
| --- | --- |
| [`normalize_vertices`](normalize_vertices.md) | Normalize into [0, 1000]^3, the internal frame |
| [`extract_surface`](extract_surface.md) | Boundary surface of a tet mesh |
| [`detect_sharp_features`](detect_sharp_features.md) | Sharp / concave edges, feature lines and corners |

## Medial spheres

| function | description |
| --- | --- |
| [`shrink_spheres`](shrink_spheres.md) | Medial spheres with the ball-shrinking algorithm |
| [`regular_triangulation`](regular_triangulation.md) | Power diagram neighbors of spheres (CGAL) |

## Medial mesh

| function | description |
| --- | --- |
| [`slab_tangent_triangles`](slab_tangent_triangles.md) | Triangles tangent to medial slabs |
| [`distance_to_medial_mesh`](distance_to_medial_mesh.md) | Signed distance to a medial mesh (CPU / CUDA) |
| [`thin_medial_mesh`](thin_medial_mesh.md) | Topology-preserving thinning |
| [`medial_mesh_from_spheres`](medial_mesh_from_spheres.md) :material-expansion-card:{ title="CUDA" } | Medial mesh from spheres via the RPD (CUDA) |

## Utilities

| function | description |
| --- | --- |
| `libmat.has_cuda()` | `True` if libmat was built with CUDA |
| `libmat.set_verbose(bool)` | print the log of the C++ modules (silenced by default) |
| `libmat.SPHERE_TYPES` | names of the sphere types returned by [`shrink_spheres`](shrink_spheres.md) |
| `libmat.__version__` | version string |
