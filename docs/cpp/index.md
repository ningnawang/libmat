---
title: C++ library
---

# C++ library

## Using libmat in a CMake project

```cmake
add_subdirectory(libmat)          # or FetchContent
target_link_libraries(my_app PRIVATE libmat::libmat)
```

| target | content |
| --- | --- |
| `libmat::core` | CPU library: functional API and all CPU modules |
| `libmat::cuda` | CUDA modules (only with `LIBMAT_WITH_CUDA=ON`) |
| `libmat::libmat` | everything that was built |

`<libmat/config.h>` defines `LIBMAT_WITH_CUDA` (0 or 1) and
`libmat::has_cuda()`.

## Functional API

Every function of the [Python API](../python/index.md) is a C++ function in
namespace `libmat`, declared in its own header (one function per file, like
[libigl](https://libigl.github.io)), taking and returning Eigen matrices:

```cpp
#include <libmat/read_surface_mesh.h>
#include <libmat/shrink_spheres.h>

Eigen::MatrixXd V, C;
Eigen::MatrixXi F;
Eigen::VectorXd R;
Eigen::VectorXi types;
libmat::read_surface_mesh("data/mesh_69058.obj", V, F);

libmat::ShrinkSpheresOptions options;
options.num_spheres = 2000;
libmat::shrink_spheres(V, F, options, C, R, types);
```

`#include <libmat/libmat.h>` includes all of them. The C++ declaration of each
function is shown at the bottom of its page in the Python API.

## Low-level modules

The functional API is a thin layer over lower-level modules, which remain
available for writing full pipelines:

```text
include/libmat/
├── <function>.h        functional API (Eigen in / out)
├── common/             common.h, common_cxx.h, common_geogram.h, predicates.h
├── inputs/             TetMesh, SurfaceMesh, AABBWrapper, Parameter,
│                       sharp / concave feature detection
├── matbase/            MedialSphere, MedialMesh, medial primitives (cones, slabs)
├── matfun/             shrinking, updating (sphere optimization), thinning
├── io/                 .msh/.tet/.vtk, .ma, .sph, surface meshes
├── rpd3d/              regular triangulation (CGAL), ConvexCellHost,
│                       power cell update from RPD output
└── cuda/               (LIBMAT_WITH_CUDA=ON)
    ├── common/         CUDA helpers
    ├── rpd3d/          restricted power diagram on the GPU (RPD3D_GPU), kNN
    ├── dist2mat/       GPU distance to medial primitives
    ├── matfun_fix/     topology / feature / geometry fixes
    └── io/             RPD export (Houdini)
src/                    implementations, same layout
```

| module | main entry points |
| --- | --- |
| `inputs` | `TetMesh`, `SurfaceMesh`, `detect_mark_sharp_features()` |
| `matbase` | `MedialSphere`, `MedialMesh`, `get_triangles_from_three_spheres()` |
| `matfun` | `init_and_shrink()`, `shrink_sphere()`, `iterate_sphere()`, `Thinning::prune()` |
| `io` | `load_tet()`, `get_surface_from_tet()`, `load_mat_clean()`, `export_ma_clean()` |
| `rpd3d` | `generate_RT_CGAL_and_mark_valid_spheres()`, `update_power_cells()` |
| `cuda/rpd3d` | `RPD3D_GPU::calculate()`, `compute_clipped_voro_diagram_GPU()` |
| `cuda/dist2mat` | `compute_closest_dist2mat()` |
| `cuda/matfun_fix` | `fix_topo_by_adding_new_sphere()`, `check_and_fix_external_feature()`, `check_and_fix_mm_geo()` |

!!! note "Conventions of the low-level modules"
    - Inputs are expected in `[0, 1000]^3` (see
      [`normalize_vertices`](../python/normalize_vertices.md)); the functional
      API normalizes internally.
    - Call `GEO::initialize()` once before using them.
    - `MedialMesh` and `RPD3D_GPU` delete the objects they point to
      (`MedialMesh::vertices`, `RPD3D_GPU::tet_mesh`, ...), so pass heap
      allocated objects.

## Viewer

`apps/viewer` is a small polyscope viewer for a tet mesh and a `.ma` file
(`-DLIBMAT_BUILD_VIEWER=ON`):

```bash
./build/apps/viewer/libmat_viewer model.msh model.ma
```
