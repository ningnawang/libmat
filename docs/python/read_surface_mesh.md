---
title: read_surface_mesh
---

::: libmat.read_surface_mesh.read_surface_mesh
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/read_surface_mesh.h>

// Read a triangle mesh from any format supported by geogram
// (.obj, .off, .ply, .stl, .mesh, .geogram, ...). Polygonal faces are
// triangulated.
//
// Inputs:
//   path  path to the surface mesh file
// Outputs:
//   V  #V by 3 vertex positions
//   F  #F by 3 triangles (0-based vertex indices)
// Returns true on success.
bool read_surface_mesh(const std::string& path, Eigen::MatrixXd& V,
                       Eigen::MatrixXi& F);
```
