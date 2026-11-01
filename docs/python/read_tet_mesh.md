---
title: read_tet_mesh
---

::: libmat.read_tet_mesh.read_tet_mesh
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/read_tet_mesh.h>

// Read a tetrahedral mesh from a .msh (gmsh v2/v4, e.g. from fTetWild),
// .tet or .vtk file.
//
// Inputs:
//   path  path to the tet mesh file
// Outputs:
//   V  #V by 3 vertex positions
//   T  #T by 4 tetrahedra (0-based vertex indices)
// Returns true on success.
bool read_tet_mesh(const std::string& path, Eigen::MatrixXd& V,
                   Eigen::MatrixXi& T);
```
