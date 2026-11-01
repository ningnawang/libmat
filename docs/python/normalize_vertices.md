---
title: normalize_vertices
---

::: libmat.normalize_vertices.normalize_vertices
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/normalize_vertices.h>

// Map vertices into the box [0, scale]^3 (uniform scaling, the longest side
// of the bounding box becomes `scale`). This is the normalization used by
// MATFP/MATTopo (scale = 1000):
//
//   VN = (V - offset) * factor
//
// Inputs:
//   V      #V by 3 vertex positions
//   scale  length of the longest side after normalization
// Outputs:
//   VN      #V by 3 normalized vertex positions
//   offset  1 by 3 translation (minimum corner of the bounding box of V)
//   factor  uniform scaling factor
void normalize_vertices(const Eigen::MatrixXd& V, double scale,
                        Eigen::MatrixXd& VN, Eigen::RowVector3d& offset,
                        double& factor);
```
