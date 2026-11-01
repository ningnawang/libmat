---
title: detect_sharp_features
---

::: libmat.detect_sharp_features.detect_sharp_features
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/detect_sharp_features.h>

// Detect sharp (convex) feature edges, concave feature edges and corners of a
// closed triangle mesh, and group the feature edges into feature lines.
//
// Inputs:
//   V  #V by 3 vertex positions
//   F  #F by 3 triangles (closed, manifold)
//   thres_concave  threshold for concave edges (smaller = more sensitive)
//   thres_convex   dihedral angle threshold in degrees for sharp edges
//                  (smaller = more sensitive)
// Outputs:
//   E            #E by 2 feature edges (indices into V)
//   E_type       #E edge type: 1 = sharp (convex) edge, 2 = concave edge
//   E_line       #E id of the feature line each edge belongs to
//   corners      #C corner vertices (indices into V)
//   corner_type  #C corner type: 1 = sharp corner (adjacent to sharp edges),
//                2 = concave corner (adjacent to > 2 concave edges only)
void detect_sharp_features(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                           double thres_concave, double thres_convex,
                           Eigen::MatrixXi& E, Eigen::VectorXi& E_type,
                           Eigen::VectorXi& E_line, Eigen::VectorXi& corners,
                           Eigen::VectorXi& corner_type);
```
