---
title: slab_tangent_triangles
---

::: libmat.slab_tangent_triangles.slab_tangent_triangles
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/slab_tangent_triangles.h>

// For every medial face (slab) spanned by three spheres, compute the two
// triangles that are tangent to all three spheres. Together with the medial
// cones they form the boundary of the medial slab, i.e. the reconstructed
// surface of a medial mesh.
//
// Inputs:
//   C  #V by 3 sphere centers
//   R  #V sphere radii
//   F  #F by 3 medial faces
// Outputs:
//   T0     #F by 9 tangent triangle on the side of the face normal (right-hand
//          rule over F), as (x0 y0 z0 x1 y1 z1 x2 y2 z2); corner k touches
//          sphere F(i, k), so T0 is oriented like F (normal pointing away
//          from the slab)
//   T1     #F by 9 tangent triangle on the other side, corner k touches
//          sphere F(i, k) (flip corners 1 and 2 for an outward normal)
//   valid  #F 1 if the tangent triangles exist, 0 otherwise (e.g. one sphere
//          contains another one); invalid rows are filled with zeros
void slab_tangent_triangles(const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
                            const Eigen::MatrixXi& F, Eigen::MatrixXd& T0,
                            Eigen::MatrixXd& T1, Eigen::VectorXi& valid);
```
