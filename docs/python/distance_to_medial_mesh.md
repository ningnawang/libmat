---
title: distance_to_medial_mesh
---

::: libmat.distance_to_medial_mesh.distance_to_medial_mesh
    options:
      show_root_heading: true
      heading_level: 1

!!! tip "CPU and CUDA"
    Runs on the CPU in every build, and on the GPU (`device="cuda"`) when libmat is built with `LIBMAT_WITH_CUDA=ON`.

## C++

```cpp
#include <libmat/distance_to_medial_mesh.h>

// Signed distance from query points to the volume represented by a medial
// mesh: the union of all medial spheres, medial cones (edges) and medial
// slabs (faces). The distance is positive outside the union and negative
// inside (inside it is the largest penetration depth over all primitives,
// not an exact signed distance).
//
// Inputs:
//   P  #P by 3 query points
//   C  #V by 3 sphere centers
//   R  #V sphere radii
//   E  #E by 2 medial edges (can be empty)
//   F  #F by 3 medial faces (can be empty)
//   use_cuda  compute on the GPU (requires a build with LIBMAT_WITH_CUDA=ON,
//             computes in single precision), otherwise multi-threaded CPU
// Outputs:
//   D          #P signed distances
//   prim_type  #P type of the closest primitive: 0 sphere, 1 cone, 2 slab
//   prim_id    #P index of the closest primitive in C, E or F
void distance_to_medial_mesh(const Eigen::MatrixXd& P,
                             const Eigen::MatrixXd& C,
                             const Eigen::VectorXd& R,
                             const Eigen::MatrixXi& E,
                             const Eigen::MatrixXi& F, bool use_cuda,
                             Eigen::VectorXd& D, Eigen::VectorXi& prim_type,
                             Eigen::VectorXi& prim_id);
```
