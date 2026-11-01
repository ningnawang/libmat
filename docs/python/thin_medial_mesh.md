---
title: thin_medial_mesh
---

::: libmat.thin_medial_mesh.thin_medial_mesh
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/thin_medial_mesh.h>

// Topology-preserving thinning of a medial mesh: medial tetrahedra (four vertices
// whose four triangles are all medial faces) are removed by collapsing
// tet-face simple pairs in order of increasing face importance, then faces of
// those tets with importance below `threshold` are pruned through their free
// edges (face-edge simple pairs). Pruning simple pairs preserves topology.
//
// Inputs:
//   C           #V by 3 sphere centers
//   R           #V sphere radii
//   E           #E by 2 medial edges
//   F           #F by 3 medial faces
//   importance  #F importance of each face (e.g. the length of the dual
//               power-diagram edge over the average sphere diameter)
//   threshold   faces with importance >= threshold are kept
// Outputs:
//   E_out  #E' by 2 remaining edges (indices into C)
//   F_out  #F' by 3 remaining faces (indices into C)
void thin_medial_mesh(const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
                      const Eigen::MatrixXi& E, const Eigen::MatrixXi& F,
                      const Eigen::VectorXd& importance, double threshold,
                      Eigen::MatrixXi& E_out, Eigen::MatrixXi& F_out);
```
