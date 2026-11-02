---
title: shrink_spheres
---

::: libmat.shrink_spheres.shrink_spheres
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/shrink_spheres.h>

struct ShrinkSpheresOptions {
  // number of random surface samples, one sphere is shrunk per sample; at
  // most one sample per triangle, so at most #F spheres
  int num_spheres = 1000;
  // detect sharp/concave features first and let spheres stop at them
  bool preserve_features = true;
  // feature detection thresholds in degrees, see detect_sharp_features()
  double thres_concave = 60.;
  double thres_convex = 55.;
  // maximum number of shrinking iterations per sphere (-1: no limit)
  int max_iter = -1;
  // seed of the random surface sampling
  unsigned int seed = 200;
};

// Sample medial spheres with the ball-shrinking algorithm [Ma et al. 2012]:
// a sphere touching the surface at a random sample p (with inward normal n)
// is shrunk until it is maximal inside the shape. Sphere types are the
// values of SphereType in libmat/matbase/medial_sphere.h.
//
// Inputs:
//   V  #V by 3 vertex positions of a closed triangle mesh
//   F  #F by 3 triangles
//   options  see ShrinkSpheresOptions
// Outputs:
//   C      #S by 3 sphere centers (in the coordinates of V)
//   R      #S sphere radii
//   types  #S sphere types (SphereType)
void shrink_spheres(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                    const ShrinkSpheresOptions& options, Eigen::MatrixXd& C,
                    Eigen::VectorXd& R, Eigen::VectorXi& types);
```
