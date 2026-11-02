---
title: write_ma
---

::: libmat.write_ma.write_ma
    options:
      show_root_heading: true
      heading_level: 1

## C++

```cpp
#include <libmat/write_ma.h>

// Write a medial mesh to a .ma file (see read_ma for the format).
//
// Inputs:
//   path  output path (overwritten)
//   V  #V by 3 medial sphere centers
//   R  #V medial sphere radii
//   E  #E by 2 medial edges
//   F  #F by 3 medial faces
// Returns true on success.
bool write_ma(const std::string& path, const Eigen::MatrixXd& V,
              const Eigen::VectorXd& R, const Eigen::MatrixXi& E,
              const Eigen::MatrixXi& F);
```
