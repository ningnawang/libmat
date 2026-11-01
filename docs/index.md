---
title: LibMAT
---

# LibMAT

**LibMAT** is a C++ library, with an optional CUDA part and Python bindings,
for computing the **3D medial axis transform** (MAT) of a shape given as a
triangle or tetrahedral mesh.

<div class="grid cards" markdown>

-   :material-language-python: **Python**

    ---

    Every function takes and returns numpy arrays, in the style of
    [gpytoolbox](https://gpytoolbox.org).

    [:octicons-arrow-right-24: Python API](python/index.md)

-   :material-language-cpp: **C++**

    ---

    The same functions in namespace `libmat` with Eigen matrices, plus the
    low-level modules they are built on, as CMake targets.

    [:octicons-arrow-right-24: C++ library](cpp/index.md)

-   :material-expansion-card: **CPU and CUDA**

    ---

    The CPU library builds everywhere. Build with `LIBMAT_WITH_CUDA=ON` to add
    the GPU restricted power diagram and GPU distance queries.

    [:octicons-arrow-right-24: Installation](installation.md)

-   :material-play-circle: **Examples**

    ---

    Python scripts covering every function, with
    [polyscope](https://polyscope.run/py/) visualizations.

    [:octicons-arrow-right-24: Examples](examples.md)

</div>

## Quick example

```python
import numpy as np
import libmat

# a closed triangle mesh with sharp and concave edges
V, F = libmat.read_surface_mesh("data/mesh_69058.obj")

# its sharp and concave features
E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(V, F)

# medial spheres with the ball-shrinking algorithm (at most one per triangle)
C, R, types = libmat.shrink_spheres(V, F)

# their neighbors in the power diagram
T, E_rt, valid = libmat.regular_triangulation(C, R)

# how well does a medial mesh approximate a surface?
SV, SF = libmat.read_surface_mesh("data/bug.obj")
MV, MR, ME, MF = libmat.read_ma("data/bug.ma")
D, prim_type, prim_id = libmat.distance_to_medial_mesh(SV, MV, MR, ME, MF)
print("one-sided Hausdorff distance:", np.abs(D).max())
```

## Credit

If LibMAT is helpful in your research, please cite us!

```bibtex
@misc{libmat,
  title = {{LibMAT}: A C++ Library for Medial Axis Transform},
  author = {Ningna Wang and others},
  howpublished = "\url{https://github.com/ningnawang/libmat}",
  year = {2024}
}
```

LibMAT is built upon the following papers:

- Ningna Wang, Rui Xu, Yibo Yin, Zichun Zhong, Taku Komura, Wenping Wang,
  Xiaohu Guo. *MATStruct: High-quality Medial Mesh Computation via
  Structure-aware Variational Optimization*. SIGGRAPH Asia 2025 Conference
  Papers.
  [Project](https://ningnawang.github.io/projects/2025_matstruct/) ·
  [arXiv](https://arxiv.org/abs/2510.10751) ·
  [Code](https://github.com/ningnawang/MATStruct)
- Ningna Wang, Hui Huang, Shibo Song, Bin Wang, Wenping Wang, Xiaohu Guo.
  *MATTopo: Topology-preserving Medial Axis Transform with Restricted Power
  Diagram*. ACM Transactions on Graphics (SIGGRAPH Asia 2024).
  [Paper](https://doi.org/10.1145/3687763) ·
  [Project](https://ningnawang.github.io/projects/2024_mattopo/) ·
  [arXiv](https://arxiv.org/abs/2403.18761) ·
  [Code](https://github.com/ningnawang/mattopo)
- Ningna Wang, Bin Wang, Wenping Wang, Xiaohu Guo. *Computing Medial Axis
  Transform with Feature Preservation via Restricted Power Diagram*. ACM
  Transactions on Graphics (SIGGRAPH Asia 2022).
  [Project](https://ningnawang.github.io/projects/2022_matfp/) ·
  [arXiv](https://arxiv.org/abs/2210.13676) ·
  [Code](https://github.com/ningnawang/MATFP)
- Jaehwan Ma, Sang Won Bae, Sunghee Choi. *3D medial axis point approximation
  using nearest neighbors and the normal field*. The Visual Computer, 2012.
  [Paper](https://doi.org/10.1007/s00371-011-0594-7)

LibMAT is released under the GPL-3.0 license.

## Contribute

We hope LibMAT is useful for your work, and we would love your help to make it
better.

LibMAT is a young library, and we want contributing to it to be as easy as
possible. You can add a new function, improve an existing one, or fix a bug,
and open a pull request on
[GitHub](https://github.com/ningnawang/libmat/pulls). A new function usually
consists of a header and a source file (`include/libmat/<name>.h`,
`src/api/<name>.cpp`), a binding (`python/bindings/`), a Python wrapper with
its docstring (`python/libmat/<name>.py`) and a test (`tests/`).

Please don't let that list stop you: we will gladly review your code with you
and adapt it to the conventions of the library, and the tests will tell us if
anything breaks. If you have code from a past project that is related to the
medial axis or to geometry processing and is not in LibMAT yet, we encourage
you to submit it as it is. Contributors are credited in the documentation of
the functions they write and on this page.

## Contributors

None yet: be the first!
