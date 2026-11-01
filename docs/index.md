---
title: LibMAT
hide:
  - navigation
---

# LibMAT

**LibMAT** is a C++ library, with an optional CUDA part and Python bindings,
for computing the **3D medial axis transform** (MAT) of a shape given as a
triangle or tetrahedral mesh.

It contains the building blocks of
[MATFP](https://github.com/ningnawang/MATFP) (*Computing Medial Axis Transform
with Feature Preservation via Restricted Power Diagram*, SIGGRAPH Asia 2022)
and [MATTopo](https://github.com/ningnawang/mattopo) (*Topology-preserving
Medial Axis Transform with Restricted Power Diagram*): sharp and concave
feature detection, medial sphere sampling, the restricted power diagram on the
GPU, medial mesh extraction, thinning, and error measurement.

<div class="grid cards" markdown>

-   :material-language-python: **Python**

    ---

    Every function takes and returns numpy arrays, in the style of
    [gpytoolbox](https://gpytoolbox.org).

    [:octicons-arrow-right-24: Python API](python/index.md)

-   :material-language-cpp: **C++**

    ---

    The same functions in namespace `libmat` with Eigen matrices, plus the
    low-level MATFP/MATTopo modules, as CMake targets.

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
V, F = libmat.read_surface_mesh("data/lblock.obj")

# features preserved by MATFP / MATTopo
E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(V, F)

# medial spheres with the ball-shrinking algorithm
C, R, types = libmat.shrink_spheres(V, F, num_spheres=2000)

# their neighbors in the power diagram
T, E_rt, valid = libmat.regular_triangulation(C, R)

# how well does a medial mesh approximate a surface?
SV, SF = libmat.read_surface_mesh("data/bug.obj")
MV, MR, ME, MF = libmat.read_ma("data/bug.ma")
D, prim_type, prim_id = libmat.distance_to_medial_mesh(SV, MV, MR, ME, MF)
print("one-sided Hausdorff distance:", np.abs(D).max())
```

## Citation

If you use LibMAT in your research, please cite:

```bibtex
@misc{libmat,
  title = {{LibMAT}: A C++ Library for Medial Axis Transform},
  author = {Ningna Wang},
  note = {https://github.com/ningnawang/libmat},
  year = {2024}
}
```

LibMAT is released under the GPL-3.0 license.
