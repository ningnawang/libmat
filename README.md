# LibMAT

A C++ library, with an optional CUDA part and Python bindings, for computing
the **3D medial axis transform** (MAT) of a shape given as a triangle or
tetrahedral mesh.

📖 **Documentation: https://ningnawang.github.io/libmat** (built from `docs/`)

Starter code is from the paper "[MATFP: Computing Medial Axis Transform with Feature Preservation via Restricted Power Diagram](https://github.com/ningnawang/MATFP)".
The extended work "[MATTopo: Topology-preserving Medial Axis Transform with Restricted Power Diagram](https://github.com/ningnawang/mattopo)" relies on this repo using tag **v0.0.1**.

```python
import numpy as np
import libmat

V, F = libmat.read_surface_mesh("data/lblock.obj")
E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(V, F)
C, R, types = libmat.shrink_spheres(V, F, num_spheres=2000)   # medial spheres
T, E_rt, valid = libmat.regular_triangulation(C, R)            # power diagram dual

MV, MR, ME, MF = libmat.read_ma("data/bug.ma")                 # a medial mesh
SV, SF = libmat.read_surface_mesh("data/bug.obj")
D, prim_type, prim_id = libmat.distance_to_medial_mesh(SV, MV, MR, ME, MF)
print("one-sided Hausdorff distance:", np.abs(D).max())
```

## Installation

```bash
python -m pip install libmat
```

installs a prebuilt CPU version (Linux x86_64, macOS Intel / Apple silicon,
Python ≥ 3.10).

To build from source you need a C++17 compiler, CMake ≥ 3.18 and the Boost
headers; CGAL, geogram, Eigen and pybind11 are downloaded if missing.

```bash
pip install .                         # CPU version
LIBMAT_WITH_CUDA=ON pip install .     # + CUDA modules (needs the CUDA toolkit)
```

C++ only:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release [-DLIBMAT_WITH_CUDA=ON]
cmake --build build -j
```

and in your project `target_link_libraries(my_app PRIVATE libmat::libmat)`.
See the [installation page](https://ningnawang.github.io/libmat/installation/)
for all options.

## Functions

Every function exists in Python (`libmat.<name>`, numpy in / out) and in C++
(`libmat::<name>` in `<libmat/<name>.h>`, Eigen in / out).

| function | | |
| --- | --- | --- |
| `read_tet_mesh`, `read_surface_mesh` | read meshes | CPU |
| `read_ma`, `write_ma` | read / write medial meshes (`.ma`) | CPU |
| `normalize_vertices` | MATFP/MATTopo normalization to `[0, 1000]^3` | CPU |
| `extract_surface` | boundary surface of a tet mesh | CPU |
| `detect_sharp_features` | sharp / concave edges, feature lines, corners | CPU |
| `shrink_spheres` | medial spheres (ball-shrinking) | CPU |
| `regular_triangulation` | power diagram neighbors of spheres (CGAL) | CPU |
| `slab_tangent_triangles` | triangles tangent to medial slabs | CPU |
| `distance_to_medial_mesh` | signed distance to a medial mesh | CPU / CUDA |
| `thin_medial_mesh` | topology-preserving thinning | CPU |
| `medial_mesh_from_spheres` | medial mesh from spheres via the restricted power diagram | CUDA |

## Repository layout

```text
include/libmat/      public headers
  <function>.h       functional API (one function per header)
  common/ inputs/ matbase/ matfun/ io/ rpd3d/     CPU modules of MATFP/MATTopo
  cuda/              CUDA modules: rpd3d (GPU RPD), dist2mat, matfun_fix, io
src/                 implementations (src/api: functional API, src/cuda: CUDA)
python/              python package (python/libmat) and pybind11 bindings
examples/            python example scripts
tests/               python tests (pytest)
docs/, mkdocs.yml    documentation website
data/                example meshes
apps/viewer/         polyscope viewer (C++, -DLIBMAT_BUILD_VIEWER=ON)
legacy/              unused code kept for reference (not built)
```

Run the tests with `pip install ".[test]" && python -m pytest tests`, and
the documentation locally with `pip install ".[docs]" && mkdocs serve`.

## Credit

If you use LibMAT in your research, please cite:

```bibtex
@misc{libmat,
  title = {{LibMAT}: A C++ Library for Medial Axis Transform},
  author = {Ningna Wang and others},
  howpublished = "\url{https://github.com/ningnawang/libmat}",
  year = {2024}
}
```

and the papers whose methods you use:

- **MATFP**: Ningna Wang, Bin Wang, Wenping Wang, Xiaohu Guo.
  *Computing Medial Axis Transform with Feature Preservation via Restricted
  Power Diagram*. ACM Transactions on Graphics (SIGGRAPH Asia 2022).
  [Project](https://ningnawang.github.io/projects/2022_matfp/) ·
  [arXiv](https://arxiv.org/abs/2210.13676) ·
  [Code](https://github.com/ningnawang/MATFP)
- **MATTopo**: Ningna Wang, Hui Huang, Shibo Song, Bin Wang, Wenping Wang,
  Xiaohu Guo. *MATTopo: Topology-preserving Medial Axis Transform with
  Restricted Power Diagram*. ACM Transactions on Graphics (SIGGRAPH Asia 2024).
  [Paper](https://doi.org/10.1145/3687763) ·
  [Project](https://ningnawang.github.io/projects/2024_mattopo/) ·
  [arXiv](https://arxiv.org/abs/2403.18761) ·
  [Code](https://github.com/ningnawang/mattopo)
- **Shrinking ball**
  ([`shrink_spheres`](https://ningnawang.github.io/libmat/python/shrink_spheres/)):
  Jaehwan Ma, Sang Won Bae, Sunghee Choi. *3D medial axis point approximation
  using nearest neighbors and the normal field*. The Visual Computer, 2012.
  [Paper](https://doi.org/10.1007/s00371-011-0594-7)

LibMAT is released under the GPL-3.0 license.

## Contribute

We hope LibMAT is useful for your work, and we would love your help to make it
better.

LibMAT is a young library, and we want contributing to it to be as easy as
possible. You can add a new function, improve an existing one, or fix a bug,
and open a pull request. A new function usually consists of a header and a
source file (`include/libmat/<name>.h`, `src/api/<name>.cpp`), a binding
(`python/bindings/`), a Python wrapper with its docstring
(`python/libmat/<name>.py`) and a test (`tests/`).

Please don't let that list stop you: we will gladly review your code with you
and adapt it to the conventions of the library, and the tests will tell us if
anything breaks. If you have code from a past project that is related to the
medial axis or to geometry processing and is not in LibMAT yet, we encourage
you to submit it as it is. Contributors are credited in the documentation of
the functions they write and on the documentation home page.

## Contributors

None yet: be the first!
