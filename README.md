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

Requires a C++17 compiler, CMake ≥ 3.18 and CGAL (`brew install cgal` /
`apt install libcgal-dev`). geogram, Eigen and pybind11 are downloaded if
missing.

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

## Attribution

If you use our library in your research paper, please cite us! You can use
the bibtex block below:

```bibtex
@misc{libmat,
  title = {{LibMAT}: A C++ Library for Medial Axis Transform},
  author = {Ningna Wang},
  note = {https://github.com/ningnawang/libmat},
  year = {2024}
}
```
