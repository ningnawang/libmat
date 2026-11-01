---
title: Installation
---

# Installation

## Dependencies

| dependency | how it is obtained |
| --- | --- |
| C++17 compiler, CMake ≥ 3.18 | system |
| [CGAL](https://www.cgal.org) ≥ 5 (with Boost, GMP, MPFR) | system: `brew install cgal` / `apt install libcgal-dev` |
| [Eigen](https://eigen.tuxfamily.org) 3.4 | system, downloaded if missing |
| [geogram](https://github.com/BrunoLevy/geogram) 1.7.9 | downloaded and built statically |
| [pybind11](https://github.com/pybind/pybind11) | pip build dependency / downloaded |
| CUDA toolkit ≥ 11 | only for the CUDA build |

## Python package

### CPU version

```bash
git clone https://github.com/ningnawang/libmat.git
cd libmat
pip install .
```

### CUDA version

On a machine with an NVIDIA GPU and the CUDA toolkit:

```bash
LIBMAT_WITH_CUDA=ON pip install .
```

```python
import libmat
libmat.has_cuda()   # True
```

The CUDA build contains everything of the CPU build, plus
[`medial_mesh_from_spheres`](python/medial_mesh_from_spheres.md) and the GPU
path of [`distance_to_medial_mesh`](python/distance_to_medial_mesh.md).
Set `CMAKE_CUDA_ARCHITECTURES` (e.g. `-C cmake.define.CMAKE_CUDA_ARCHITECTURES=86`)
to build for a GPU other than the one of the build machine.

### Optional extras

```bash
pip install ".[viz]"    # polyscope, used by the examples
pip install ".[test]"   # pytest:  python -m pytest tests
pip install ".[docs]"   # mkdocs:  mkdocs serve
```

### Development build

Build the extension in place and point `PYTHONPATH` to the sources:

```bash
cmake -B build -DLIBMAT_BUILD_PYTHON=ON
cmake --build build -j
export PYTHONPATH=$PWD/python
python -m pytest tests
```

## C++ library

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release            # CPU only
cmake -B build -DCMAKE_BUILD_TYPE=Release -DLIBMAT_WITH_CUDA=ON
cmake --build build -j
```

| CMake option | default | |
| --- | --- | --- |
| `LIBMAT_WITH_CUDA` | `OFF` | build the CUDA modules (`libmat::cuda`) |
| `LIBMAT_BUILD_PYTHON` | `OFF` | build the python extension module |
| `LIBMAT_BUILD_VIEWER` | `OFF` | build the polyscope viewer in `apps/viewer` |
| `LIBMAT_GEOGRAM_DIR` | empty | use an existing geogram checkout instead of downloading it |
| `LIBMAT_AUTO_DOWNLOAD` | `ON` | download missing dependencies |

The old option names `LIBMAT_IS_USE_GPU` and `LIBMAT_IS_BUILD_TEST` are still
accepted. See [C++ library](cpp/index.md) to use libmat from your own project.
