#pragma once
// Shared helpers for the python bindings.

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdio>

namespace py = pybind11;

namespace libmat_py {

// The legacy C++ modules print a lot to stdout. Unless libmat.set_verbose(True)
// was called, stdout is redirected to /dev/null while a binding runs.
bool& verbose();

class ScopedSilence {
 public:
  ScopedSilence();
  ~ScopedSilence();

 private:
  int saved_fd_ = -1;
};

// Release the GIL and silence stdout for the duration of a C++ call.
struct CallGuard {
  ScopedSilence silence;
  py::gil_scoped_release release;
};

}  // namespace libmat_py

// each binding_<name>.cpp defines one of these
void bind_detect_sharp_features(py::module_& m);
void bind_distance_to_medial_mesh(py::module_& m);
void bind_extract_surface(py::module_& m);
void bind_io(py::module_& m);
void bind_medial_mesh_from_spheres(py::module_& m);
void bind_normalize_vertices(py::module_& m);
void bind_regular_triangulation(py::module_& m);
void bind_shrink_spheres(py::module_& m);
void bind_slab_tangent_triangles(py::module_& m);
void bind_thin_medial_mesh(py::module_& m);
