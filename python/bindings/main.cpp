#include "common.h"

#include <cstdio>

#ifdef _WIN32
#include <io.h>
#define LIBMAT_DUP _dup
#define LIBMAT_DUP2 _dup2
#define LIBMAT_CLOSE _close
#define LIBMAT_FILENO _fileno
#define LIBMAT_DEVNULL "NUL"
#else
#include <unistd.h>
#define LIBMAT_DUP dup
#define LIBMAT_DUP2 dup2
#define LIBMAT_CLOSE close
#define LIBMAT_FILENO fileno
#define LIBMAT_DEVNULL "/dev/null"
#endif

#include "libmat/config.h"

namespace libmat_py {

bool& verbose() {
  static bool v = false;
  return v;
}

ScopedSilence::ScopedSilence() {
  if (verbose()) return;
  std::fflush(stdout);
  saved_fd_ = LIBMAT_DUP(LIBMAT_FILENO(stdout));
  FILE* devnull = std::fopen(LIBMAT_DEVNULL, "w");
  if (devnull == nullptr) return;
  LIBMAT_DUP2(LIBMAT_FILENO(devnull), LIBMAT_FILENO(stdout));
  std::fclose(devnull);
}

ScopedSilence::~ScopedSilence() {
  if (saved_fd_ < 0) return;
  std::fflush(stdout);
  LIBMAT_DUP2(saved_fd_, LIBMAT_FILENO(stdout));
  LIBMAT_CLOSE(saved_fd_);
}

}  // namespace libmat_py

PYBIND11_MODULE(_libmat, m) {
  m.doc() = "LibMAT C++ bindings (use the `libmat` package instead)";
  m.attr("__version__") = LIBMAT_VERSION;
  m.def("has_cuda", &libmat::has_cuda,
        "True if libmat was built with LIBMAT_WITH_CUDA=ON");
  m.def(
      "set_verbose", [](bool v) { libmat_py::verbose() = v; }, py::arg("verbose"));
  m.def("get_verbose", []() { return libmat_py::verbose(); });

  bind_detect_sharp_features(m);
  bind_distance_to_medial_mesh(m);
  bind_extract_surface(m);
  bind_io(m);
  bind_medial_mesh_from_spheres(m);
  bind_normalize_vertices(m);
  bind_regular_triangulation(m);
  bind_shrink_spheres(m);
  bind_slab_tangent_triangles(m);
  bind_thin_medial_mesh(m);
}
