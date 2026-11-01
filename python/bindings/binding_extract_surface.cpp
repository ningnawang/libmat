#include "common.h"
#include "libmat/extract_surface.h"

void bind_extract_surface(py::module_& m) {
  m.def("extract_surface", [](const Eigen::MatrixXd& V, const Eigen::MatrixXi& T) {
    Eigen::MatrixXd SV;
    Eigen::MatrixXi SF;
    Eigen::VectorXi sf2tet;
    {
      libmat_py::CallGuard guard;
      libmat::extract_surface(V, T, SV, SF, sf2tet);
    }
    return py::make_tuple(SV, SF, sf2tet);
  });
}
