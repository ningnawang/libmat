#include "common.h"
#include "libmat/normalize_vertices.h"

void bind_normalize_vertices(py::module_& m) {
  m.def("normalize_vertices", [](const Eigen::MatrixXd& V, double scale) {
    Eigen::MatrixXd VN;
    Eigen::RowVector3d offset;
    double factor;
    libmat::normalize_vertices(V, scale, VN, offset, factor);
    return py::make_tuple(VN, offset, factor);
  });
}
