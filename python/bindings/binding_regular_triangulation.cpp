#include "common.h"
#include "libmat/regular_triangulation.h"

void bind_regular_triangulation(py::module_& m) {
  m.def("regular_triangulation",
        [](const Eigen::MatrixXd& C, const Eigen::VectorXd& R) {
          Eigen::MatrixXi T, E;
          Eigen::VectorXi valid;
          {
            libmat_py::CallGuard guard;
            libmat::regular_triangulation(C, R, T, E, valid);
          }
          return py::make_tuple(T, E, valid);
        });
}
