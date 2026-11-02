#include "common.h"
#include "libmat/thin_medial_mesh.h"

void bind_thin_medial_mesh(py::module_& m) {
  m.def("thin_medial_mesh",
        [](const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
           const Eigen::MatrixXi& E, const Eigen::MatrixXi& F,
           const Eigen::VectorXd& importance, double threshold) {
          Eigen::MatrixXi E_out, F_out;
          {
            libmat_py::CallGuard guard;
            libmat::thin_medial_mesh(C, R, E, F, importance, threshold, E_out,
                                     F_out);
          }
          return py::make_tuple(E_out, F_out);
        });
}
