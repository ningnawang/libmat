#include "common.h"
#include "libmat/slab_tangent_triangles.h"

void bind_slab_tangent_triangles(py::module_& m) {
  m.def("slab_tangent_triangles", [](const Eigen::MatrixXd& C,
                                     const Eigen::VectorXd& R,
                                     const Eigen::MatrixXi& F) {
    Eigen::MatrixXd T0, T1;
    Eigen::VectorXi valid;
    {
      libmat_py::CallGuard guard;
      libmat::slab_tangent_triangles(C, R, F, T0, T1, valid);
    }
    return py::make_tuple(T0, T1, valid);
  });
}
