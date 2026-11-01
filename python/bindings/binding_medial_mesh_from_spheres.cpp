#include "common.h"
#include "libmat/medial_mesh_from_spheres.h"

void bind_medial_mesh_from_spheres(py::module_& m) {
  m.def("medial_mesh_from_spheres",
        [](const Eigen::MatrixXd& V, const Eigen::MatrixXi& T,
           const Eigen::MatrixXd& C, const Eigen::VectorXd& R) {
          Eigen::MatrixXi E, F;
          Eigen::VectorXd importance;
          {
            libmat_py::CallGuard guard;
            libmat::medial_mesh_from_spheres(V, T, C, R, E, F, importance);
          }
          return py::make_tuple(E, F, importance);
        });
}
