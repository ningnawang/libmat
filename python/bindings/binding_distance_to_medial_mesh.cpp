#include "common.h"
#include "libmat/distance_to_medial_mesh.h"

void bind_distance_to_medial_mesh(py::module_& m) {
  m.def("distance_to_medial_mesh",
        [](const Eigen::MatrixXd& P, const Eigen::MatrixXd& C,
           const Eigen::VectorXd& R, const Eigen::MatrixXi& E,
           const Eigen::MatrixXi& F, bool use_cuda) {
          Eigen::VectorXd D;
          Eigen::VectorXi prim_type, prim_id;
          {
            libmat_py::CallGuard guard;
            libmat::distance_to_medial_mesh(P, C, R, E, F, use_cuda, D,
                                            prim_type, prim_id);
          }
          return py::make_tuple(D, prim_type, prim_id);
        });
}
