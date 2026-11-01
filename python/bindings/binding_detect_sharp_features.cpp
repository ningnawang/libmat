#include "common.h"
#include "libmat/detect_sharp_features.h"

void bind_detect_sharp_features(py::module_& m) {
  m.def("detect_sharp_features",
        [](const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
           double thres_concave, double thres_convex) {
          Eigen::MatrixXi E;
          Eigen::VectorXi E_type, E_line, corners, corner_type;
          {
            libmat_py::CallGuard guard;
            libmat::detect_sharp_features(V, F, thres_concave, thres_convex,
                                          E, E_type, E_line, corners,
                                          corner_type);
          }
          return py::make_tuple(E, E_type, E_line, corners, corner_type);
        });
}
