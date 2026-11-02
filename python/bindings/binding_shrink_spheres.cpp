#include "common.h"
#include "libmat/shrink_spheres.h"

void bind_shrink_spheres(py::module_& m) {
  m.def("shrink_spheres",
        [](const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int num_spheres,
           bool preserve_features, double thres_concave, double thres_convex,
           int max_iter, unsigned int seed) {
          libmat::ShrinkSpheresOptions opt;
          opt.num_spheres = num_spheres;
          opt.preserve_features = preserve_features;
          opt.thres_concave = thres_concave;
          opt.thres_convex = thres_convex;
          opt.max_iter = max_iter;
          opt.seed = seed;
          Eigen::MatrixXd C;
          Eigen::VectorXd R;
          Eigen::VectorXi types;
          {
            libmat_py::CallGuard guard;
            libmat::shrink_spheres(V, F, opt, C, R, types);
          }
          return py::make_tuple(C, R, types);
        });
}
