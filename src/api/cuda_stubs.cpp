// Compiled only when LIBMAT_WITH_CUDA=OFF: CUDA-only entry points of the
// functional API throw a clear error.

#include <stdexcept>

#include "api/internal/cuda_backend.h"
#include "libmat/medial_mesh_from_spheres.h"

namespace libmat {

namespace {
[[noreturn]] void throw_no_cuda(const char* name) {
  throw std::runtime_error(
      std::string(name) +
      " requires CUDA, but libmat was built with LIBMAT_WITH_CUDA=OFF");
}
}  // namespace

void medial_mesh_from_spheres(const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                              const Eigen::MatrixXd&, const Eigen::VectorXd&,
                              Eigen::MatrixXi&, Eigen::MatrixXi&,
                              Eigen::VectorXd&) {
  throw_no_cuda("medial_mesh_from_spheres");
}

namespace internal {
void distance_to_medial_mesh_cuda(const Eigen::MatrixXd&,
                                  const Eigen::MatrixXd&,
                                  const Eigen::VectorXd&,
                                  const Eigen::MatrixXi&,
                                  const Eigen::MatrixXi&, Eigen::VectorXd&,
                                  Eigen::VectorXi&, Eigen::VectorXi&) {
  throw_no_cuda("distance_to_medial_mesh(use_cuda=true)");
}
}  // namespace internal

}  // namespace libmat
