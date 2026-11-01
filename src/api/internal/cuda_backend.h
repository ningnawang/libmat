#pragma once
// CUDA implementations of functional API entry points. Defined in
// src/cuda/api/ (LIBMAT_WITH_CUDA=ON) or in src/api/cuda_stubs.cpp, which
// throws a std::runtime_error.

#include <Eigen/Core>

namespace libmat {
namespace internal {

void distance_to_medial_mesh_cuda(const Eigen::MatrixXd& P,
                                  const Eigen::MatrixXd& C,
                                  const Eigen::VectorXd& R,
                                  const Eigen::MatrixXi& E,
                                  const Eigen::MatrixXi& F, Eigen::VectorXd& D,
                                  Eigen::VectorXi& prim_type,
                                  Eigen::VectorXi& prim_id);

}  // namespace internal
}  // namespace libmat
