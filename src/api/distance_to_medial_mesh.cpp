#include "libmat/distance_to_medial_mesh.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

#include "api/internal/cuda_backend.h"
#include "api/internal/dist2mat_cpu.h"
#include "api/internal/surface_context.h"

namespace libmat {

void distance_to_medial_mesh(const Eigen::MatrixXd& P,
                             const Eigen::MatrixXd& C,
                             const Eigen::VectorXd& R,
                             const Eigen::MatrixXi& E,
                             const Eigen::MatrixXi& F, bool use_cuda,
                             Eigen::VectorXd& D, Eigen::VectorXi& prim_type,
                             Eigen::VectorXi& prim_id) {
  internal::check_cols(P, 3, "P");
  internal::check_cols(C, 3, "C");
  if (R.size() != C.rows())
    throw std::invalid_argument("R must have one entry per row of C");
  if (E.size() > 0) internal::check_cols(E, 2, "E");
  if (F.size() > 0) internal::check_cols(F, 3, "F");
  internal::check_indices(E, C.rows(), "E");
  internal::check_indices(F, C.rows(), "F");
  if (C.rows() == 0) throw std::invalid_argument("C is empty");

  if (use_cuda) {
    internal::distance_to_medial_mesh_cuda(P, C, R, E, F, D, prim_type,
                                           prim_id);
    return;
  }

  std::vector<internal::Sphere> spheres(C.rows());
  for (int i = 0; i < C.rows(); ++i)
    spheres[i] << C(i, 0), C(i, 1), C(i, 2), R(i);
  const int ne = E.size() > 0 ? E.rows() : 0;
  const int nf = F.size() > 0 ? F.rows() : 0;

  D.resize(P.rows());
  prim_type.resize(P.rows());
  prim_id.resize(P.rows());
  auto work = [&](int begin, int end) {
    for (int i = begin; i < end; ++i) {
      const Eigen::Vector3d p = P.row(i).transpose();
      double best = std::numeric_limits<double>::infinity();
      int best_type = -1, best_id = -1;
      auto update = [&](double d, int type, int id) {
        if (d < best) best = d, best_type = type, best_id = id;
      };
      for (int s = 0; s < (int)spheres.size(); ++s)
        update(internal::distance_to_sphere(p, spheres[s]), 0, s);
      for (int e = 0; e < ne; ++e)
        update(internal::distance_to_cone(p, spheres[E(e, 0)],
                                          spheres[E(e, 1)]),
               1, e);
      for (int f = 0; f < nf; ++f)
        update(internal::distance_to_slab(p, spheres[F(f, 0)],
                                          spheres[F(f, 1)], spheres[F(f, 2)]),
               2, f);
      D(i) = best;
      prim_type(i) = best_type;
      prim_id(i) = best_id;
    }
  };

  const int n = P.rows();
  const int num_threads = std::max(
      1, std::min<int>(std::thread::hardware_concurrency(), n / 64 + 1));
  std::vector<std::thread> threads;
  const int chunk = (n + num_threads - 1) / num_threads;
  for (int t = 0; t < num_threads; ++t) {
    const int begin = t * chunk, end = std::min(n, begin + chunk);
    if (begin >= end) break;
    threads.emplace_back(work, begin, end);
  }
  for (auto& th : threads) th.join();
}

}  // namespace libmat
