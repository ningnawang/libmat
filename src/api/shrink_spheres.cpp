#include "libmat/shrink_spheres.h"

#include <cstdlib>
#include <stdexcept>

#include "api/internal/surface_context.h"
#include "libmat/matfun/shrinking.h"

namespace libmat {

void shrink_spheres(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                    const ShrinkSpheresOptions& options, Eigen::MatrixXd& C,
                    Eigen::VectorXd& R, Eigen::VectorXi& types) {
  if (options.num_spheres <= 0)
    throw std::invalid_argument("num_spheres must be positive");
  internal::SurfaceContext ctx;
  internal::build_surface_context(V, F, options.preserve_features,
                                  options.thres_concave, options.thres_convex,
                                  ctx);

  std::srand(options.seed);
  std::vector<MedialSphere> spheres;
  init_and_shrink(*ctx.sf_mesh, *ctx.tet_mesh, spheres, options.num_spheres,
                  options.max_iter, false /*is_debug*/);

  int n = 0;
  for (const auto& s : spheres) n += !s.is_deleted;
  C.resize(n, 3);
  R.resize(n);
  types.resize(n);
  int i = 0;
  for (const auto& s : spheres) {
    if (s.is_deleted) continue;
    C.row(i) = ctx.normalization.from_normalized(
        Eigen::RowVector3d(s.center[0], s.center[1], s.center[2]));
    R(i) = s.radius / ctx.normalization.factor;
    types(i) = static_cast<int>(s.type);
    ++i;
  }
}

}  // namespace libmat
