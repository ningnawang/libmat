#include "libmat/normalize_vertices.h"

#include <stdexcept>

#include "api/internal/surface_context.h"

namespace libmat {

void normalize_vertices(const Eigen::MatrixXd& V, double scale,
                        Eigen::MatrixXd& VN, Eigen::RowVector3d& offset,
                        double& factor) {
  internal::check_cols(V, 3, "V");
  if (!(scale > 0)) throw std::invalid_argument("scale must be positive");
  Parameter params;
  params.scale_max = scale;
  const internal::Normalization n =
      internal::normalize_and_fill_params(V, params, VN);
  offset = n.offset;
  factor = n.factor;
}

}  // namespace libmat
