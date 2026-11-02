#include "libmat/read_tet_mesh.h"

#include "api/internal/surface_context.h"
#include "libmat/io/io.h"

namespace libmat {

bool read_tet_mesh(const std::string& path, Eigen::MatrixXd& V,
                   Eigen::MatrixXi& T) {
  internal::ensure_geogram_initialized();
  std::vector<float> vertices;
  std::vector<int> indices;
  Parameter params;
  if (!load_tet(path, vertices, indices, false /*normalize*/, params))
    return false;
  V.resize(vertices.size() / 3, 3);
  for (int i = 0; i < V.rows(); ++i)
    for (int j = 0; j < 3; ++j) V(i, j) = vertices[3 * i + j];
  T.resize(indices.size() / 4, 4);
  for (int i = 0; i < T.rows(); ++i)
    for (int j = 0; j < 4; ++j) T(i, j) = indices[4 * i + j];
  return true;
}

}  // namespace libmat
