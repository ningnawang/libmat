#include "libmat/extract_surface.h"

#include "api/internal/surface_context.h"
#include "libmat/io/io.h"

namespace libmat {

void extract_surface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T,
                     Eigen::MatrixXd& SV, Eigen::MatrixXi& SF,
                     Eigen::VectorXi& sf2tet) {
  internal::ensure_geogram_initialized();
  internal::check_cols(V, 3, "V");
  internal::check_cols(T, 4, "T");
  internal::check_indices(T, V.rows(), "T");

  std::vector<float> tet_vertices(3 * V.rows());
  for (int i = 0; i < V.rows(); ++i)
    for (int j = 0; j < 3; ++j) tet_vertices[3 * i + j] = V(i, j);
  std::vector<int> tet_indices(4 * T.rows());
  for (int i = 0; i < T.rows(); ++i)
    for (int j = 0; j < 4; ++j) tet_indices[4 * i + j] = T(i, j);

  std::vector<std::array<float, 3>> surf_vertices;
  std::vector<std::array<int, 3>> surf_faces;
  std::vector<int> mapping;
  get_surface_from_tet(tet_vertices, tet_indices, surf_vertices, surf_faces,
                       mapping);

  // copy positions from V to keep double precision
  sf2tet.resize(mapping.size());
  SV.resize(mapping.size(), 3);
  for (int i = 0; i < (int)mapping.size(); ++i) {
    sf2tet(i) = mapping[i];
    SV.row(i) = V.row(mapping[i]);
  }
  SF.resize(surf_faces.size(), 3);
  for (int i = 0; i < (int)surf_faces.size(); ++i)
    for (int j = 0; j < 3; ++j) SF(i, j) = surf_faces[i][j];
}

}  // namespace libmat
