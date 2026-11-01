#include "libmat/read_surface_mesh.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <fstream>

#include "api/internal/surface_context.h"

namespace libmat {

bool read_surface_mesh(const std::string& path, Eigen::MatrixXd& V,
                       Eigen::MatrixXi& F) {
  internal::ensure_geogram_initialized();
  if (!std::ifstream(path).good()) return false;
  GEO::Mesh mesh;
  GEO::MeshIOFlags flags;
  flags.reset_element(GEO::MESH_CELLS);
  if (!GEO::mesh_load(path, mesh, flags)) return false;
  if (!mesh.facets.are_simplices()) mesh.facets.triangulate();

  V.resize(mesh.vertices.nb(), 3);
  for (GEO::index_t v = 0; v < mesh.vertices.nb(); ++v) {
    const double* p = mesh.vertices.point_ptr(v);
    V(v, 0) = p[0], V(v, 1) = p[1], V(v, 2) = p[2];
  }
  F.resize(mesh.facets.nb(), 3);
  for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f)
    for (GEO::index_t lv = 0; lv < 3; ++lv)
      F(f, lv) = mesh.facets.vertex(f, lv);
  return true;
}

}  // namespace libmat
