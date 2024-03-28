#include "../params.h"
#include "Delaunay_psm.h"
#include "convex_cell.h"
#include "voronoi_fallback.h"

void initialize_geogram(int argc, char** argv) {
  GEO::initialize();
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("log:quiet", true);
  GEO::CmdLine::set_arg("algo:delaunay", "BDEL");
  GEO::CmdLine::parse(argc, argv);
  // Use lexicographic order for symbolic perturbations.
  // We do that instead of point address, because our
  // points are local variables, generated on the fly,
  // with randomly varying addresses.
  GEO::PCK::set_SOS_mode(GEO::PCK::SOS_LEXICO);
}

void fallback_voro_diagram_CPU(std::vector<float>& pts,
                               std::vector<Status>& stat,
                               std::vector<float>& bary,
                               std::vector<int>& KNN) {
  int nb_pts = pts.size() / 3;
  geo_assert(pts.size() == nb_pts * 3);
  geo_assert(stat.size() == nb_pts);
  geo_assert(bary.size() == nb_pts * 3);
  geo_assert(KNN.size() == nb_pts * _K_);

  GEO::Delaunay_var delaunay;
  std::vector<double> points(3 * (_K_ + 1));
  VBW::ConvexCell C;
  GEO::vector<GEO::index_t> neighbors;

  int nb_recompute = 0;
  for (int v = 0; v < nb_pts; ++v) {
    if (stat[v] != success) {
      ++nb_recompute;

      // Create a Delaunay triangulation if we do not have
      // one already.
      if (delaunay.is_nil()) {
        delaunay = GEO::Delaunay::create(3);
        //   We use circular incident cell lists to
        // accelerate Voronoi cell constructions.
        delaunay->set_stores_cicl(true);
        //   Keep infinite vertex, we need it to determine
        // infinite Voronoi cells / infinite Voronoi facets.
        delaunay->set_keeps_infinite(true);
      }

      // Copy current point and its neighbors.
      points[0] = double(pts[3 * v]);
      points[1] = double(pts[3 * v + 1]);
      points[2] = double(pts[3 * v + 2]);

      for (int j = 0; j < _K_; ++j) {
        int neigh = KNN[v * _K_ + j];
        points[3 + 3 * j] = double(pts[3 * neigh]);
        points[3 + 3 * j + 1] = double(pts[3 * neigh + 1]);
        points[3 + 3 * j + 2] = double(pts[3 * neigh + 2]);
      }

      // Compute the Delaunay triangulation of the
      // point and its neighbors.
      delaunay->set_vertices(_K_ + 1, points.data());

      // Extract the Voronoi cell of v from the Delaunay
      // triangulation.
      VBW::copy_Voronoi_cell_from_Delaunay(delaunay, 0, C, neighbors);

      // Clip the cell with the box.
      const double cube_edge_len = 1000.0;
      C.clip_by_plane(VBW::make_vec4(1.0, 0.0, 0.0, 0.0));
      C.clip_by_plane(VBW::make_vec4(-1.0, 0.0, 0.0, cube_edge_len));
      C.clip_by_plane(VBW::make_vec4(0.0, 1.0, 0.0, 0.0));
      C.clip_by_plane(VBW::make_vec4(0.0, -1.0, 0.0, cube_edge_len));
      C.clip_by_plane(VBW::make_vec4(0.0, 0.0, 1.0, 0.0));
      C.clip_by_plane(VBW::make_vec4(0.0, 0.0, -1.0, cube_edge_len));
      C.compute_geometry();

      // For debugging purposes, save the cell in .obj file format.
      if (false) {
        char name[500];
        sprintf(name, "cell_%d.obj", v);
        C.save(name);
      }

      // Compute the barycenter of the Voronoi cell.
      VBW::vec3 g = C.barycenter();
      bary[3 * v] = float(g.x);
      bary[3 * v + 1] = float(g.y);
      bary[3 * v + 2] = float(g.z);
    }
  }
  if (nb_recompute != 0) {
    std::cout << "CPU fallback: recomputed " << nb_recompute
              << " cells on the CPU" << std::endl;
  }
}
