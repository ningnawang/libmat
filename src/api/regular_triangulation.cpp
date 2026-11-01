#include "libmat/regular_triangulation.h"

#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "api/internal/surface_context.h"
#include "libmat/rpd3d/triangulation.h"

namespace libmat {

void regular_triangulation(const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
                           Eigen::MatrixXi& T, Eigen::MatrixXi& E,
                           Eigen::VectorXi& valid) {
  internal::check_cols(C, 3, "C");
  if (R.size() != C.rows())
    throw std::invalid_argument("R must have one entry per row of C");

  // same CGAL regular triangulation as generate_RT_CGAL_given_spheres(), in
  // double precision and without the 8 bounding box points
  std::vector<std::pair<Weighted_point, RVI>> points(C.rows());
  for (int i = 0; i < C.rows(); ++i) {
    points[i].first = Weighted_point(Point_rt(C(i, 0), C(i, 1), C(i, 2)),
                                      R(i) * R(i));
    points[i].second.all_id = i;
  }
  Rt rt(points.begin(), points.end());

  valid = Eigen::VectorXi::Zero(C.rows());
  for (auto vit = rt.finite_vertices_begin(); vit != rt.finite_vertices_end();
       ++vit)
    valid(vit->info().all_id) = 1;

  std::vector<std::array<int, 4>> tets;
  for (auto cit = rt.finite_cells_begin(); cit != rt.finite_cells_end(); ++cit)
    tets.push_back({{cit->vertex(0)->info().all_id,
                     cit->vertex(1)->info().all_id,
                     cit->vertex(2)->info().all_id,
                     cit->vertex(3)->info().all_id}});
  T.resize(tets.size(), 4);
  for (int i = 0; i < (int)tets.size(); ++i)
    for (int j = 0; j < 4; ++j) T(i, j) = tets[i][j];

  std::vector<aint2> edges;
  for (auto eit = rt.finite_edges_begin(); eit != rt.finite_edges_end();
       ++eit) {
    int a = eit->first->vertex(eit->second)->info().all_id;
    int b = eit->first->vertex(eit->third)->info().all_id;
    if (a > b) std::swap(a, b);
    edges.push_back({{a, b}});
  }
  std::sort(edges.begin(), edges.end());
  E.resize(edges.size(), 2);
  for (int i = 0; i < (int)edges.size(); ++i)
    E(i, 0) = edges[i][0], E(i, 1) = edges[i][1];
}

}  // namespace libmat
