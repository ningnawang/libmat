#pragma once
#include <Eigen/Core>

namespace libmat {

// Regular (weighted Delaunay) triangulation of a set of spheres, i.e. the
// dual of their power diagram, computed with CGAL. The sphere with center c
// and radius r has weight r^2.
//
// Inputs:
//   C  #S by 3 sphere centers
//   R  #S sphere radii
// Outputs:
//   T      #T by 4 finite tetrahedra (indices into C)
//   E      #E by 2 finite edges, i.e. pairs of neighboring power cells
//   valid  #S 1 if the sphere is a vertex of the triangulation, 0 if it is
//          hidden (its power cell is empty)
void regular_triangulation(const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
                           Eigen::MatrixXi& T, Eigen::MatrixXi& E,
                           Eigen::VectorXi& valid);

}  // namespace libmat
