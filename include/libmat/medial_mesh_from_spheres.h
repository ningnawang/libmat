#pragma once
#include <Eigen/Core>

namespace libmat {

// [CUDA] Compute the medial mesh connecting a set of medial spheres, as the
// dual of their power diagram restricted to the input tet mesh (RPD): two
// spheres are connected by an edge if their
// restricted power cells share a face, and three spheres span a medial face
// if their restricted power cells share an edge.
//
// Requires a build with LIBMAT_WITH_CUDA=ON, throws std::runtime_error
// otherwise.
//
// Inputs:
//   V  #V by 3 tet mesh vertex positions
//   T  #T by 4 tetrahedra
//   C  #S by 3 medial sphere centers (same coordinates as V)
//   R  #S medial sphere radii
// Outputs:
//   E           #E by 2 medial edges (indices into C)
//   F           #F by 3 medial faces (indices into C)
//   importance  #F face importance (dual RPD edge length over the average
//               sphere diameter), to be used with thin_medial_mesh()
void medial_mesh_from_spheres(const Eigen::MatrixXd& V,
                              const Eigen::MatrixXi& T,
                              const Eigen::MatrixXd& C,
                              const Eigen::VectorXd& R, Eigen::MatrixXi& E,
                              Eigen::MatrixXi& F, Eigen::VectorXd& importance);

}  // namespace libmat
