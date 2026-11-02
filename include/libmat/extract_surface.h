#pragma once
#include <Eigen/Core>

namespace libmat {

// Extract the boundary triangle mesh of a tetrahedral mesh, with faces
// oriented outwards.
//
// Inputs:
//   V  #V by 3 tet mesh vertex positions
//   T  #T by 4 tetrahedra
// Outputs:
//   SV      #SV by 3 surface vertex positions
//   SF      #SF by 3 surface triangles (indices into SV)
//   sf2tet  #SV map from surface vertex to tet mesh vertex: SV = V(sf2tet,:)
void extract_surface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T,
                     Eigen::MatrixXd& SV, Eigen::MatrixXi& SF,
                     Eigen::VectorXi& sf2tet);

}  // namespace libmat
