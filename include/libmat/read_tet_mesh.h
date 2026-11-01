#pragma once
#include <Eigen/Core>
#include <string>

namespace libmat {

// Read a tetrahedral mesh from a .msh (gmsh v2/v4, e.g. from fTetWild),
// .tet or .vtk file.
//
// Inputs:
//   path  path to the tet mesh file
// Outputs:
//   V  #V by 3 vertex positions
//   T  #T by 4 tetrahedra (0-based vertex indices)
// Returns true on success.
bool read_tet_mesh(const std::string& path, Eigen::MatrixXd& V,
                   Eigen::MatrixXi& T);

}  // namespace libmat
