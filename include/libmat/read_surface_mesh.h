#pragma once
#include <Eigen/Core>
#include <string>

namespace libmat {

// Read a triangle mesh from any format supported by geogram
// (.obj, .off, .ply, .stl, .mesh, .geogram, ...). Polygonal faces are
// triangulated.
//
// Inputs:
//   path  path to the surface mesh file
// Outputs:
//   V  #V by 3 vertex positions
//   F  #F by 3 triangles (0-based vertex indices)
// Returns true on success.
bool read_surface_mesh(const std::string& path, Eigen::MatrixXd& V,
                       Eigen::MatrixXi& F);

}  // namespace libmat
