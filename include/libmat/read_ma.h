#pragma once
#include <Eigen/Core>
#include <string>

namespace libmat {

// Read a medial mesh from a .ma file (the format of MATFP, MATTopo, Q-MAT):
//
//   #V #E #F
//   v x y z r      (#V lines)
//   e i j          (#E lines)
//   f i j k        (#F lines)
//
// Inputs:
//   path  path to the .ma file
// Outputs:
//   V  #V by 3 medial sphere centers
//   R  #V medial sphere radii
//   E  #E by 2 medial edges (medial cones)
//   F  #F by 3 medial faces (medial slabs)
// Returns true on success.
bool read_ma(const std::string& path, Eigen::MatrixXd& V, Eigen::VectorXd& R,
             Eigen::MatrixXi& E, Eigen::MatrixXi& F);

}  // namespace libmat
