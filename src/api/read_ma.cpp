#include "libmat/read_ma.h"

#include <fstream>
#include <sstream>
#include <vector>

namespace libmat {

bool read_ma(const std::string& path, Eigen::MatrixXd& V, Eigen::VectorXd& R,
             Eigen::MatrixXi& E, Eigen::MatrixXi& F) {
  std::ifstream in(path);
  if (!in.good()) return false;
  int nv = 0, ne = 0, nf = 0;
  if (!(in >> nv >> ne >> nf) || nv < 0 || ne < 0 || nf < 0) return false;
  V.resize(nv, 3);
  R.resize(nv);
  E.resize(ne, 2);
  F.resize(nf, 3);
  int iv = 0, ie = 0, iff = 0;
  std::string tag;
  while (in >> tag) {
    if (tag == "v" && iv < nv) {
      if (!(in >> V(iv, 0) >> V(iv, 1) >> V(iv, 2) >> R(iv))) return false;
      ++iv;
    } else if (tag == "e" && ie < ne) {
      if (!(in >> E(ie, 0) >> E(ie, 1))) return false;
      ++ie;
    } else if (tag == "f" && iff < nf) {
      if (!(in >> F(iff, 0) >> F(iff, 1) >> F(iff, 2))) return false;
      ++iff;
    } else {
      // unknown line or extra entries: skip the rest of the line
      std::string rest;
      std::getline(in, rest);
    }
  }
  return iv == nv && ie == ne && iff == nf;
}

}  // namespace libmat
