#include "libmat/write_ma.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "api/internal/surface_context.h"

namespace libmat {

bool write_ma(const std::string& path, const Eigen::MatrixXd& V,
              const Eigen::VectorXd& R, const Eigen::MatrixXi& E,
              const Eigen::MatrixXi& F) {
  internal::check_cols(V, 3, "V");
  if (R.size() != V.rows())
    throw std::invalid_argument("R must have one entry per row of V");
  if (E.size() > 0) internal::check_cols(E, 2, "E");
  if (F.size() > 0) internal::check_cols(F, 3, "F");
  internal::check_indices(E, V.rows(), "E");
  internal::check_indices(F, V.rows(), "F");

  std::ofstream out(path);
  if (!out.good()) return false;
  out << V.rows() << " " << E.rows() << " " << F.rows() << "\n";
  out << std::setiosflags(std::ios::fixed) << std::setprecision(15);
  for (int i = 0; i < V.rows(); ++i)
    out << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << " " << R(i)
        << "\n";
  for (int i = 0; i < E.rows(); ++i)
    out << "e " << E(i, 0) << " " << E(i, 1) << "\n";
  for (int i = 0; i < F.rows(); ++i)
    out << "f " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\n";
  return out.good();
}

}  // namespace libmat
