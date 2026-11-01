#include "libmat/slab_tangent_triangles.h"

#include <cmath>
#include <stdexcept>

#include "api/internal/surface_context.h"
#include "libmat/matbase/medial_primitives.h"

namespace libmat {

void slab_tangent_triangles(const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
                            const Eigen::MatrixXi& F, Eigen::MatrixXd& T0,
                            Eigen::MatrixXd& T1, Eigen::VectorXi& valid) {
  internal::check_cols(C, 3, "C");
  internal::check_cols(F, 3, "F");
  internal::check_indices(F, C.rows(), "F");
  if (R.size() != C.rows())
    throw std::invalid_argument("R must have one entry per row of C");

  T0 = Eigen::MatrixXd::Zero(F.rows(), 9);
  T1 = Eigen::MatrixXd::Zero(F.rows(), 9);
  valid = Eigen::VectorXi::Zero(F.rows());
  auto vec = [&](int i) { return Vector3(C(i, 0), C(i, 1), C(i, 2)); };
  for (int f = 0; f < F.rows(); ++f) {
    const int a = F(f, 0), b = F(f, 1), c = F(f, 2);
    SimpleTriangle st0, st1;
    const int type = get_triangles_from_three_spheres(
        vec(a), R(a), vec(b), R(b), vec(c), R(c), st0, st1);
    if (type == 1) continue;  // invalid
    valid(f) = 1;

    // get_triangles_from_three_spheres() may swap corners 1 and 2 to orient
    // a triangle; undo it so that corner k is tangent to sphere F(f, k)
    auto tangent_error = [&](const Vector3& p, int i) {
      return std::abs((p - vec(i)).length() - R(i));
    };
    for (SimpleTriangle* st : {&st0, &st1})
      if (tangent_error(st->v[1], b) + tangent_error(st->v[2], c) >
          tangent_error(st->v[1], c) + tangent_error(st->v[2], b))
        std::swap(st->v[1], st->v[2]);

    // T0 is on the side of the face normal (right-hand rule over F(f, :))
    const Vector3 n = GEO::cross(vec(b) - vec(a), vec(c) - vec(a));
    const Vector3 center = (vec(a) + vec(b) + vec(c)) / 3.;
    const Vector3 st0_center = (st0.v[0] + st0.v[1] + st0.v[2]) / 3.;
    const bool is_st0_first = GEO::dot(st0_center - center, n) >= 0;
    const SimpleTriangle& first = is_st0_first ? st0 : st1;
    const SimpleTriangle& second = is_st0_first ? st1 : st0;
    for (int k = 0; k < 3; ++k)
      for (int d = 0; d < 3; ++d) {
        T0(f, 3 * k + d) = first.v[k][d];
        T1(f, 3 * k + d) = second.v[k][d];
      }
  }
}

}  // namespace libmat
