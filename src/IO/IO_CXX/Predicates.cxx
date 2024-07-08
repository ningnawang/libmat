#include "Predicates.hpp"

#include <geogram/delaunay/delaunay_3d.h>
// #include <igl/predicates/predicates.h>

#define GEO_PREDICATES true
const int Predicates::ORI_POSITIVE;
const int Predicates::ORI_ZERO;
const int Predicates::ORI_NEGATIVE;
const int Predicates::ORI_UNKNOWN;

int Predicates::orient_3d(const Vector3& p1, const Vector3& p2,
                          const Vector3& p3, const Vector3& p4) {
#if GEO_PREDICATES
  const int result =
      -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
  igl::predicates::exactinit();
  auto res = igl::predicates::orient3d(p1, p2, p3, p4);
  Scalar result;
  if (res == igl::predicates::Orientation::POSITIVE)
    result = 1;
  else if (res == igl::predicates::Orientation::NEGATIVE)
    result = -1;
  else
    result = 0;
#endif

  if (result > 0)
    return ORI_POSITIVE;
  else if (result < 0)
    return ORI_NEGATIVE;
  else
    return ORI_ZERO;
}
