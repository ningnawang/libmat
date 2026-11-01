#include "libmat/common/predicates.h"

#include <geogram/delaunay/delaunay_3d.h>

const int Predicates::ORI_POSITIVE;
const int Predicates::ORI_ZERO;
const int Predicates::ORI_NEGATIVE;
const int Predicates::ORI_UNKNOWN;

int Predicates::orient_3d(const Vector3& p1, const Vector3& p2,
                          const Vector3& p3, const Vector3& p4) {
  // geogram's exact predicate (sign convention flipped to match libigl)
  const int result =
      -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());

  if (result > 0)
    return ORI_POSITIVE;
  else if (result < 0)
    return ORI_NEGATIVE;
  else
    return ORI_ZERO;
}
