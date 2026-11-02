#pragma once
// CPU (double precision) port of the distance functions of
// src/cuda/dist2mat/dist2mat.cu. Keep both in sync.

#include <Eigen/Core>
#include <cmath>

namespace libmat {
namespace internal {

using Sphere = Eigen::Vector4d;  // (x, y, z, r)

// same semantics as CUDA's clamp(): fmaxf(a, fminf(f, b)), so NaN -> b
inline double clamp_cuda(double f, double a, double b) {
  return std::fmax(a, std::fmin(f, b));
}

inline double distance_to_sphere(const Eigen::Vector3d& p, const Sphere& s) {
  return (p - s.head<3>()).norm() - s[3];
}

inline Sphere bary_lerp(const Sphere& m1, const Sphere& m2, const Sphere& m3,
                        double t1, double t2) {
  return m1 * t1 + m2 * t2 + m3 * (1. - t1 - t2);
}

inline bool solve_quadratic(double A, double B, double C, double* root) {
  root[0] = -1.;
  root[1] = -1.;
  if (A == 0.) {  // degenerated case
    if (B != 0.) {
      root[0] = -C / B;
      root[1] = root[0];
    }
  } else {
    const double delta = B * B - 4. * A * C;
    if (delta < 0.) return false;
    root[0] = (-B - std::sqrt(delta)) / (2. * A);
    root[1] = (-B + std::sqrt(delta)) / (2. * A);
  }
  return true;
}

// distance to the medial cone (convex hull) of two spheres
inline double distance_to_cone(const Eigen::Vector3d& pos, Sphere m1,
                               Sphere m2) {
  if (m1[3] > m2[3]) std::swap(m1, m2);
  const Eigen::Vector3d c1 = m1.head<3>(), c2 = m2.head<3>();
  const double r1 = m1[3], r2 = m2[3];

  // the closest center lies on c(t) = c2 + t (c1 - c2), t in [0, 1]
  const Eigen::Vector3d c21 = c1 - c2;
  const Eigen::Vector3d cq2 = c2 - pos;
  const double A = c21.dot(c21);
  const double D = 2. * c21.dot(cq2);
  const double F = cq2.dot(cq2);
  const double R1 = r1 - r2;

  double t = -(A * D - R1 * R1 * D) -
             std::sqrt((D * D - 4.0 * A * F) * (R1 * R1 - A) * R1 * R1);
  t /= 2. * (A * A - A * R1 * R1);
  t = clamp_cuda(t, 0., 1.);

  const Eigen::Vector3d sp = c2 + t * c21;
  const double r = r2 + t * R1;
  return (pos - sp).norm() - r;
}

// distance to the medial slab (convex hull) of three spheres
inline double distance_to_slab(const Eigen::Vector3d& pos, const Sphere& m1,
                               const Sphere& m2, const Sphere& m3) {
  const Eigen::Vector3d c31 = (m1 - m3).head<3>();
  const Eigen::Vector3d c32 = (m2 - m3).head<3>();
  const Eigen::Vector3d cm3 = m3.head<3>() - pos;

  const double R1 = m1[3] - m3[3];
  const double R2 = m2[3] - m3[3];
  const double A = c31.dot(c31);
  const double B = 2. * c31.dot(c32);
  const double C = c32.dot(c32);
  const double D = 2. * c31.dot(cm3);
  const double E = 2. * c32.dot(cm3);
  const double F = cm3.dot(cm3);
  auto sq = [](double x) { return x * x; };
  auto dist = [&](double t1, double t2) {
    return distance_to_sphere(pos, bary_lerp(m1, m2, m3, t1, t2));
  };

  double t1 = -1., t2 = -1.;
  if (R1 == 0. && R2 == 0.) {
    const double denom = 4. * A * C - B * B;
    t1 = (B * E - 2.0 * C * D) / denom;
    t2 = (B * D - 2.0 * A * E) / denom;
  } else if (R1 != 0. && R2 == 0.) {
    const double H2 = -B / (2. * C);
    const double K2 = -E / (2. * C);
    const double W1 =
        sq(2. * A + B * H2) - 4. * R1 * R1 * (A + B * H2 + C * H2 * H2);
    const double W2 = 2. * (2. * A + B * H2) * (B * K2 + D) -
                      4. * R1 * R1 * (B * K2 + 2. * C * H2 * K2 + D + E * H2);
    const double W3 =
        sq(B * K2 + D) - 4. * R1 * R1 * (C * K2 * K2 + E * K2 + F);
    double root[2];  // {t11, t12}
    solve_quadratic(W1, W2, W3, root);
    const double t21 = H2 * root[0] + K2;
    const double t22 = H2 * root[1] + K2;
    t1 = root[0], t2 = t21;
    if (dist(root[1], t22) < dist(root[0], t21)) t1 = root[1], t2 = t22;
  } else if (R1 == 0. && R2 != 0.) {
    const double H1 = -B / (2. * A);
    const double K1 = -D / (2. * A);
    const double W1 =
        sq(2. * C + B * H1) - 4. * R2 * R2 * (C + B * H1 + A * H1 * H1);
    const double W2 = 2. * (2. * C + B * H1) * (B * K1 + E) -
                      4. * R2 * R2 * (B * K1 + 2. * A * H1 * K1 + E + D * H1);
    const double W3 =
        sq(B * K1 + E) - 4. * R2 * R2 * (A * K1 * K1 + D * K1 + F);
    double root[2];  // {t21, t22}
    solve_quadratic(W1, W2, W3, root);
    const double t11 = H1 * root[0] + K1;
    const double t12 = H1 * root[1] + K1;
    t1 = t11, t2 = root[0];
    if (dist(t12, root[1]) < dist(t11, root[0])) t1 = t12, t2 = root[1];
  } else {
    const double L1 = 2. * A * R2 - B * R1;
    const double L2 = 2. * C * R1 - B * R2;
    const double L3 = E * R1 - D * R2;
    if (L1 == 0. && L2 != 0.) {
      t2 = -L3 / L2;
      const double W1 = 4. * A * A - 4. * R1 * R1 * A;
      const double W2 = 4. * A * (B * t2 + D) - 4. * R1 * R1 * (B * t2 + D);
      const double W3 = sq(B * t2 + D) - (C * t2 * t2 + E * t2 + F);
      double root[2];  // t11, t12
      solve_quadratic(W1, W2, W3, root);
      t1 = root[0];
      if (dist(root[1], t2) < dist(root[0], t2)) t1 = root[1];
    } else if (L1 != 0. && L2 == 0.) {
      t1 = L3 / L1;
      const double W1 = 4. * C * C - 4. * R2 * R2 * C;
      const double W2 = 4. * C * (B * t1 + E) - 4. * R2 * R2 * (B * t1 + E);
      const double W3 = sq(B * t1 + E) - (A * t1 * t1 + D * t1 + F);
      double root[2];  // t21, t22
      solve_quadratic(W1, W2, W3, root);
      t2 = root[0];
      if (dist(t1, root[1]) < dist(t1, root[0])) t2 = root[1];
    } else {
      const double H3 = L2 / L1;
      const double K3 = L3 / L1;
      const double W1 =
          sq(2. * C + B * H3) - 4. * R2 * R2 * (A * H3 * H3 + B * H3 + C);
      const double W2 =
          2. * (2. * C + B * H3) * (B * K3 + E) -
          4. * R2 * R2 * (2. * A * H3 * K3 + B * K3 + D * H3 + E);
      const double W3 =
          sq(B * K3 + E) - 4. * R2 * R2 * (A * K3 * K3 + D * K3 + F);
      double root[2];  // t21, t22
      solve_quadratic(W1, W2, W3, root);
      const double t11 = H3 * root[0] + K3;
      const double t12 = H3 * root[1] + K3;
      t1 = t11, t2 = root[0];
      if (dist(t12, root[1]) < dist(t11, root[0])) t1 = t12, t2 = root[1];
    }
  }

  if ((t1 + t2) < 1. && t1 >= 0. && t1 <= 1. && t2 >= 0. && t2 <= 1.)
    return dist(t1, t2);  // valid solution inside the triangle

  const double dis1 = distance_to_cone(pos, m1, m3);
  const double dis2 = distance_to_cone(pos, m2, m3);
  const double dis3 = distance_to_cone(pos, m1, m2);
  return std::fmin(dis1, std::fmin(dis2, dis3));
}

}  // namespace internal
}  // namespace libmat
