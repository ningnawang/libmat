#pragma once
#include <cuda_runtime.h>

#include <set>

#define cuda_check(x) \
  if (x != cudaSuccess) exit(1);

// NOT built-in type
struct __align__(16) float5 {
  float x, y, z, w, h;
};
inline __host__ __device__ float5 make_float5(float x, float y, float z,
                                              float w, float h) {
  float5 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  t.h = h;
  return t;
}
inline __host__ __device__ float5 make_float5(float4 orig, float h) {
  float5 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  t.h = h;
  return t;
}
struct __align__(16) float7 {
  float x, y, z, w, h, g, k;
};
inline __host__ __device__ float7 make_float7(float x, float y, float z,
                                              float w, float h, float g,
                                              float k) {
  float7 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  t.h = h;
  t.g = g;
  t.k = k;
  return t;
}
inline __host__ __device__ float7 make_float7(float4 orig, float h, float g,
                                              float k) {
  float7 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  t.h = h;
  t.g = g;
  t.k = k;
  return t;
}
inline __host__ __device__ uchar4 make_uchar4(uchar3 orig, unsigned char w) {
  uchar4 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = w;
  return t;
}
// remove last element
inline __host__ __device__ uchar3 make_uchar3(uchar4 orig) {
  uchar3 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  return t;
}
inline __host__ __device__ float4 make_float4(float5 orig) {
  float4 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  return t;
}
inline __host__ __device__ float4 make_float4(float7 orig) {
  float4 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  return t;
}
inline __host__ __device__ float5 make_float5(float7 orig) {
  float5 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  t.h = orig.h;
  return t;
}

// ###################  Common   ######################
inline __host__ __device__ float4 minus4(float4 A, float4 B) {
  return make_float4(A.x - B.x, A.y - B.y, A.z - B.z, A.w - B.w);
}
inline __host__ __device__ float3 minus3(float3 A, float3 B) {
  return make_float3(A.x - B.x, A.y - B.y, A.z - B.z);
}
inline __host__ __device__ float4 plus4(float4 A, float4 B) {
  return make_float4(A.x + B.x, A.y + B.y, A.z + B.z, A.w + B.w);
}
inline __host__ __device__ float3 plus3(float4 A, float4 B) {
  return make_float3(A.x + B.x, A.y + B.y, A.z + B.z);
}
inline __host__ __device__ float3 plus3(float3 A, float3 B) {
  return make_float3(A.x + B.x, A.y + B.y, A.z + B.z);
}
inline __host__ __device__ float4 divide4(float4 A, float s) {
  return make_float4(A.x / s, A.y / s, A.z / s, A.w / s);
}
inline __host__ __device__ float3 divide3(float3 A, float s) {
  return make_float3(A.x / s, A.y / s, A.z / s);
}
inline __device__ float dot4(float4 A, float4 B) {
  return A.x * B.x + A.y * B.y + A.z * B.z + A.w * B.w;
}
inline __device__ float dot3(float4 A, float4 B) {
  return A.x * B.x + A.y * B.y + A.z * B.z;
}
inline __device__ float dot3(float3 A, float3 B) {
  return A.x * B.x + A.y * B.y + A.z * B.z;
}
inline __device__ float4 mul3(float s, float4 A) {
  return make_float4(s * A.x, s * A.y, s * A.z, 1.);
}
inline __device__ float4 mul4(float s, float4 A) {
  return make_float4(s * A.x, s * A.y, s * A.z, s * A.w);
}
inline __device__ float4 cross3(float4 A, float4 B) {
  return make_float4(A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z,
                     A.x * B.y - A.y * B.x, 0);
}
inline __host__ __device__ float3 cross3(float3 A, float3 B) {
  return make_float3(A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z,
                     A.x * B.y - A.y * B.x);
}
inline __device__ float4 plane_from_point_and_normal(float4 P, float4 n) {
  return make_float4(n.x, n.y, n.z, -dot3(P, n));
}
inline __device__ float4 plane_from_point_and_normal(float3 P, float3 n) {
  return make_float4(n.x, n.y, n.z, -dot3(P, n));
}
inline __host__ __device__ float det2x2(float a11, float a12, float a21,
                                        float a22) {
  return a11 * a22 - a12 * a21;
}
inline __host__ __device__ float det3x3(float a11, float a12, float a13,
                                        float a21, float a22, float a23,
                                        float a31, float a32, float a33) {
  return a11 * det2x2(a22, a23, a32, a33) - a21 * det2x2(a12, a13, a32, a33) +
         a31 * det2x2(a12, a13, a22, a23);
}

inline __device__ float det4x4(float a11, float a12, float a13, float a14,
                               float a21, float a22, float a23, float a24,
                               float a31, float a32, float a33, float a34,
                               float a41, float a42, float a43, float a44) {
  float m12 = a21 * a12 - a11 * a22;
  float m13 = a31 * a12 - a11 * a32;
  float m14 = a41 * a12 - a11 * a42;
  float m23 = a31 * a22 - a21 * a32;
  float m24 = a41 * a22 - a21 * a42;
  float m34 = a41 * a32 - a31 * a42;

  float m123 = m23 * a13 - m13 * a23 + m12 * a33;
  float m124 = m24 * a13 - m14 * a23 + m12 * a43;
  float m134 = m34 * a13 - m14 * a33 + m13 * a43;
  float m234 = m34 * a23 - m24 * a33 + m23 * a43;

  return (m234 * a14 - m134 * a24 + m124 * a34 - m123 * a44);
}

inline __device__ double det2x2(double a11, double a12, double a21,
                                double a22) {
  return a11 * a22 - a12 * a21;
}

inline __device__ double det3x3(double a11, double a12, double a13, double a21,
                                double a22, double a23, double a31, double a32,
                                double a33) {
  return a11 * det2x2(a22, a23, a32, a33) - a21 * det2x2(a12, a13, a32, a33) +
         a31 * det2x2(a12, a13, a22, a23);
}

inline __device__ double det4x4(double a11, double a12, double a13, double a14,
                                double a21, double a22, double a23, double a24,
                                double a31, double a32, double a33, double a34,
                                double a41, double a42, double a43,
                                double a44) {
  double m12 = a21 * a12 - a11 * a22;
  double m13 = a31 * a12 - a11 * a32;
  double m14 = a41 * a12 - a11 * a42;
  double m23 = a31 * a22 - a21 * a32;
  double m24 = a41 * a22 - a21 * a42;
  double m34 = a41 * a32 - a31 * a42;

  double m123 = m23 * a13 - m13 * a23 + m12 * a33;
  double m124 = m24 * a13 - m14 * a23 + m12 * a43;
  double m134 = m34 * a13 - m14 * a33 + m13 * a43;
  double m234 = m34 * a23 - m24 * a33 + m23 * a43;

  return (m234 * a14 - m134 * a24 + m124 * a34 - m123 * a44);
}

// compute the power distance from p to seed (pos, radius)
inline __device__ float get_power_dist(float3 p, float4 seed) {
  float3 pos_seed = make_float3(seed.x, seed.y, seed.z);
  float3 diff = minus3(p, pos_seed);
  float pd = dot3(diff, diff) - seed.w * seed.w;
  return pd;
}

inline __device__ float get_tet_volume(float4 A, float4 B, float4 C) {
  return -det3x3(A.x, A.y, A.z, B.x, B.y, B.z, C.x, C.y, C.z) / 6.;
}
inline __device__ void get_tet_volume_and_barycenter(float4& bary,
                                                     float& volume, float4 A,
                                                     float4 B, float4 C,
                                                     float4 D) {
  volume = get_tet_volume(minus4(A, D), minus4(B, D), minus4(C, D));
  bary = make_float4(.25f * (A.x + B.x + C.x + D.x),
                     .25f * (A.y + B.y + C.y + D.y),
                     .25f * (A.z + B.z + C.z + D.z), 1.0f);
}
inline __device__ float4 project_on_plane(float4 P, float4 plane) {
  float4 n = make_float4(plane.x, plane.y, plane.z, 0);
  float n_2 = dot4(n, n);
  float lambda = n_2 > 1e-2 ? (dot4(n, P) + plane.w) / n_2 : 0.0f;
  return plus4(P, mul3(-lambda, n));
}
template <typename T>
inline __device__ void swap_cuda(T& a, T& b) {
  T c(a);
  a = b;
  b = c;
}

inline __device__ float4 tri2plane(float3* vertices) {
  float3 normal = cross3(minus3(vertices[1], vertices[0]),
                         minus3(vertices[2], vertices[0]));
  // plane ax+by+cz+d=0 NO need to use unit normal
  // a normal vector is enough
  return plane_from_point_and_normal(vertices[0], normal);
}

inline __device__ float max4(float a, float b, float c, float d) {
  return fmaxf(fmaxf(a, b), fmaxf(c, d));
}

inline __device__ void get_minmax3(float& m, float& M, float x1, float x2,
                                   float x3) {
  m = fminf(fminf(x1, x2), x3);
  M = fmaxf(fmaxf(x1, x2), x3);
}

inline __device__ double max4(double a, double b, double c, double d) {
  return fmax(fmax(a, b), fmax(c, d));
}

inline __device__ void get_minmax3(double& m, double& M, double x1, double x2,
                                   double x3) {
  m = fmin(fmin(x1, x2), x3);
  M = fmax(fmax(x1, x2), x3);
}
