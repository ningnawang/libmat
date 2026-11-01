// CUDA backend of libmat::distance_to_medial_mesh(), reusing the dist2mat
// kernel. All samples share a single primitive list (offset 0, #prims each).

#include <limits>

#include "api/internal/cuda_backend.h"
#include "libmat/cuda/dist2mat/dist2mat.h"

namespace libmat {
namespace internal {

void distance_to_medial_mesh_cuda(const Eigen::MatrixXd& P,
                                  const Eigen::MatrixXd& C,
                                  const Eigen::VectorXd& R,
                                  const Eigen::MatrixXi& E,
                                  const Eigen::MatrixXi& F, Eigen::VectorXd& D,
                                  Eigen::VectorXi& prim_type,
                                  Eigen::VectorXi& prim_id) {
  const int num_samples = P.rows();
  const int ns = C.rows();
  const int ne = E.size() > 0 ? E.rows() : 0;
  const int nf = F.size() > 0 ? F.rows() : 0;
  const int num_prims = ns + ne + nf;
  D.resize(num_samples);
  prim_type.resize(num_samples);
  prim_id.resize(num_samples);
  if (num_samples == 0) return;

  GpuBuffer<float4> spheres;
  GpuBuffer<float3> samples;
  GpuBuffer<uint> offset;
  GpuBuffer<uint> num_per_sample;
  GpuBuffer<int3> prims;
  GpuBuffer<float> results;
  GpuBuffer<int> closest_mat_id;

  spheres.HResize(ns);
  for (int i = 0; i < ns; ++i)
    spheres.HPtr()[i] = make_float4(C(i, 0), C(i, 1), C(i, 2), R(i));

  // primitive encoding of the dist2mat kernel:
  // sphere (-1,-1,s), cone (-1,s0,s1), slab (s0,s1,s2)
  prims.HResize(num_prims);
  for (int i = 0; i < ns; ++i) prims.HPtr()[i] = make_int3(-1, -1, i);
  for (int e = 0; e < ne; ++e)
    prims.HPtr()[ns + e] = make_int3(-1, E(e, 0), E(e, 1));
  for (int f = 0; f < nf; ++f)
    prims.HPtr()[ns + ne + f] = make_int3(F(f, 0), F(f, 1), F(f, 2));

  samples.HResize(num_samples);
  offset.HResize(num_samples);
  num_per_sample.HResize(num_samples);
  results.HResize(num_samples);
  closest_mat_id.HResize(num_samples);
  for (int i = 0; i < num_samples; ++i) {
    samples.HPtr()[i] = make_float3(P(i, 0), P(i, 1), P(i, 2));
    offset.HPtr()[i] = 0;
    num_per_sample.HPtr()[i] = num_prims;
    results.HPtr()[i] = 1e28f;
    closest_mat_id.HPtr()[i] = -1;
  }

  compute_closest_dist2mat(spheres, num_samples, samples, offset,
                           num_per_sample, prims, results, closest_mat_id);

  for (int i = 0; i < num_samples; ++i) {
    const int id = closest_mat_id.HPtr()[i];
    D(i) = results.HPtr()[i];
    if (id < 0) {
      prim_type(i) = -1, prim_id(i) = -1;
    } else if (id < ns) {
      prim_type(i) = 0, prim_id(i) = id;
    } else if (id < ns + ne) {
      prim_type(i) = 1, prim_id(i) = id - ns;
    } else {
      prim_type(i) = 2, prim_id(i) = id - ns - ne;
    }
  }
}

}  // namespace internal
}  // namespace libmat
