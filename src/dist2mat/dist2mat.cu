#include <cub/cub.cuh>

#include "dist2mat.h"

__host__ __device__ __forceinline__ float distance_to_sphere(const float3 p,
                                                             const float4 sp) {
  return length(p - make_float3(sp.x, sp.y, sp.z)) - sp.w;
}

__host__ __device__ __forceinline__ float4 bary_lerp(const float4 m1,
                                                     const float4 m2,
                                                     const float4 m3, float t1,
                                                     float t2) {
  return m1 * t1 + m2 * t2 + m3 * (1.f - t1 - t2);
}

__host__ __device__ __forceinline__ bool solve_quadratic(float A, float B,
                                                         float C, float* root) {
  root[0] = -1.f;
  root[1] = -1.f;
  if (A == 0.f) {  // degnerated case
    if (B != 0.f) {
      root[0] = -C / B;
      root[1] = root[0];
    }
  } else {
    float delta = B * B - 4.f * A * C;
    if (delta < 0.f) {
      return false;
    } else {
      root[0] = (-B - sqrtf(delta)) / (2 * A);
      root[1] = (-B + sqrtf(delta)) / (2 * A);
    }
  }
  return true;
}

__host__ __device__ __forceinline__ float compute_distance_to_cone(
    const float3 pos, const float4 m1, const float4 m2) {
  float3 c1 = make_float3(m1.x, m1.y, m1.z);
  float3 c2 = make_float3(m2.x, m2.y, m2.z);
  float r1 = m1.w, r2 = m2.w;

  // bool inversed = false;
  if (r1 > r2) {
    float3 tmp = c1;
    c1 = c2;
    c2 = tmp;
    float tmp_flt = r1;
    r1 = r2;
    r2 = tmp_flt;
  }

  float3 c21 = c1 - c2;
  float3 cq2 = c2 - pos;
  float A = dot(c21, c21);
  float D = 2.f * dot(c21, cq2);
  float F = dot(cq2, cq2);
  float R1 = r1 - r2;

  float t = -(A * D - R1 * R1 * D) -
            sqrtf((D * D - 4.0 * A * F) * (R1 * R1 - A) * R1 * R1);
  t /= 2.f * (A * A - A * R1 * R1);
  t = clamp(t, 0.f, 1.f);

  float3 lerp_sp = lerp(c1, c2, t);
  float lerp_r = lerp(r1, r2, t);
  return length(pos - lerp_sp) - lerp_r;
}

__host__ __device__ __forceinline__ float compute_distance_to_slab(
    const float3 pos, const float4 m1, const float4 m2, const float4 m3) {
  float3 c31 = make_float3(m1.x - m3.x, m1.y - m3.y, m1.z - m3.z);
  float3 c32 = make_float3(m2.x - m3.x, m2.y - m3.y, m2.z - m3.z);
  float3 cm3 = make_float3(m3.x - pos.x, m3.y - pos.y, m3.z - pos.z);

  float R1 = m1.w - m3.w;
  float R2 = m2.w - m3.w;
  float A = dot(c31, c31);
  float B = 2.f * dot(c31, c32);
  float C = dot(c32, c32);
  float D = 2.f * dot(c31, cm3);
  float E = 2.f * dot(c32, cm3);
  float F = dot(cm3, cm3);

  float t1 = -1.f, t2 = -1.f;
  if (R1 == 0.f && R2 == 0.f) {
    float denom = 4.f * A * C - B * B;
    t1 = (B * E - 2.0 * C * D) / denom;
    t2 = (B * D - 2.0 * A * E) / denom;
  } else if (R1 != 0.f && R2 == 0.f) {
    float H2 = -B / (2.f * C);
    float K2 = -E / (2.f * C);
    float W1 = powf(2.f * A + B * H2, 2.f) -
               4.f * R1 * R1 * (A + B * H2 + C * H2 * H2);
    float W2 = 2.f * (2.f * A + B * H2) * (B * K2 + D) -
               4.f * R1 * R1 * (B * K2 + 2.f * C * H2 * K2 + D + E * H2);
    float W3 =
        powf(B * K2 + D, 2.f) - 4.f * R1 * R1 * (C * K2 * K2 + E * K2 + F);
    float root[2];  // {t11, t12}
    solve_quadratic(W1, W2, W3, root);
    float t21 = H2 * root[0] + K2;
    float t22 = H2 * root[1] + K2;
    float dis = distance_to_sphere(pos, bary_lerp(m1, m2, m3, root[0], t21));
    t1 = root[0];
    t2 = t21;
    float dis2 = distance_to_sphere(pos, bary_lerp(m1, m2, m3, root[1], t22));
    if (dis2 < dis) {
      t1 = root[1];
      t2 = t22;
    }
  } else if (R1 == 0.f && R2 != 0.f) {
    float H1 = -B / (2.f * A);
    float K1 = -D / (2.f * A);
    float W1 = powf(2.f * C + B * H1, 2.f) -
               4.f * R2 * R2 * (C + B * H1 + A * H1 * H1);
    float W2 = 2.f * (2.f * C + B * H1) * (B * K1 + E) -
               4.f * R2 * R2 * (B * K1 + 2.f * A * H1 * K1 + E + D * H1);
    float W3 =
        powf(B * K1 + E, 2.f) - 4.f * R2 * R2 * (A * K1 * K1 + D * K1 + F);
    float root[2];  // {t21, t22}
    solve_quadratic(W1, W2, W3, root);
    float t11 = H1 * root[0] + K1;
    float t12 = H1 * root[1] + K1;
    float dis = distance_to_sphere(pos, bary_lerp(m1, m2, m3, t11, root[0]));
    t1 = t11;
    t2 = root[0];
    float dis2 = distance_to_sphere(pos, bary_lerp(m1, m2, m3, t12, root[1]));
    if (dis2 < dis) {
      t1 = t12;
      t2 = root[1];
    }
  } else {
    float L1 = 2.f * A * R2 - B * R1;
    float L2 = 2.f * C * R1 - B * R2;
    float L3 = E * R1 - D * R2;
    if (L1 == 0.f && L2 != 0.f) {
      t2 = -L3 / L2;
      float W1 = 4.f * A * A - 4.f * R1 * R1 * A;
      float W2 = 4.f * A * (B * t2 + D) - 4.f * R1 * R1 * (B * t2 + D);
      float W3 = powf(B * t2 + D, 2.f) - (C * t2 * t2 + E * t2 + F);
      float root[2];  // t11, t12
      solve_quadratic(W1, W2, W3, root);
      float dis = distance_to_sphere(pos, bary_lerp(m1, m2, m3, root[0], t2));
      t1 = root[0];
      if (distance_to_sphere(pos, bary_lerp(m1, m2, m3, root[1], t2)) < dis) {
        t1 = root[1];
      }
    } else if (L1 != 0.f && L2 == 0.f) {
      t1 = L3 / L1;
      float W1 = 4.f * C * C - 4.f * R2 * R2 * C;
      float W2 = 4.f * C * (B * t1 + E) - 4.f * R2 * R2 * (B * t1 + E);
      float W3 = powf(B * t1 + E, 2.f) - (A * t1 * t1 + D * t1 + F);
      float root[2];  // t21, t22
      solve_quadratic(W1, W2, W3, root);
      float dis = distance_to_sphere(pos, bary_lerp(m1, m2, m3, t1, root[0]));
      t2 = root[0];
      if (distance_to_sphere(pos, bary_lerp(m1, m2, m3, t1, root[1])) < dis) {
        t2 = root[1];
      }
    } else {
      float H3 = L2 / L1;
      float K3 = L3 / L1;
      float W1 = powf(2.f * C + B * H3, 2.f) -
                 4.f * R2 * R2 * (A * H3 * H3 + B * H3 + C);
      float W2 = 2.f * (2.f * C + B * H3) * (B * K3 + E) -
                 4.f * R2 * R2 * (2.f * A * H3 * K3 + B * K3 + D * H3 + E);
      float W3 =
          powf(B * K3 + E, 2.f) - 4.f * R2 * R2 * (A * K3 * K3 + D * K3 + F);
      float root[2];  // t21, t22
      solve_quadratic(W1, W2, W3, root);
      float t11 = H3 * root[0] + K3;
      float t12 = H3 * root[1] + K3;
      float dis = distance_to_sphere(pos, bary_lerp(m1, m2, m3, t11, root[0]));
      t1 = t11;
      t2 = root[0];
      if (distance_to_sphere(pos, bary_lerp(m1, m2, m3, t12, root[1])) < dis) {
        t1 = t12;
        t2 = root[1];
      }
    }
  }

  if ((t1 + t2) < 1.f && t1 >= 0.f && t1 <= 1.f && t2 >= 0.f && t2 <= 1.f) {
    return distance_to_sphere(pos,
                              bary_lerp(m1, m2, m3, t1, t2));  // valid solution
  }

  float dis1 = compute_distance_to_cone(pos, m1, m3);
  float dis2 = compute_distance_to_cone(pos, m2, m3);
  float dis3 = compute_distance_to_cone(pos, m1, m2);
  return fminf(dis1, fminf(dis2, dis3));
}

__global__ void ClosestDistanceToLocalMat(const float3* samples,
                                          const float4* spheres,
                                          const int3* prims,
                                          const uint* offsets,
                                          const uint* num_per_sample,
                                          float* result, int* closest_id) {
  const int tid = threadIdx.x;
  const int blockid = blockIdx.x;  // block id, each block process a sample
  // max = num_sample*kWarpSize(32)
  // may > num_prims when num_prim per sample < 32
  // [no use]
  // const int gid = blockIdx.x * kWarpSize + threadIdx.x;

  // block-level data
  int num_prim = num_per_sample[blockIdx.x];
  int offset = offsets[blockIdx.x];
  float3 pos = samples[blockIdx.x];
  float local_closest_dist = 1e16f;
  int local_closest_id = -1;

  // if (threadIdx.x == 0) {
  //   printf("tid %d num_prim %d, offset %d, pos (%f,%f,%f)\n", tid, num_prim,
  //          offset, pos.x, pos.y, pos.z);
  // }

  typedef cub::BlockReduce<float, kWarpSize> BlockReduce;
  __shared__ typename BlockReduce::TempStorage tmp_storage;
  __shared__ float min_distance[kWarpSize];
  __shared__ int min_id[kWarpSize];

  // put range check here instead of early return to make sure BlockReduce
  // works correctly
  for (int i = tid; i < num_prim && tid < num_prim;
       i += kWarpSize)  // block stride for-loop
  {
    // if (blockid == 2) printf("i: %d, gid %d, tid: %d\n", i, gid, tid);
    int3 prim = prims[offset + i];
    float dist = 1e16f;
    if (prim.x == -1 && prim.y == -1)  // to sphere
    {
      dist = distance_to_sphere(pos, spheres[prim.z]);
    } else if (prim.x == -1 && prim.y != -1)  // to cone
    {
      dist = compute_distance_to_cone(pos, spheres[prim.y], spheres[prim.z]);
    } else if (prim.x != -1)  // to slab
    {
      dist = compute_distance_to_slab(pos, spheres[prim.x], spheres[prim.y],
                                      spheres[prim.z]);
    } else  // to unknown
    {
      printf("Input data error: [%d, %d, %d]\n", prim.x, prim.y, prim.z);
    }

    if (dist < local_closest_dist) {
      local_closest_dist = fminf(local_closest_dist, dist);
      local_closest_id = i;
    }

    // uncomment to check if each block solves correct number of mat primitives
    // if (blockid == 2)
    //   printf("i: %d, gid %d, local closest dist: %f, local_closest_id: %d\n",
    //   i,
    //          gid, local_closest_dist, local_closest_id);
  }
  min_distance[tid] = local_closest_dist;
  min_id[tid] = local_closest_id;
  __syncthreads();

  const float reduced =
      BlockReduce(tmp_storage).Reduce(local_closest_dist, cub::Min());

  if (tid == 0) {
    // printf("gid: %d, blockid %d\n", gid, blockid);
    result[blockid] = reduced;  // write back to global buffer
    for (int i = 0; i < kWarpSize; ++i) {
      if (fabsf(min_distance[i] - reduced) < 1e-10f) {
        // if (blockid == 0)
        //   printf("final closest dist: %f, i %d, min_id %d\n", reduced, i,
        //          min_id[i]);
        closest_id[blockid] = min_id[i];
      }
    }
  }
}

void compute_closest_dist2mat(GpuBuffer<float4>& spheres, const int num_samples,
                              GpuBuffer<float3>& samples,
                              GpuBuffer<uint>& offset,
                              GpuBuffer<uint>& num_per_sample,
                              GpuBuffer<int3>& prims, GpuBuffer<float>& results,
                              GpuBuffer<int>& closest_mat_id) {
  // host to device
  spheres.H2D();
  samples.H2D();
  offset.H2D();
  num_per_sample.H2D();
  prims.H2D();
  results.H2D();
  closest_mat_id.H2D();

  // actual launch
  // int block_num = (num_samples + kWarpSize - 1) / kWarpSize;

  // kernel是按照一个block处理一个sample点，一个thread处理一个sample的一个MAT
  // primitive或者多个,
  // 一个sample对应的MAT如果小于32个的话，一个thread就处理一个，如果超过32个话，第i个thread会处理第i个，第i+32个这样子
  int block_num = num_samples;
  ClosestDistanceToLocalMat<<<block_num, kWarpSize>>>(
      samples.DPtr(), spheres.DPtr(), prims.DPtr(), offset.DPtr(),
      num_per_sample.DPtr(), results.DPtr(), closest_mat_id.DPtr());

  results.D2H();
  closest_mat_id.D2H();

  // for (int i = 0; i < num_samples; i++) {
  //   int clost_prim_id = closest_mat_id.HPtr()[i];
  //   float clost_dist = results.HPtr()[i];
  //   printf("sample %d has closest prim %d with dist %f\n", i, clost_prim_id,
  //          clost_dist);
  // }
}

// int main() {
//   const std::string sample_path = "../sample.samples";
//   std::ifstream stream;
//   stream.open(sample_path);

//   int num_spheres;
//   int num_samples;

//   int total_prim_num = 0;
//   GpuBuffer<float4> spheres;
//   GpuBuffer<float3> samples;
//   GpuBuffer<uint> offset;
//   GpuBuffer<uint> num_per_sample;
//   GpuBuffer<int3> prims;
//   GpuBuffer<int> closest_mat_id;
//   GpuBuffer<float> results;

//   std::vector<int3> prim_per_samples;  // temporary

//   if (stream.is_open()) {
//     stream >> num_spheres;
//     spheres.HResize(num_spheres);
//     for (int i = 0; i < num_spheres; ++i) {
//       float4 s;
//       stream >> s.x >> s.y >> s.z >> s.w;
//       spheres.HPtr()[i] = s;
//     }

//     stream >> num_samples;
//     samples.HResize(num_samples);
//     results.HResize(num_samples);
//     offset.HResize(num_samples);
//     closest_mat_id.HResize(num_samples);
//     num_per_sample.HResize(num_samples);
//     for (int i = 0; i < num_samples; ++i) {
//       float3 pos;
//       stream >> pos.x >> pos.y >> pos.z;
//       samples.HPtr()[i] = pos;
//       results.HPtr()[i] = 1e28f;
//       uint num_prims;
//       stream >> num_prims;
//       num_per_sample.HPtr()[i] = num_prims;
//       offset.HPtr()[i] = total_prim_num;
//       total_prim_num += num_prims;

//       for (int j = 0; j < num_prims; ++j) {
//         int3 prim;
//         stream >> prim.x >> prim.y >> prim.z;
//         prim_per_samples.push_back(prim);
//       }

//       closest_mat_id.HPtr()[i] = -1;
//     }
//   }
//   stream.close();

//   prims.HResize(prim_per_samples.size());
//   for (size_t i = 0; i < prim_per_samples.size(); ++i) {
//     prims.HPtr()[i] = prim_per_samples[i];
//   }

//   // host to device
//   spheres.H2D();
//   samples.H2D();
//   offset.H2D();
//   num_per_sample.H2D();
//   prims.H2D();
//   results.H2D();
//   closest_mat_id.H2D();

//   // actual launch
//   int block_num = (num_samples + kWarpSize - 1) / kWarpSize;
//   ClosestDistanceToLocalMat<<<block_num, kWarpSize>>>(
//       total_prim_num, samples.DPtr(), spheres.DPtr(), prims.DPtr(),
//       offset.DPtr(), num_per_sample.DPtr(), results.DPtr(),
//       closest_mat_id.DPtr());

//   /*
//   * unit test
//   float3 pos = make_float3(0.f, 0.7f, -0.35f);
//   float4 m2 = make_float4(0.5f, 0.f,0.f, 0.35f);
//   float4 m3 = make_float4(-0.5f, 0.0f, 0.0f, 0.5f);
//   float4 m4 = make_float4(0.f, 0.0f, -0.6f, 0.2f);

//   float dis = compute_distance_to_slab(pos, m2, m3, m4);
//   printf("distance: %f\n", dis);
//   */

//   return 0;
// }