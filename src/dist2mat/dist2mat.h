
#ifndef H_DIST2MAT_H
#define H_DIST2MAT_H

// #include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#include "cuda_helper_math.h"
#include "cuda_utils.h"

#define kWarpSize 32

void compute_closest_dist2mat(GpuBuffer<float4>& spheres, const int num_samples,
                              GpuBuffer<float3>& samples,
                              GpuBuffer<uint>& offset,
                              GpuBuffer<uint>& num_per_sample,
                              GpuBuffer<int3>& prims, GpuBuffer<float>& results,
                              GpuBuffer<int>& closest_mat_id);

#endif  // __DIST2MAT_H__