
#ifndef H_VORONOI_H
#define H_VORONOI_H

#include <assert.h>
#include <cuda_runtime.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common_cuda.h"
#include "voronoi_common.h"
#include "voronoi_defs.h"

//----------------------------------WRAPPER
template <class T>
struct GPUBuffer {
  void init(T* data) {
    IF_VERBOSE(std::cerr << "GPU: " << size * sizeof(T) / 1048576 << " Mb used"
                         << std::endl);
    cpu_data = data;
    cuda_check(cudaMalloc((void**)&gpu_data, size * sizeof(T)));
    cpu2gpu();
  }
  GPUBuffer(std::vector<T>& v) {
    size = v.size();
    init(v.data());
  }
  ~GPUBuffer() { cuda_check(cudaFree(gpu_data)); }

  void cpu2gpu() {
    cuda_check(cudaMemcpy(gpu_data, cpu_data, size * sizeof(T),
                          cudaMemcpyHostToDevice));
  }
  void gpu2cpu() {
    cuda_check(cudaMemcpy(cpu_data, gpu_data, size * sizeof(T),
                          cudaMemcpyDeviceToHost));
  }

  T* cpu_data;
  T* gpu_data;
  int size;
};

std::vector<ConvexCellHost> compute_clipped_voro_diagram_GPU(
    const int num_itr_global, const std::vector<float>& vertices,
    const std::vector<int>& indices, const std::map<int, std::set<int>>& v2tets,
    const std::vector<int>& v_adjs, const std::vector<int>& e_adjs,
    const std::vector<int>& f_adjs, const std::vector<int>& f_ids,
    std::vector<float>& site, const int n_site,
    const std::vector<float>& site_weights, const std::vector<uint>& site_flags,
    const std::vector<int>& site_knn, const int site_k,
    std::vector<float>& site_cell_vol, const bool site_is_transposed,
    int nb_Lloyd_iter = 1, int preferred_tet_k = 0);

#endif  // __VORONOI_H__