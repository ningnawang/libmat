#include <math.h>

#include <algorithm>
#include <array>

#include "convex_cell.h"
#include "kNN-CUDA/knncuda.h"
#include "stopwatch.h"
#include "voronoi.h"

// ###################  Status   ######################
char StatusStr[7][128] = {"triangle_overflow",
                          "vertex_overflow",
                          "inconsistent_boundary",
                          "security_radius_not_reached",
                          "success",
                          "needs_exact_predicates",
                          "no_intersection"};

void show_status_stats(std::vector<Status>& stat) {
  IF_VERBOSE(
      std::cerr << " \n\n\n---------Summary of success/failure------------\n");
  std::vector<int> nb_statuss(7, 0);
  FOR(i, stat.size()) nb_statuss[stat[i]]++;
  IF_VERBOSE(FOR(r, 7) std::cerr << " " << StatusStr[r] << "   "
                                 << nb_statuss[r] << "\n";)
  std::cerr << " " << StatusStr[4] << "   " << nb_statuss[4] << " /  "
            << stat.size() << "\n";
}

void cuda_check_error() {
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "Failed (1) (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

// ###################  Functions   ######################
__global__ void compute_new_site(float* site, int n_site, size_t site_pitch,
                                 float* cell_bary_sum,
                                 const size_t cell_bary_sum_pitch,
                                 float* cell_vol) {
  int seed = blockIdx.x * blockDim.x + threadIdx.x;
  if (seed >= n_site) return;

  if (cell_vol[seed] != 0) {
    if (site_pitch) {
      site[seed] = cell_bary_sum[seed] / cell_vol[seed];
      site[seed + site_pitch] =
          cell_bary_sum[seed + cell_bary_sum_pitch] / cell_vol[seed];
      site[seed + (site_pitch << 1)] =
          cell_bary_sum[seed + (cell_bary_sum_pitch << 1)] / cell_vol[seed];
    } else {
      site[3 * seed] = cell_bary_sum[3 * seed] / cell_vol[seed];
      site[3 * seed + 1] = cell_bary_sum[3 * seed + 1] / cell_vol[seed];
      site[3 * seed + 2] = cell_bary_sum[3 * seed + 2] / cell_vol[seed];
    }

    // if (seed == 588) {
    // printf(
    //     "seed %d has cell_bary_sum: (%f,%f,%f), cell_vol: %f, site: "
    //     "(%f,%f,%f) \n",
    //     seed, cell_bary_sum[seed], cell_bary_sum[seed + cell_bary_sum_pitch],
    //     cell_bary_sum[seed + (cell_bary_sum_pitch << 1)], cell_vol[seed],
    //     site[seed], site[seed + site_pitch], site[seed + (site_pitch << 1)]);
    // }
  }
}

__global__ void compute_tet_centroid(const float* vert, const int n_vert,
                                     const size_t vert_pitch, const int* idx,
                                     const int n_tet, const size_t idx_pitch,
                                     float* tet_centroid,
                                     const size_t tet_centroid_pitch) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n_tet) return;

  float3 centroid = {0.0f, 0.0f, 0.0f};

  FOR(i, 4) {
    centroid.x += vert[idx[tid + i * idx_pitch]];
    centroid.y += vert[idx[tid + i * idx_pitch] + vert_pitch];
    centroid.z += vert[idx[tid + i * idx_pitch] + (vert_pitch << 1)];
  }

  centroid.x *= 0.25f;
  centroid.y *= 0.25f;
  centroid.z *= 0.25f;

  tet_centroid[tid] = centroid.x;
  tet_centroid[tid + tet_centroid_pitch] = centroid.y;
  tet_centroid[tid + (tet_centroid_pitch << 1)] = centroid.z;
}

__global__ void transpose_site(const float* site, const int n_site,
                               float* site_transposed,
                               const size_t site_transposed_pitch) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  site_transposed[idx] = site[3 * idx];
  site_transposed[idx + site_transposed_pitch] = site[3 * idx + 1];
  site_transposed[idx + (site_transposed_pitch << 1)] = site[3 * idx + 2];
}

// ninwang:
// for each tet centroid, find n_tet number of nearest weighte sites
void compute_tet_weighted_knn_dev(const float* vert_dev, const int n_vert,
                                  const size_t vert_pitch, const int* idx_dev,
                                  const int n_tet, const size_t idx_pitch,
                                  const float* site_dev, const int n_site,
                                  const size_t site_pitch,
                                  const float* site_weights_dev,
                                  int* tet_knn_dev, const size_t tet_knn_pitch,
                                  const int tet_k) {
  float* tet_centroid_dev = nullptr;
  size_t tet_centroid_pitch_in_bytes;
  cudaMallocPitch((void**)&tet_centroid_dev, &tet_centroid_pitch_in_bytes,
                  n_tet * sizeof(float), 3);
  cuda_check_error();
  size_t tet_centroid_pitch = tet_centroid_pitch_in_bytes / sizeof(float);

  compute_tet_centroid<<<n_tet / VORO_BLOCK_SIZE + 1, VORO_BLOCK_SIZE>>>(
      vert_dev, n_vert, vert_pitch, idx_dev, n_tet, idx_pitch, tet_centroid_dev,
      tet_centroid_pitch);

  if (site_pitch)  // has been transposed
    knn_weighted_cuda_global_dev(site_dev, n_site, site_pitch, site_weights_dev,
                                 tet_centroid_dev, n_tet, tet_centroid_pitch, 3,
                                 tet_k, tet_knn_dev, tet_knn_pitch);
  else {
    float* site_transposed_dev = nullptr;
    size_t site_transposed_pitch_in_bytes;
    cudaMallocPitch((void**)&site_transposed_dev,
                    &site_transposed_pitch_in_bytes, n_site * sizeof(float), 3);
    cuda_check_error();
    size_t site_transposed_pitch =
        site_transposed_pitch_in_bytes / sizeof(float);
    transpose_site<<<n_site / KNN_BLOCK_SIZE + 1, KNN_BLOCK_SIZE>>>(
        site_dev, n_site, site_transposed_dev, site_transposed_pitch);

    knn_weighted_cuda_global_dev(site_transposed_dev, n_site,
                                 site_transposed_pitch, site_weights_dev,
                                 tet_centroid_dev, n_tet, tet_centroid_pitch, 3,
                                 tet_k, tet_knn_dev, tet_knn_pitch);

    cudaFree(site_transposed_dev);
  }

  cudaFree(tet_centroid_dev);
}

// return tet_sphere_relate_dev
// size #spheres x #tets
__global__ void tet_sphere_relations_dev(
    const int n_vert, const int* idx_dev, const int n_tet,
    const size_t idx_pitch, const int n_site, const uint* site_flags_dev,
    const int* site_knn_dev, const size_t site_knn_pitch, const int site_k,
    const float* pdist_dev, const size_t pdist_pitch,
    int* tet_sphere_relate_dev, size_t tet_sphere_relate_pitch) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n_tet) return;

  // each thread handles all n_site
  FOR(si_idx, n_site) {
    // case 1: if si is not selected, then tet is not related to si
    if (site_flags_dev[si_idx] != SiteFlag::is_selected) {
      tet_sphere_relate_dev[tid + si_idx * tet_sphere_relate_pitch] = 0;
      continue;
    }
    // case 2: find all neighboring site of sphere i (see note)
    int num_relate_hp = 0;
    int site_real_k = 0;  // some may be -1
    FOR(sm, site_k) {
      int sm_idx = site_knn_dev[si_idx + sm * site_knn_pitch];
      if (sm_idx == -1) continue;
      site_real_k += 1;
      FOR(t, 4) {
        int v_idx = idx_dev[tid + t * idx_pitch];
        // get power distances to sphere i and sphere m
        float pd_i = pdist_dev[v_idx + si_idx * pdist_pitch];
        float pd_m = pdist_dev[v_idx + sm_idx * pdist_pitch];
        // tet vertex v_idx closer to sphere i than j
        if (pd_m > pd_i) {
          num_relate_hp += 1;
          break;
        }
      }  // 4 tet vertices
    }  // for site_k (all neighbor spheres)

    int is_relate = num_relate_hp == site_real_k ? 1 : 0;
    tet_sphere_relate_dev[tid + si_idx * tet_sphere_relate_pitch] = is_relate;
  }  // for n_site
}

// ninwang
// update tet_knn and tet_k
// use power distance to find the potential relations between tet and sphere
void compute_tet_sphere_relation(
    const float* vert_dev, const int n_vert, const size_t vert_pitch,
    const int* idx_dev, const int n_tet, const size_t idx_pitch,
    const float* site_dev, const uint* site_flags_dev, const int n_site,
    const size_t site_pitch, const float* site_weights_dev,
    const int* site_knn_dev, const int site_k, const size_t site_knn_pitch,
    std::vector<int>& tet_knn, int& tet_k) {
  // printf("calling compute_tet_sphere_relation... \n");

  assert(site_pitch > 0);  // always transposed
  cudaError_t err0;
  // step 1: compute power distancese for each tet vertex
  // allocate global memory for power distance matrix
  // size: #spheres x #tet_vertices
  float* tet_pdist_dev = nullptr;
  size_t tet_pdist_pitch_in_bytes;
  err0 = cudaMallocPitch((void**)&tet_pdist_dev, &tet_pdist_pitch_in_bytes,
                         n_vert * sizeof(float), n_site);
  if (err0 != cudaSuccess) {
    printf("ERROR: Memory allocation error\n");
    cudaFree(tet_pdist_dev);
  }
  size_t tet_pdist_pitch = tet_pdist_pitch_in_bytes / sizeof(float);
  power_dist_cuda_global_dev(site_dev, n_site, site_pitch, site_weights_dev,
                             vert_dev, n_vert, vert_pitch, 3, tet_pdist_dev,
                             tet_pdist_pitch);

  // ninwang: debug
  // copy knn back to host
  cudaStreamSynchronize(0);
  std::vector<float> tet_pdist(n_site * n_vert);
  cudaMemcpy2D(tet_pdist.data(), n_vert * sizeof(float), tet_pdist_dev,
               tet_pdist_pitch_in_bytes, n_vert * sizeof(float), n_site,
               cudaMemcpyDeviceToHost);
  cuda_check_error();

  // printf("tet_pdist matrix: \n\t");
  // for (uint vid = 0; vid < n_vert; vid++) {    // column
  //   for (uint sid = 0; sid < n_site; sid++) {  // row
  //     printf("%f, ", tet_pdist[vid + sid * n_vert]);
  //   }
  //   printf("\n\t ");
  // }
  // printf("\n");
  // printf("done power_dist_cuda_global_dev \n");

  // step 2: compute tet-sphere relation matrix
  // size: #spheres x #tets
  // value:
  // 1: tet relates to sphere
  // 0: not relate
  int* tet_sphere_relate_dev = nullptr;
  size_t tet_sphere_relate_pitch_in_bytes;
  err0 = cudaMallocPitch((void**)&tet_sphere_relate_dev,
                         &tet_sphere_relate_pitch_in_bytes, n_tet * sizeof(int),
                         n_site);
  if (err0 != cudaSuccess) {
    printf("ERROR: Memory allocation error\n");
    cudaFree(tet_sphere_relate_dev);
  }
  size_t tet_sphere_relate_pitch =
      tet_sphere_relate_pitch_in_bytes / sizeof(int);
  tet_sphere_relations_dev<<<n_tet / VORO_BLOCK_SIZE + 1, VORO_BLOCK_SIZE>>>(
      n_vert, idx_dev, n_tet, idx_pitch, n_site, site_flags_dev, site_knn_dev,
      site_knn_pitch, site_k, tet_pdist_dev, tet_pdist_pitch,
      tet_sphere_relate_dev, tet_sphere_relate_pitch);

  // copy tet_sphere_relate_dev back to CPU
  cudaStreamSynchronize(0);
  std::vector<int> tet_sphere(n_site * n_tet);
  cudaMemcpy2D(tet_sphere.data(), n_tet * sizeof(int), tet_sphere_relate_dev,
               tet_sphere_relate_pitch_in_bytes, n_tet * sizeof(int), n_site,
               cudaMemcpyDeviceToHost);
  cuda_check_error();

  // printf("tet_sphere matrix: \n\t");
  // for (uint tid = 0; tid < n_tet; tid++) {  // column
  //   // for (uint sid = 0; sid < n_site; sid++) {  // row
  //   for (uint sid = 16; sid < 17; sid++) {  // row
  //     printf("%d ", tet_sphere[tid + sid * n_tet]);
  //   }
  //   printf("\n\t ");
  // }
  // printf("\n");

  // update tet_k = maximum number of related sphere per tet
  tet_k = 0;
  std::vector<std::vector<int>> tet_spheres_all(n_tet);
  for (uint tid = 0; tid < n_tet; tid++) {
    for (uint sid = 0; sid < n_site; sid++) {
      if (tet_sphere[tid + sid * n_tet] == 1) {
        tet_spheres_all[tid].push_back(sid);  // store related spheres
      }
    }
    if (tet_spheres_all[tid].size() > tet_k) {
      tet_k = tet_spheres_all[tid].size();

      // printf("tet_id %d has %d neighbors: [", tid,
      // tet_spheres_all[tid].size()); for (const auto& n :
      // tet_spheres_all[tid]) {
      //   printf("%d, ", n);
      // }
      // printf("]\n");
    }
  }
  printf("updated tet_k: %d\n", tet_k);

  // convert to flat 2D matrix
  // 1. size: (tet_k) x n_tet
  // 2. each column j store all related sphere of tet j
  // 3. init as -1
  tet_knn.clear();
  tet_knn.resize(tet_k * n_tet, -1);
  for (int tid = 0; tid < n_tet; tid++) {
    const auto& related_spheres = tet_spheres_all.at(tid);
    for (int sid = 0; sid < related_spheres.size(); sid++) {
      assert(sid <= tet_k);
      tet_knn[tid + sid * n_tet] = related_spheres[sid];
    }
  }

  // Memory clean-up
  cudaFree(tet_pdist_dev);
  cudaFree(tet_sphere_relate_dev);
}

void copy_tet_data(const std::vector<float>& vertices,
                   const std::vector<int>& indices, float*& vert_dev,
                   size_t& vert_pitch, int*& idx_dev, size_t& idx_pitch) {
  size_t n_vert = vertices.size() / 3, n_tet = (indices.size() >> 2);

  // transpose
  float* vert_T = new float[3 * n_vert];
  int* idx_T = new int[n_tet << 2];
  FOR(i, n_vert) {
    vert_T[i] = vertices[3 * i];
    vert_T[i + n_vert] = vertices[3 * i + 1];
    vert_T[i + (n_vert << 1)] = vertices[3 * i + 2];
  }
  FOR(i, n_tet) {
    idx_T[i] = indices[(i << 2)];
    idx_T[i + n_tet] = indices[(i << 2) + 1];
    idx_T[i + (n_tet << 1)] = indices[(i << 2) + 2];
    idx_T[i + n_tet * 3] = indices[(i << 2) + 3];
  }

  size_t vert_pitch_in_bytes, idx_pitch_in_bytes;
  cudaMallocPitch((void**)&vert_dev, &vert_pitch_in_bytes,
                  n_vert * sizeof(float), 3);
  cuda_check_error();
  cudaMallocPitch((void**)&idx_dev, &idx_pitch_in_bytes, n_tet * sizeof(int),
                  4);
  cuda_check_error();
  vert_pitch = vert_pitch_in_bytes / sizeof(float);
  idx_pitch = idx_pitch_in_bytes / sizeof(int);
  cudaMemcpy2D(vert_dev, vert_pitch_in_bytes, vert_T, n_vert * sizeof(float),
               n_vert * sizeof(float), 3, cudaMemcpyHostToDevice);
  cuda_check_error();
  cudaMemcpy2D(idx_dev, idx_pitch_in_bytes, idx_T, n_tet * sizeof(int),
               n_tet * sizeof(int), 4, cudaMemcpyHostToDevice);
  cuda_check_error();

  delete[] vert_T;
  delete[] idx_T;
}

/**
 * See load_tet_adj_info()
 *
 * 1. We store the number of shared adjacent cells for furture calculation of
 * Euler.
 * 2. We also stores unique id (int2) for each face, defined by (fid, -1)
 *
 * Vertex Eulers are stored as indices
 * Euler(v) = 1. / (#adjacent cells)
 *
 * Edge Eulers are stored as a diagonal adjacency matrix
 *
 * Face Eulers (originally from tets) are always 1/2 (adjacent to 2 cells)
 * except boundary
 */
void load_num_adjacent_cells_and_ids(const std::vector<int>& v_adjs,
                                     const std::vector<int>& e_adjs,
                                     const std::vector<int>& f_adjs,
                                     const std::vector<int>& f_ids,
                                     int*& v_adjs_dev, int*& e_adjs_dev,
                                     int*& f_adjs_dev, int*& f_ids_dev) {
  assert(!v_adjs.empty() && !e_adjs.empty() && !f_adjs.empty() &&
         !f_ids.empty());
  std::cout << "loaded #adjacent cells for v: " << v_adjs.size()
            << ", and e: " << e_adjs.size() << ", and f " << f_adjs.size()
            << std::endl;

  // ninwang: cuda need a pointer to pointer
  cudaMalloc((void**)&v_adjs_dev, v_adjs.size() * sizeof(int));
  cuda_check_error();
  cudaMemcpy(v_adjs_dev, v_adjs.data(), v_adjs.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cuda_check_error();

  cudaMalloc((void**)&e_adjs_dev, e_adjs.size() * sizeof(int));
  cuda_check_error();
  cudaMemcpy(e_adjs_dev, e_adjs.data(), e_adjs.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cuda_check_error();

  cudaMalloc((void**)&f_adjs_dev, f_adjs.size() * sizeof(int));
  cuda_check_error();
  cudaMemcpy(f_adjs_dev, f_adjs.data(), f_adjs.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cuda_check_error();

  cudaMalloc((void**)&f_ids_dev, f_ids.size() * sizeof(int));
  cuda_check_error();
  cudaMemcpy(f_ids_dev, f_ids.data(), f_ids.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cuda_check_error();
}

__host__ __device__ cuchar3 convert(uchar3 orig) {
  return cmake_uchar3(orig.x, orig.y, orig.z);
}
__host__ __device__ cuchar4 convert(uchar4 orig) {
  return cmake_uchar4(orig.x, orig.y, orig.z, orig.w);
}
__host__ __device__ cint2 convert(int2 orig) {
  return cmake_int2(orig.x, orig.y);
}
__host__ __device__ cfloat4 convert(float4 orig) {
  return cmake_float4(orig.x, orig.y, orig.z, orig.w);
}
__host__ __device__ cfloat5 convert(float5 orig) {
  return cmake_float5(orig.x, orig.y, orig.z, orig.w, orig.h);
}

void copy_cc(const ConvexCellTransfer& cc, ConvexCellHost& cc_trans) {
  cc_trans.is_active = true;
  cc_trans.status = cc.status;
  cc_trans.thread_id = cc.thread_id;
  cc_trans.voro_id = cc.voro_id;
  cc_trans.tet_id = cc.tet_id;
  cc_trans.euler = cc.euler;
  cc_trans.weight = cc.weight;
  cc_trans.nb_v = cc.nb_v;
  cc_trans.nb_p = cc.nb_p;
  cc_trans.nb_e = cc.nb_e;
  FOR(i, cc.nb_v) cc_trans.ver_data_trans[i] = convert(cc.ver_data_trans[i]);
  FOR(i, cc.nb_p) cc_trans.clip_data_trans[i] = convert(cc.clip_data_trans[i]);
  FOR(i, cc.nb_p)
  cc_trans.clip_id2_data_trans[i] = convert(cc.clip_id2_data_trans[i]);
  FOR(i, cc.nb_e) cc_trans.edge_data[i] = convert(cc.edge_data[i]);
}

/**
 * vertices: tet vertices
 * indices: tet 4 indices of vertices [can be parital indices]
 */
std::vector<ConvexCellHost> compute_clipped_voro_diagram_GPU(
    const int num_itr_global, const std::vector<float>& vertices,
    const std::vector<int>& indices, const std::map<int, std::set<int>>& v2tets,
    const std::vector<int>& v_adjs, const std::vector<int>& e_adjs,
    const std::vector<int>& f_adjs, const std::vector<int>& f_ids,
    std::vector<float>& site, const int n_site,
    const std::vector<float>& site_weights, const std::vector<uint>& site_flags,
    const std::vector<int>& site_knn, const int site_k,
    std::vector<float>& site_cell_vol, const bool site_is_transposed,
    int nb_Lloyd_iter, int preferred_tet_k) {
  cudaSetDevice(0);  // specify a device to be used for GPU computation
  int n_vert = vertices.size() / 3;
  int n_tet = (indices.size() >> 2);

  // copy tet vertices and indices to device
  float* vert_dev = nullptr;
  int* idx_dev = nullptr;
  size_t vert_pitch, idx_pitch;
  copy_tet_data(vertices, indices, vert_dev, vert_pitch, idx_dev, idx_pitch);

  // load adjacencies for Euler
  int* v_adjs_dev = nullptr;
  int* e_adjs_dev = nullptr;
  int* f_adjs_dev = nullptr;
  int* f_ids_dev = nullptr;
  load_num_adjacent_cells_and_ids(v_adjs, e_adjs, f_adjs, f_ids, v_adjs_dev,
                                  e_adjs_dev, f_adjs_dev, f_ids_dev);
  assert(v_adjs.size() == n_vert);

  // allocate memory for voronoi cell
  VoronoiCell* voronoi_cells_dev = nullptr;
  cudaMalloc((void**)&voronoi_cells_dev, n_site * sizeof(VoronoiCell));
  cuda_check_error();

  // allocate memory forninwang:   output points
  float* cell_bary_sum_dev = nullptr;
  size_t cell_bary_sum_pitch_in_bytes = 0, cell_bary_sum_pitch = 0;
  if (site_is_transposed) {
    cudaMallocPitch((void**)&cell_bary_sum_dev, &cell_bary_sum_pitch_in_bytes,
                    n_site * sizeof(float), 3);
    cell_bary_sum_pitch = cell_bary_sum_pitch_in_bytes / sizeof(float);
  } else
    cudaMalloc((void**)&cell_bary_sum_dev, 3 * n_site * sizeof(float));
  cuda_check_error();

  // allocate memory for cell volume
  site_cell_vol.clear();
  site_cell_vol.resize(n_site);
  float* cell_vol_dev = nullptr;
  cudaMalloc((void**)&cell_vol_dev, n_site * sizeof(float));
  cuda_check_error();

  // ninwang: allocate memory for site weights
  assert(site_weights.size() == n_site);
  float* site_weights_dev = nullptr;
  cudaMalloc((void**)&site_weights_dev, n_site * sizeof(float));
  cuda_check_error();
  cudaMemcpy(site_weights_dev, site_weights.data(), n_site * sizeof(float),
             cudaMemcpyHostToDevice);
  cuda_check_error();

  // ninwang: allocate memory for site flag
  assert(site_flags.size() == n_site);
  uint* site_flags_dev = nullptr;
  cudaMalloc((void**)&site_flags_dev, n_site * sizeof(uint));
  cuda_check_error();
  cudaMemcpy(site_flags_dev, site_flags.data(), n_site * sizeof(uint),
             cudaMemcpyHostToDevice);
  cuda_check_error();

  // allocate memory for site and site knn
  assert(site_knn.size() == (site_k + 1) * n_site);
  float* site_transposed_dev = nullptr;
  int* site_knn_dev = nullptr;
  size_t site_pitch = 0, site_pitch_in_bytes = 0, site_knn_pitch = 0,
         site_knn_pitch_in_bytes = 0;

  //////////////////////////////
  // Load Site and Site Neighbors
  {
    // allocate memory for site and site knn
    cudaMallocPitch((void**)&site_transposed_dev, &site_pitch_in_bytes,
                    n_site * sizeof(float), 3);
    cudaMallocPitch((void**)&site_knn_dev, &site_knn_pitch_in_bytes,
                    n_site * sizeof(int), site_k + 1);
    cuda_check_error();

    site_pitch = site_pitch_in_bytes / sizeof(float);
    site_knn_pitch = site_knn_pitch_in_bytes / sizeof(int);
    // printf("------------------site_pitch: %d, n_site: %d\n", site_pitch,
    //        n_site);
    assert(site_pitch != 0);

    // copy sites to device
    cudaMemcpy2D(site_transposed_dev, site_pitch_in_bytes, site.data(),
                 n_site * sizeof(float), n_site * sizeof(float), 3,
                 cudaMemcpyHostToDevice);
    cuda_check_error();

    printf("done copy site to device \n");

    // copy site_knn to device, init site_knn as -1
    // site_knn is 2d flat matrix
    // dim: (site_k+1) x n_site
    // each column j store all neighbors of sphere all_medial_spheres.at(j)
    cudaMemcpy2D(site_knn_dev, site_knn_pitch_in_bytes, site_knn.data(),
                 n_site * sizeof(int), n_site * sizeof(int), site_k + 1,
                 cudaMemcpyHostToDevice);
    cuda_check_error();

    printf("done copy site_knn to device \n");

  }  // Site and Site Neighbors

  //////////////////////////////
  // Store records
  std::ofstream record("record.csv", std::ios::app);
  record << "n_site, n_tet, site_k, tet_k, Tet_Sphere, Compute_RPD, "
            "GPU2CPU, Non_Dup_RPCs\n";
  record << std::setprecision(5) << std::setiosflags(std::ios::fixed);

  Stopwatch sw("Iteration");
  double start_time = 0.0, stop_time = 0.0;
  start_time = sw.now();
  record << n_site << ", " << n_tet << ", " << site_k << ", ";

  //////////////////////////////
  // Tet-Sphere
  cudaStreamSynchronize(0);  // for printf in device code
  int tet_k = -1;
  std::vector<int> tet_knn;  // init as -1
  int* tet_knn_dev = nullptr;
  size_t tet_knn_pitch, tet_knn_pitch_in_bytes;
  {
    // update tet_knn and tet_k
    compute_tet_sphere_relation(
        vert_dev, n_vert, vert_pitch, idx_dev, n_tet, idx_pitch,
        site_transposed_dev, site_flags_dev, n_site, site_pitch,
        site_weights_dev, site_knn_dev, site_k, site_knn_pitch, tet_knn, tet_k);

    // allocate memory for tet knn (2D array)
    cudaMallocPitch((void**)&tet_knn_dev, &tet_knn_pitch_in_bytes,
                    n_tet * sizeof(int), tet_k);
    cuda_check_error();
    tet_knn_pitch = tet_knn_pitch_in_bytes / sizeof(int);
    // copy tet_knn to device
    cudaMemcpy2D(tet_knn_dev, tet_knn_pitch_in_bytes, tet_knn.data(),
                 n_tet * sizeof(int), n_tet * sizeof(int), tet_k,
                 cudaMemcpyHostToDevice);
    cuda_check_error();
  }
  // // ninwang: debug
  // // copy knn back to host
  // cudaStreamSynchronize(0);
  // tet_knn.resize(n_tet * tet_k);
  // cudaMemcpy2D(tet_knn.data(), n_tet * sizeof(int), tet_knn_dev,
  //              tet_knn_pitch_in_bytes, n_tet * sizeof(int), tet_k,
  //              cudaMemcpyDeviceToHost);
  // cuda_check_error();

  // printf("tet_knn matrix after: \n\t");
  // // for (uint tid = 0; tid < n_tet; tid++) {  // column
  // // for (uint tid = 145; tid < 146; tid++) {  // column
  // uint tid = 33;
  // for (uint i = 0; i < tet_k; i++) {  // row
  //   printf("%d ", tet_knn[i * n_tet + tid]);
  // }
  // printf("\n\t ");
  // // }
  // printf("\n");
  // printf("done compute_tet_weighted_knn_dev \n");

  stop_time = sw.now();  // record Tet_Sphere
  record << tet_k << ", " << stop_time - start_time << ", ";
  start_time = sw.now();

  //////////////////////////////
  // Barycenters
  {
    if (site_is_transposed)
      cudaMemset2D(cell_bary_sum_dev, cell_bary_sum_pitch_in_bytes, 0,
                   n_site * sizeof(float), 3);
    else
      cudaMemset(cell_bary_sum_dev, 0, 3 * n_site * sizeof(float));
    cudaMemset(cell_vol_dev, 0, n_site * sizeof(float));
    cuda_check_error();
  }

  ////////////////////////////////////////////
  // Compute RPD
  //
  // GPU: total thread size = n_grids * n_blocks
  // Each thread compute a ConvexCell defined by
  // one tet and one nearby seed.
  //
  // Note: will be updated later by tet_k
  int n_grids = n_tet * tet_k / VORO_BLOCK_SIZE + 1;
  int n_blocks = VORO_BLOCK_SIZE;
  printf("n_vert: %d, n_tet: %d, n_site: %d, tet_k: %d, site_k: %d\n", n_vert,
         n_tet, n_site, tet_k, site_k);
  printf("n_grids: %d, n_blocks: %d, #threads: %d \n", n_grids, n_blocks,
         n_grids * n_blocks);

  // allocate more, after tet_k been updated
  // by function vcompute_tet_sphere_relation()
  // ninwang: allocate memory for all convex cells
  std::vector<ConvexCellTransfer> convex_cells_host(n_grids * n_blocks);
  ConvexCellTransfer* convex_cells_dev = nullptr;
  cudaMalloc((void**)&convex_cells_dev,
             n_grids * n_blocks * sizeof(ConvexCellTransfer));
  cuda_check_error();

  // allocate memory for stats
  std::vector<Status> stat(n_tet * tet_k, security_radius_not_reached);
  GPUBuffer<Status> gpu_stat(stat);

  {  // GPU voro kernel only
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    cuda_check_error();

    printf("calling clipped_voro_cell_test_GPU_param_tet \n");
    printf("n_grids: %d, n_blocks: %d, #threads: %d \n", n_grids, n_blocks,
           n_grids * n_blocks);

    // clip tet-cell (init as tet)
    // "voronoi_cells_dev" is not used
    cudaStreamSynchronize(0);  // for printf in device code
    clipped_voro_cell_test_GPU_param_tet<<<n_grids, n_blocks>>>(
        site_transposed_dev, n_site, site_pitch, site_weights_dev,
        site_flags_dev, site_knn_dev, site_knn_pitch, site_k, vert_dev, n_vert,
        vert_pitch, idx_dev, n_tet, idx_pitch, v_adjs_dev, e_adjs_dev,
        f_adjs_dev, f_ids_dev, tet_knn_dev, tet_knn_pitch, tet_k,
        gpu_stat.gpu_data, voronoi_cells_dev, convex_cells_dev,
        cell_bary_sum_dev, cell_bary_sum_pitch, cell_vol_dev);
    cuda_check_error();

    cudaStreamSynchronize(0);  // for printf in device code
    cudaEventRecord(stop);
    cudaEventSynchronize(start);
    cudaEventSynchronize(stop);
    // printf("done clipped_voro_cell_test_GPU_param_tet \n");

    stop_time = sw.now();  // record Compute_RPD
    record << stop_time - start_time << ", ";
    start_time = sw.now();

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }  // GPU voro kernel only

  ////////////////////////////////////////////
  // GPU2CPU: copy data back to the cpu
  {
    start_time = sw.now();
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    cudaMemcpy2D(site.data(), n_site * sizeof(float), site_transposed_dev,
                 site_pitch_in_bytes, n_site * sizeof(float), 3,
                 cudaMemcpyDeviceToHost);

    // ninwang: copy all concave cells
    cudaMemcpy(convex_cells_host.data(), convex_cells_dev,
               n_grids * n_blocks * sizeof(ConvexCellTransfer),
               cudaMemcpyDeviceToHost);

    cudaEventRecord(stop);
    cudaEventSynchronize(start);
    cudaEventSynchronize(stop);

    stop_time = sw.now();  // gpu2cpu
    record << stop_time - start_time << ", ";

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }  // copy data back to the cpu

  ////////////////////////////////////////////
  // GPU2CPU: copy data back to the cpu
  //
  // stores non-duplicated convex_cells_host
  // seed -> tet_ids, do not process duplicates
  printf("collecting convex_cells_host_non_dup ... \n");
  start_time = sw.now();
  std::vector<ConvexCellHost> convex_cells_host_non_dup;
  std::map<int, std::set<int>> tets2seed;  // mostly orderd by tet_ids
  for (auto& cc_trans : convex_cells_host) {
    // this is important!
    // to avoid random value assigned in Status::early_return
    if (!is_convex_cell_valid(cc_trans)) continue;
    // check if this seed&tet pair has been stored
    // each seed&tet pair should be unique but we might calculate
    // multiple times because of multi-thread, same idea used in
    // get_voro_cell_euler()
    auto& seed_set = tets2seed[cc_trans.tet_id];
    if (seed_set.find(cc_trans.voro_id) != seed_set.end()) continue;
    seed_set.insert(cc_trans.voro_id);
    // if (cc_trans.voro_id == 4 && cc_trans.tet_id == 33)
    //   printf("cc_trans.tet_id %d has cc_trans.voro_id: %d\n ",
    //   cc_trans.tet_id,
    //          cc_trans.voro_id);

    ConvexCellHost cc_new;
    copy_cc(cc_trans, cc_new);
    convex_cells_host_non_dup.push_back(cc_new);
    // easier for debug
    // assign id for each convex cell as ConvexCellHost::id
    // matching index in convex_cells_host_non_dup
    convex_cells_host_non_dup.back().id = convex_cells_host_non_dup.size() - 1;
  }
  stop_time = sw.now();  // Non_Dup_RPCs
  record << stop_time - start_time << ", ";
  printf("saved %zu/%zu unduplicated convex cells \n",
         convex_cells_host_non_dup.size(), convex_cells_host.size());

  record << std::endl;

  cudaFree(vert_dev);
  cudaFree(idx_dev);
  cudaFree(v_adjs_dev);
  cudaFree(e_adjs_dev);
  cudaFree(f_adjs_dev);
  cudaFree(f_ids_dev);
  cudaFree(tet_knn_dev);
  cudaFree(cell_vol_dev);
  cudaFree(voronoi_cells_dev);
  cudaFree(convex_cells_dev);
  cudaFree(cell_bary_sum_dev);
  cudaFree(site_transposed_dev);
  cudaFree(site_knn_dev);
  cudaFree(site_weights_dev);

  record.close();

  return convex_cells_host_non_dup;
}
