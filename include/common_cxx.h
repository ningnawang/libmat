#pragma once

#include <limits.h>

#include <algorithm>
#include <array>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "assert.h"
#include "common.h"

//----------------------------------------------------------------------------

struct __attribute__((aligned(8))) cint2 {
  int x, y;
};
struct __attribute__((aligned(2))) cuchar2 {
  unsigned char x, y;
};
struct __attribute__((aligned(4))) cuchar3 {
  unsigned char x, y, z;
};
struct __attribute__((aligned(4))) cuchar4 {
  unsigned char x, y, z, w;
};
struct __attribute__((aligned(16))) cfloat3 {
  float x, y, z;
};
struct __attribute__((aligned(16))) cfloat4 {
  float x, y, z, w;
};
struct __attribute__((aligned(32))) cfloat5 {
  float x, y, z, w, h;
};
inline cuchar2 cmake_uchar2(uchar x, uchar y) {
  cuchar2 t;
  t.x = x;
  t.y = y;
  return t;
}
inline cuchar3 cmake_uchar3(cuchar4 orig) {
  cuchar3 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  return t;
}
inline cuchar3 cmake_uchar3(uchar x, uchar y, uchar z) {
  cuchar3 t;
  t.x = x;
  t.y = y;
  t.z = z;
  return t;
}
inline cuchar4 cmake_uchar4(uchar x, uchar y, uchar z, uchar w) {
  cuchar4 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  return t;
}
inline cint2 cmake_int2(int x, int y) {
  cint2 t;
  t.x = x;
  t.y = y;
  return t;
}
inline cfloat3 cmake_float3(float x, float y, float z) {
  cfloat3 t;
  t.x = x;
  t.y = y;
  t.z = z;
  return t;
}
inline cfloat4 cmake_float4(float x, float y, float z, float w) {
  cfloat4 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  return t;
}
inline cfloat4 cmake_float4(cfloat5 orig) {
  cfloat4 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  return t;
}
inline cfloat5 cmake_float5(float x, float y, float z, float w, float h) {
  cfloat5 t;
  t.x = x;
  t.y = y;
  t.z = z;
  t.w = w;
  t.h = h;
  return t;
}
inline cfloat5 cmake_float5(cfloat4 orig, float h) {
  cfloat5 t;
  t.x = orig.x;
  t.y = orig.y;
  t.z = orig.z;
  t.w = orig.w;
  t.h = h;
  return t;
}
inline cfloat3 cplus3(cfloat4 A, cfloat4 B) {
  return cmake_float3(A.x + B.x, A.y + B.y, A.z + B.z);
}
inline cfloat3 cplus3(cfloat3 A, cfloat3 B) {
  return cmake_float3(A.x + B.x, A.y + B.y, A.z + B.z);
}
inline cfloat3 cdivide3(cfloat3 A, float s) {
  return cmake_float3(A.x / s, A.y / s, A.z / s);
}
inline float cdet2x2(float a11, float a12, float a21, float a22) {
  return a11 * a22 - a12 * a21;
}
inline float cdet3x3(float a11, float a12, float a13, float a21, float a22,
                     float a23, float a31, float a32, float a33) {
  return a11 * cdet2x2(a22, a23, a32, a33) - a21 * cdet2x2(a12, a13, a32, a33) +
         a31 * cdet2x2(a12, a13, a22, a23);
}

//----------------------------------------------------------------------------
// 4 faces of tet abcd: cbd acd bad abc (in vids)
constexpr int tet_faces_lvid_host[4][3] = {
    {2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
// 6 edges of tet (in vids)
constexpr int tet_edges_lvid_host[6][2] = {{2, 3}, {1, 3}, {1, 2},
                                           {0, 3}, {0, 2}, {0, 1}};
// 4 vertices of tet (in lfids)
constexpr int tet_vs_lfid_host[4][3] = {
    {1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};
// 6 edges of tet (in lfids)
constexpr int tet_edges_lfid_host[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                           {1, 2}, {1, 3}, {2, 3}};

// keep the same code as get_edge_idx() in convex_cell.h
inline int get_edge_idx_copy(int v1, int v2, int n) {
  int vmin = v1;
  int vmax = v2;
  if (v1 > v2) {
    vmin = v2;
    vmax = v1;
  }
  int idx = 0;
  for (uint j = 0; j <= vmin; j++) {
    idx += (n - j);
  }
  idx -= (n - vmax);
  return idx;
}

//----------------------------------------------------------------------------

template <typename T>
struct is_array_or_vector {
  enum { value = false };
};

template <typename T, typename A>
struct is_array_or_vector<std::vector<T, A>> {
  enum { value = true };
};

template <typename T, std::size_t N>
struct is_array_or_vector<std::array<T, N>> {
  enum { value = true };
};

template <class T>
inline void swap_in_place(T& x, T& y) {
  T temp;
  temp = x;
  x = y;
  y = temp;
}

// binary search function using template
// x: value to search
// the fucntion returns -1 if x is not found in arr
// otherwise it returns index of x
template <typename T>
inline int binary_search(std::vector<T> arr, T x) {
  int n = arr.size();
  int start = 0;
  int end = n - 1;
  while (start <= end) {
    int mid = (start + end) / 2;
    if (arr[mid] == x)
      return mid;
    else if (arr[mid] < x)
      start = mid + 1;
    else
      end = mid - 1;
  }
  return -1;
}

inline std::string get_timestamp() {
  auto now = std::time(nullptr);
  std::ostringstream os;
  // os << std::put_time(std::gmtime(&now),"%F  %T");
  os << std::put_time(std::localtime(&now), "%F_%T");
  return os.str();
}

inline bool is_file_exist(const std::string filePath) {
  std::ifstream infile(filePath);
  return infile.good();
}

inline std::string get_file_ext(const std::string filePath) {
  return filePath.substr(filePath.find_last_of(".") + 1);
}

inline std::string get_file_no_ext(std::string filePath) {
  std::string filename = filePath.substr(0, filePath.find_last_of("."));
  // ftetwild save format xxx.ply_.msh
  filename = filePath.substr(0, filename.find_last_of("."));
  return filename;
}

inline std::string get_only_file_name(std::string filePath, bool withExtension,
                                      char seperator = '/') {
  std::string filename_ext =
      filePath.substr(filePath.find_last_of(seperator) + 1);
  if (withExtension) return filename_ext;
  size_t lastindex = filename_ext.find_last_of(".");
  return filename_ext.substr(0, lastindex);
}

// namespace fs = std::filesystem;
// inline bool create_dir(const std::string& dir) {
//   // if dir not exists, create it
//   if (!fs::is_directory(dir.c_str())) {
//     return fs::create_directories(dir.c_str());
//   }
//   return true;
// }

// using .geogram for saving vertex attributes
//
// (NO) matching ftetwild output (https://github.com/wildmeshing/fTetWild)
//
// type:
// 0 -> .geogram
// 1 -> _sf.obj
// 2 -> _pts.xyz
// 3 -> _sf_01_scaled.obj // matching matfp?
// 4 -> _extf.ma
// 5 -> _pts_fid.xyz
inline std::string get_other_file_path(std::string filePath, int type) {
  std::string file_path = get_file_no_ext(filePath);
  switch (type) {
    case 0:
      return file_path + "_sf.geogram";
      break;
    case 1:
      return file_path + "_sf.obj";
    case 2:
      return file_path + "_pts.xyz";
    case 3:  // no use
      return file_path + "_sf_01_scaled.obj";
    case 4:
      return file_path + "_extf.ma";
    case 5:
      return file_path + "_pts_fid.xyz";
    default:
      break;
  }
  assert(false);
}

template <typename T>
inline void vector_unique(std::vector<T>& v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

template <typename T>
inline T get_sorted(const T& in) {
  T out = in;
  std::sort(out.begin(), out.end());
  return out;
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      std::set<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  for (uint x : s1) {
    if (s2.count(x)) {
      v.insert(x);
    }
  }
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  if (s1.empty() || s2.empty()) return;
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      const std::set<T>& s3, std::set<T>& v) {
  v.clear();
  std::set<T> neigh_cells_tmp;
  set_intersection<T>(s1, s2, neigh_cells_tmp);
  if (neigh_cells_tmp.empty()) return;
  set_intersection<T>(neigh_cells_tmp, s3, v);
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2, std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2,
                      const std::unordered_set<T>& s3, std::vector<T>& v) {
  if (s2.size() < s1.size() && s2.size() < s1.size()) {
    set_intersection(s2, s1, s3, v);
    return;
  }

  if (s3.size() < s1.size() && s3.size() < s2.size()) {
    set_intersection(s3, s1, s2, v);
    return;
  }

  assert(s1.size() <= s2.size());
  assert(s1.size() <= s3.size());

  v.clear();
  v.reserve(s1.size());
  for (int x : s1) {
    if (s2.count(x) && s3.count(x)) {
      v.push_back(x);
      if (v.size() == 2) break;
    }
  }
}

template <typename T>
bool has_intersection(const std::set<T>& s1, const std::set<T>& s2) {
  if (s1.empty() || s2.empty()) return false;
  for (uint x : s1) {
    if (s2.count(x)) {
      return true;
    }
  }
  return false;
}

template <typename T>
std::set<T> to_set(const std::vector<T>& in) {
  std::set<T> out;
  for (int x : in) {
    out.insert(x);
  }
  return out;
}

inline int mod3(int j) { return j % 3; }
inline int mod4(int j) { return j % 4; }
inline int modk(int j, int k) { return j % k; }

/**
 * @brief Given cells/others to trace, and return the number of connected
 * components
 *
 * @param cell_to_visit
 * @param cell_neighbors
 * @param cc_cells return, grouped in each CC
 * @return int number of CC
 */
template <typename T>
inline int get_CC_given_neighbors(
    const std::set<T>& cell_to_visit,
    const std::map<T, std::set<T>>& cell_neighbors,
    std::vector<std::set<T>>& cc_cells) {
  std::set<T> cell_unvisited = cell_to_visit;  // copy
  cc_cells.clear();

  // printf("cell_unvisited: %ld, cell_neighbors size: %ld \n",
  //        cell_unvisited.size(), cell_neighbors.size());
  int num_cc = 0;
  std::set<T> cell_visited, one_cc_visited;
  // calculate the number of CCs
  std::queue<T> cc_queue;
  while (!cell_unvisited.empty()) {
    cc_queue.push(*cell_unvisited.begin());
    while (!cc_queue.empty()) {
      T cell_id = cc_queue.front();
      cc_queue.pop();

      if (cell_visited.find(cell_id) != cell_visited.end()) continue;
      cell_visited.insert(cell_id);    // all visited cells
      cell_unvisited.erase(cell_id);   // all unvisited cells
      one_cc_visited.insert(cell_id);  // visited cells in one CC
      // push all neighbors
      if (cell_neighbors.find(cell_id) == cell_neighbors.end()) {
        // assert(false);
        continue;
      }
      auto& neighbors = cell_neighbors.at(cell_id);
      for (auto& neigh_id : neighbors) {
        // if not in cell_to_visit, do not add to queue
        if (cell_to_visit.find(neigh_id) != cell_to_visit.end())
          cc_queue.push(neigh_id);
      }
    }
    // printf("found one CC with cell_visited: %zu/%zu, cell_unvisited:
    // %zu/%zu\n",
    //        cell_visited.size(), cell_to_visit.size(), cell_unvisited.size(),
    //        cell_to_visit.size());
    num_cc++;
    cc_cells.push_back(one_cc_visited);  // store cell ids in one CC
    one_cc_visited.clear();
  }

  return num_cc;
}

template <typename T>
inline void print_set(const std::set<T>& to_print,
                      const std::string name = "") {
  printf("%s: [", name.c_str());
  for (const auto& v : to_print) {
    std::cout << v << ", ";
  }
  printf("]\n");
}

template <typename T>
inline void print_vec(const std::vector<T>& to_print) {
  printf("[");
  for (const auto& v : to_print) {
    std::cout << v << ", ";
  }
  printf("]\n");
}

// Define a upper triangular matrix with size n.
// matrix index -> (v1_min, v2_max)
// n -> #vertices
// idx = idx(v1_min, v2_max)
//     = n + (n-1) + ... + (n-vid_min) - (n - vid_max)
inline int get_upper_tri_matrix_size(const int n) {
  return n * (1 + n) / 2 + 1;
}
inline int get_upper_tri_matrix_idx(int v1, int v2, int n) {
  int vmin = v1;
  int vmax = v2;
  if (v1 > v2) {
    vmin = v2;
    vmax = v1;
  }
  int idx = 0;
  for (uint j = 0; j <= vmin; j++) {
    idx += (n - j);
  }
  idx -= (n - vmax);
  return idx;
}
