#pragma once

#include "common.h"

#define FOR(I, UPPERBND) for (int I = 0; I < int(UPPERBND); ++I)

// must be uint, matching site_flag
enum SiteFlag { no_flag = 0, is_selected = 1 };

enum Status {
  // For cases:
  // 1. tid >= n_tet
  // 2. seed NOT in [0, n_site)
  // 3. un-success clip
  early_return = -1,  // ninwang,
  triangle_overflow = 0,
  vertex_overflow = 1,
  inconsistent_boundary = 2,
  security_radius_not_reached = 3,
  success = 4,
  needs_exact_predicates = 5,
  no_intersection = 6,
  edge_overflow = 7,
  needs_perturb = 8  // ninwang
};

// Uncomment to activate arithmetic filters.
//   If arithmetic filters are activated,
//   status is set to needs_exact_predicates
//   whenever predicates could not be evaluated
//   using floating points on the GPU
// #define USE_ARITHMETIC_FILTER

#define PRESET 3

#if 0 == PRESET  // conservative settings (white noise)
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 32
#define _K_ 90  // k_site nearest sites of each site
// #define _TET_K_         70
#define _MAX_P_ 64
#define _MAX_T_ 96
#define _MAX_E_ 152
#elif 1 == PRESET  // perturbed grid settings
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 32
#define _K_ 90
#define _MAX_P_ 50
#define _MAX_T_ 96
#elif 2 == PRESET  // blue noise settings
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 64
#define _K_ 35
#define _MAX_P_ 32
#define _MAX_T_ 96
#elif 3 == PRESET  // ninwang setting
#define VORO_BLOCK_SIZE 16
#define KNN_BLOCK_SIZE 32  // no use
#define _K_ 150            // for testing
#define _MAX_P_ 64
#define _MAX_T_ 96
#define _MAX_E_ 152
#endif

#define IF_VERBOSE(x) x