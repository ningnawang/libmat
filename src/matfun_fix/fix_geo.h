#pragma once

#include "common_geogram.h"
#include "fix_topo.h"
#include "medial_mesh.h"
#include "medial_primitives.h"
#include "medial_sphere.h"
#include "shrinking.h"

// No use
void check_mm_geo(const GEO::Mesh& sf_mesh, const Parameter& params,
                  MedialMesh& mmesh,
                  std::vector<MedialSphere>& all_medial_spheres,
                  std::set<int>& slabs_to_fix, bool is_debug);
// No use
void fix_mm_by_adding_new_sphere(const SurfaceMesh& sf_mesh,
                                 const std::set<int>& slabs_to_fix,
                                 const MedialMesh& mmesh,
                                 std::vector<MedialSphere>& all_medial_spheres,
                                 bool is_debug);

void check_and_fix_mm_geo(const int num_itr_global, const SurfaceMesh& sf_mesh,
                          const TetMesh& tet_mesh, const Parameter& params,
                          const MedialMesh& mmesh,
                          std::vector<MedialSphere>& all_medial_spheres,
                          std::vector<double>& samples,
                          std::vector<float>& samples_dist2mat,
                          std::vector<aint3>& samples_clostprim, bool is_debug);