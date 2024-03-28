#pragma once

#include "common_geogram.h"
#include "input_types.h"
#include "medial_mesh.h"
#include "medial_primitives.h"
#include "medial_sphere.h"
#include "updating.h"

// only add one sphere if any to fix
bool random_check_edges_and_fix_internal_feature(
    const int num_itr_global, const SurfaceMesh& sf_mesh,
    const TetMesh& tet_mesh, const MedialMesh& mmesh,
    std::vector<MedialSphere>& all_medial_spheres,
    std::set<aint2>& checked_mmesh_edges, bool is_debug);

// [no use]
int delete_T2_if_TN(const int num_itr_global, const SurfaceMesh& sf_mesh,
                    const TetMesh& tet_mesh,
                    std::vector<MedialSphere>& all_medial_spheres,
                    bool is_debug);