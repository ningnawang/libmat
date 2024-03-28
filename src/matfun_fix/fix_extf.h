#pragma once

#include "common_geogram.h"
#include "fix_common.h"
#include "medial_mesh.h"
#include "voronoi_defs.h"

int check_and_fix_external_feature(
    const int num_itr_global, const Parameter& param,
    const std::vector<float>& tet_vertices,
    const std::vector<ConvexCellHost>& convex_cells_host,
    const GEO::Mesh& sf_mesh, const std::map<aint4, int>& tet_vs_lfs2tvs_map,
    const std::map<int, std::set<int>>& fl2corner_sphere,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug);