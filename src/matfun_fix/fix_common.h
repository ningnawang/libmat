#pragma once

#include "common_geogram.h"
#include "medial_sphere.h"
#include "sharp_feature_detection.h"
#include "voronoi_defs.h"

void convert_fes_from_ccell_to_tet(
    const MedialSphere& msphere, const std::vector<float>& tet_vertices,
    const std::vector<ConvexCellHost>& convex_cells_host,
    const std::map<aint4, int>& tet_vs_lfs2tvs_map,
    const std::set<aint5>& se_covered_lvids, std::set<aint4>& fe_tvs_all,
    std::map<int, Vector3>& tvs_pos_map, bool is_debug = false);

int get_grouped_fes_in_order(const std::set<aint4>& fe_to_visit,
                             std::vector<std::vector<aint4>>& fe_groups);