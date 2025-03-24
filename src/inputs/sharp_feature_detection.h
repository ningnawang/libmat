#pragma once

#include <unordered_set>
#include <vector>

#include "common_geogram.h"
#include "input_types.h"
#include "params.h"

enum EdgeConvexConcave { FLAT = -1, CONVEX = 1, CONCAVE = 2 };

int get_grouped_feature_edges(const std::set<aint2>& fe_to_visit,
                              const std::set<int>& corners,
                              std::vector<std::vector<aint2>>& fe_groups);

void detect_mark_sharp_features(const Parameter& args, SurfaceMesh& sf_mesh,
                                TetMesh& tet_mesh, bool is_sf_file_exist);