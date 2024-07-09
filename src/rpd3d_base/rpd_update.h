#pragma once

#include "input_types.h"
#include "medial_sphere.h"
#include "voronoi_defs.h"

void update_power_cells(const SurfaceMesh& sf_mesh,
                        std::vector<ConvexCellHost>& convex_cells_host,
                        std::vector<MedialSphere>& all_medial_spheres,
                        const std::map<aint3, aint3>& tet_es2fe_map,
                        bool is_debug);