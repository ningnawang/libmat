#pragma once
#include <assert.h>
#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "medial_sphere.h"
#include "params.h"
#include "voronoi_defs.h"

bool save_convex_cells_houdini(
    const Parameter params, const std::vector<MedialSphere>& all_medial_spheres,
    const std::vector<ConvexCellHost>& convex_cells_returned,
    std::string rpd_name, const int max_sf_fid, const bool is_boundary_only,
    bool is_slice_plane = false);
