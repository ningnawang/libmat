#pragma once

#include "medial_mesh.h"
#include "medial_sphere.h"

// for writing to file
// [no use]
void sample2slab_cone_and_write(
    const SurfaceMesh& sf_mesh, const double sampling_step,
    const std::vector<MedialSphere>& all_medial_spheres,
    const MedialMesh& mmesh, const std::string filename);

void sample2prims_and_compute_dist2mat(
    const SurfaceMesh& sf_mesh, const double sampling_step,
    const std::vector<MedialSphere>& all_medial_spheres,
    const MedialMesh& mmesh, std::vector<double>& samples,
    std::vector<int>& sample_fids, std::vector<float>& samples_dist2mat,
    std::vector<aint3>& samples_clostprim, bool is_filter_se);