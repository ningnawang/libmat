
#include "fix_geo_error.h"

#include <geogram/mesh/mesh_sampling.h>

#include <algorithm>

#include "dist2mat.h"

namespace {

using namespace GEO;

/**
 * \brief Generates a set of random samples over a surfacic mesh.
 * \param[in] mesh the mesh
 * \param[out] p pointer to an array of generated samples, of size
 *   \p nb_points times DIM. To be allocated by the caller.
 * \param[out] f pointer to an array of fids, of size
 *   \p nb_points. To be allocated by the caller.
 * \param[in] nb_points number of points to generate
 * \param[in] weight a reference to a vertex attribute. If bound, it
 *  is taken into account.
 * \param[in] facets_begin_in if specified, first index of the facet
 *  sequence in which points should be generated. If left unspecified (-1),
 *  points are generated over all the facets of the mesh.
 * \param[in] facets_end_in if specified, one position past the last
 *  index of the facet sequence in which points should be generated.
 *  If left unspecified (-1), points are generated over all the facets
 *  of the mesh.
 * \tparam DIM dimension of the points, specified as a template argument
 *  for efficiency reasons
 * \return true if everything went OK, false otherwise. Whenever all the
 *  points land in the same facet, the function returns false to notify
 *  a potential numerical problem.
 */
template <GEO::index_t DIM>
inline bool mesh_generate_random_samples_on_surface(
    const GEO::Mesh& mesh, double* p, int* f, GEO::index_t nb_points,
    GEO::Attribute<double>& weight, GEO::signed_index_t facets_begin_in = -1,
    GEO::signed_index_t facets_end_in = -1) {
  geo_assert(mesh.facets.are_simplices());
  geo_assert(mesh.vertices.dimension() >= DIM);
  geo_assert(mesh.facets.nb() > 0);

  index_t facets_begin = 0;
  index_t facets_end = mesh.facets.nb();
  if (facets_begin_in != -1) {
    facets_begin = index_t(facets_begin_in);
  }
  if (facets_end_in != -1) {
    facets_end = index_t(facets_end_in);
  }

  typedef vecng<DIM, double> Point;

  // To ensure reproducibility accros successive
  // runs, reset the random number generator.
  Numeric::random_reset();

  vector<double> s(nb_points);
  for (index_t i = 0; i < nb_points; i++) {
    s[i] = Numeric::random_float64();
  }
  std::sort(s.begin(), s.end());

  double Atot = 0.0;
  for (index_t t = facets_begin; t < facets_end; ++t) {
    double At = mesh_facet_mass<DIM>(mesh, t, weight);
    Atot += At;
  }

  signed_index_t first_t = -1;
  signed_index_t last_t = 0;

  index_t cur_t = facets_begin;
  double cur_s = mesh_facet_mass<DIM>(mesh, facets_begin, weight) / Atot;
  for (index_t i = 0; i < nb_points; i++) {
    geo_debug_assert(i < s.size());
    while (s[i] > cur_s && cur_t < facets_end - 1) {
      cur_t++;
      geo_debug_assert(cur_t < facets_end);
      cur_s += mesh_facet_mass<DIM>(mesh, cur_t, weight) / Atot;
    }
    if (first_t == -1) {
      first_t = signed_index_t(cur_t);
    }
    last_t = std::max(last_t, signed_index_t(cur_t));

    // TODO: take weights into account
    //  with a new random_point_in_triangle_weighted()
    //  function.
    index_t v1 = mesh.facets.vertex(cur_t, 0);
    index_t v2 = mesh.facets.vertex(cur_t, 1);
    index_t v3 = mesh.facets.vertex(cur_t, 2);
    Point cur_p = Geom::random_point_in_triangle(
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v1)),
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v2)),
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v3)));
    for (coord_index_t coord = 0; coord < DIM; coord++) {
      p[i * DIM + coord] = cur_p[coord];
    }
    f[i] = cur_t;
  }
  if (mesh.facets.nb() > 1 && last_t == first_t) {
    Logger::warn("Sampler")
        << "Did put all the points in the same triangle" << std::endl;
    return false;
  }
  return true;
}
}  // namespace

/****************************************************************************/
void sample_points_on_mesh(const SurfaceMesh& sf_mesh,
                           const double sampling_step,
                           std::vector<double>& samples,
                           std::vector<int>& sample_fids, bool is_filter_se) {
  GEO::index_t nb_samples =
      GEO::index_t(Geom::mesh_area(sf_mesh, 3) / GEO::geo_sqr(sampling_step));
  nb_samples = std::max(nb_samples, GEO_SAMPLE_MAX);

  printf("Sampling %d surface samples...\n", nb_samples);
  samples.resize(nb_samples * 3);
  sample_fids.resize(nb_samples);
  GEO::Attribute<double> weight;  // left unbound
  mesh_generate_random_samples_on_surface<3>(
      sf_mesh, samples.data(), sample_fids.data(), nb_samples, weight);

  if (is_filter_se) {
    std::vector<double> filterd_samples;
    std::vector<int> filterd_samples_fids;
    for (int i = 0; i < samples.size() / 3; i++) {
      Vector3 p(samples.at(i * 3), samples.at(i * 3 + 1),
                samples.at(i * 3 + 2));
      double sq_dist2se = sf_mesh.aabb_wrapper.get_sq_dist_to_se(p);
      // if (sq_dist2se < sampling_step) continue;  // skip
      if (sq_dist2se < SCALAR_1) continue;  // skip
      filterd_samples.push_back(p[0]);
      filterd_samples.push_back(p[1]);
      filterd_samples.push_back(p[2]);
      filterd_samples_fids.push_back(sample_fids.at(i));
    }
    samples = filterd_samples;
    sample_fids = filterd_samples_fids;
  }
}

void gather_point_to_sites(const std::vector<MedialSphere>& all_medial_spheres,
                           const std::vector<double>& samples,
                           const std::vector<int>& sample_fids,
                           std::map<int, std::set<int>>& sample2sites) {
  int dim = 3;
  // load from site2fids to fid2sites
  std::map<int, std::set<int>> fid2sites;
  for (const auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) continue;
    for (const auto& pair : msphere.pcell.cell_to_surfv2fid) {
      for (const auto& v2fid : pair.second) {
        fid2sites[v2fid.second].insert(msphere.id);
      }
    }
  }

  // load from fid2site to sample2sites
  sample2sites.clear();
  for (int sample_id = 0; sample_id < samples.size() / dim; sample_id++) {
    int sfid = sample_fids.at(sample_id);
    // TOOD: why? needs to debug
    if (fid2sites.find(sfid) == fid2sites.end())
      sample2sites[sample_id] = std::set<int>();  // empty
    else
      sample2sites[sample_id] = fid2sites.at(sfid);
    // printf("[Sample2SCS] sample %d has fid %d clipped by %zu sites\n",
    //        sample_id, sfid, sample2sites.at(sample_id).size());
  }
}

void gather_point_to_slab_and_cone(
    const std::vector<MedialSphere>& all_medial_spheres,
    const MedialMesh& mmesh, const std::map<int, std::set<int>>& sample2sites,
    std::map<int, std::vector<aint3>>& sample2slab_cone) {
  auto fill_sample_slab = [&](const bool& is_cone, const int id) {
    if (is_cone) {
      aint2 vs = mmesh.edges.at(id).vertices_;
      aint3 result = {{-1, vs[0], vs[1]}};
      return result;
    }
    return mmesh.faces.at(id).vertices_;
  };

  sample2slab_cone.clear();
  std::vector<aint3> tmp;
  for (const auto& pair : sample2sites) {
    int sample_id = pair.first;
    const auto& sample_sites = pair.second;
    tmp.clear();
    // add slab + cone
    for (const auto& site_id : sample_sites) {
      for (const auto& mfid : all_medial_spheres.at(site_id).faces_)
        tmp.push_back(fill_sample_slab(false, mfid));
      for (const auto& meid : all_medial_spheres.at(site_id).edges_)
        tmp.push_back(fill_sample_slab(true, meid));
      // add sphere
      tmp.push_back({{-1, -1, site_id}});
    }
    // add to sample2slab_cone
    sample2slab_cone[sample_id] = tmp;
    // printf("sample %d has prims: [");
    // for (const auto& prim : tmp) {
    //   printf("(%d,%d,%d), ", prim[0], prim[1], prim[2]);
    // }
    // printf("]\n");
  }
}

// [format]
// #spheres
// x y z r
// #samples
// x y x #slab+cone+sphere
// v1 v2 v3 (slab)
// -1 v1 v2 (cone)
// -1 -1 v1 (sphere)
void write_sample2slab_cone_to_file(
    const std::vector<MedialSphere>& all_medial_spheres,
    const std::vector<double>& samples,
    const std::map<int, std::vector<aint3>>& sample2slab_cone,
    const std::string filename) {
  std::string cells_path =
      "../out/geo/sample_" + filename + "_" + get_timestamp() + ".samples";
  int dim = 3;
  int n_site = all_medial_spheres.size();
  int n_samples = samples.size() / dim;

  std::fstream file;
  file.open(cells_path, std::ios_base::out);
  // [format]
  // #spheres
  // x y z r
  file << n_site << std::endl;
  for (const auto& msphere : all_medial_spheres) {
    file << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << msphere.center[0] << " " << msphere.center[1] << " "
         << msphere.center[2] << " " << msphere.radius << std::endl;
  }
  // [format]
  // #samples
  // x y x #slab+cone+sphere
  file << n_samples << std::endl;
  FOR(sid, n_samples) {
    const auto& slab_cone_sphere = sample2slab_cone.at(sid);
    file << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << samples[sid * dim] << " " << samples[sid * dim + 1] << " "
         << samples[sid * dim + 2] << " " << slab_cone_sphere.size()
         << std::endl;
    // [format]
    // v1 v2 v3 (slab)
    // -1 v1 v2 (cone)
    // -1 -1 v1 (sphere)
    for (const auto& scs : slab_cone_sphere) {
      file << scs[0] << " " << scs[1] << " " << scs[2] << std::endl;
    }
  }
}

// for saving '.sample' file
void sample2slab_cone_and_write(
    const SurfaceMesh& sf_mesh, const double sampling_step,
    const std::vector<MedialSphere>& all_medial_spheres,
    const MedialMesh& mmesh, const std::string filename) {
  std::vector<double> samples;
  std::vector<int> sample_fids;
  sample_points_on_mesh(sf_mesh, sampling_step, samples, sample_fids,
                        false /*is_filter_se*/);
  int n_samples = samples.size() / 3;
  assert(n_samples == sample_fids.size());
  printf("[Sample2SCS] sampled %d on sf_mesh \n", n_samples);

  std::map<int, std::set<int>> sample2sites;
  gather_point_to_sites(all_medial_spheres, samples, sample_fids, sample2sites);
  assert(sample2sites.size() == n_samples);

  // stores mspheres whose pcells covers sample fid
  std::map<int, std::vector<aint3>> sample2slab_cone;
  gather_point_to_slab_and_cone(all_medial_spheres, mmesh, sample2sites,
                                sample2slab_cone);
  assert(sample2slab_cone.size() == n_samples);

  printf("[Sample2SCS] sample2slab_cone %zu \n", sample2slab_cone.size());

  write_sample2slab_cone_to_file(all_medial_spheres, samples, sample2slab_cone,
                                 filename);
}

// -------------------------------------------------------------------------------------------------

#include "dist2mat_cuda_utils.h"
// call GPU to compute the distance
void load_and_compute_sample_dist2mat_gpubuffer(
    const std::vector<MedialSphere>& all_medial_spheres,
    const std::vector<double>& samples_double /*dim=3*/,
    const std::map<int, std::vector<aint3>>& sample2slab_cone,
    std::vector<float>& samples_dist2mat, std::vector<aint3>& samples_clostprim,
    bool is_debug = false) {
  int num_spheres = all_medial_spheres.size();
  int num_samples = samples_double.size() / 3;

  int total_prim_num = 0;
  GpuBuffer<float4> spheres;
  GpuBuffer<float3> samples;
  GpuBuffer<uint> offset;
  GpuBuffer<uint> num_per_sample;
  GpuBuffer<int3> prims;
  GpuBuffer<int> closest_mat_id;
  GpuBuffer<float> results;

  std::vector<int3> prim_per_samples;  // temporary

  // load spheres
  spheres.HResize(num_spheres);
  for (int i = 0; i < num_spheres; ++i) {
    const auto& mspheres = all_medial_spheres.at(i);
    float4 s;
    s.x = mspheres.center[0];
    s.y = mspheres.center[1];
    s.z = mspheres.center[2];
    s.w = mspheres.radius;
    spheres.HPtr()[i] = s;
  }

  samples.HResize(num_samples);
  results.HResize(num_samples);
  offset.HResize(num_samples);
  closest_mat_id.HResize(num_samples);
  num_per_sample.HResize(num_samples);
  for (int i = 0; i < num_samples; ++i) {
    float3 pos;
    pos.x = samples_double.at(i * 3);
    pos.y = samples_double.at(i * 3 + 1);
    pos.z = samples_double.at(i * 3 + 2);
    samples.HPtr()[i] = pos;
    results.HPtr()[i] = 1e28f;
    const auto& sample_prims = sample2slab_cone.at(i);
    uint num_prims = sample_prims.size();
    num_per_sample.HPtr()[i] = num_prims;
    offset.HPtr()[i] = total_prim_num;
    total_prim_num += num_prims;

    for (int j = 0; j < num_prims; ++j) {
      int3 prim;
      prim.x = sample_prims.at(j)[0];
      prim.y = sample_prims.at(j)[1];
      prim.z = sample_prims.at(j)[2];
      prim_per_samples.push_back(prim);
    }
    closest_mat_id.HPtr()[i] = -1;
  }

  prims.HResize(prim_per_samples.size());
  for (size_t i = 0; i < prim_per_samples.size(); ++i) {
    prims.HPtr()[i] = prim_per_samples[i];
  }

  compute_closest_dist2mat(spheres, num_samples, samples, offset,
                           num_per_sample, prims, results, closest_mat_id);

  samples_dist2mat.resize(num_samples);
  samples_clostprim.resize(num_samples);
  for (int i = 0; i < num_samples; i++) {
    int clost_prim_id = closest_mat_id.HPtr()[i];
    float clost_dist = results.HPtr()[i];
    // TODO: error? why this happens？
    if (clost_prim_id == -1) {
      samples_dist2mat[i] = 0.f;
      samples_clostprim[i] = {{-1, -1, -1}};
      continue;
    }
    samples_dist2mat[i] = clost_dist;
    const aint3& prim = sample2slab_cone.at(i).at(clost_prim_id);
    samples_clostprim[i] = prim;
    if (is_debug)
      printf("sample %d has closest prim %d (%d,%d,%d) with dist %f\n", i,
             clost_prim_id, prim[0], prim[1], prim[2], clost_dist);
  }
  if (is_debug) printf("done load_and_compute_sample_dist2mat_gpubuffer\n");
}

void sample2prims_and_compute_dist2mat(
    const SurfaceMesh& sf_mesh, const double sampling_step,
    const std::vector<MedialSphere>& all_medial_spheres,
    const MedialMesh& mmesh, std::vector<double>& samples,
    std::vector<int>& sample_fids, std::vector<float>& samples_dist2mat,
    std::vector<aint3>& samples_clostprim, bool is_filter_se) {
  samples.clear();
  sample_fids.clear();
  sample_points_on_mesh(sf_mesh, sampling_step, samples, sample_fids,
                        is_filter_se);
  int n_samples = samples.size() / 3;
  assert(n_samples == sample_fids.size());
  printf("[Sample2SCS] sampled %d on sf_mesh \n", n_samples);

  std::map<int, std::set<int>> sample2sites;
  gather_point_to_sites(all_medial_spheres, samples, sample_fids, sample2sites);
  assert(sample2sites.size() == n_samples);
  printf("[Sample2SCS] sample2sites: %zu\n", sample2sites.size());

  // stores mspheres whose pcells covers sample fid
  std::map<int, std::vector<aint3>> sample2slab_cone;
  gather_point_to_slab_and_cone(all_medial_spheres, mmesh, sample2sites,
                                sample2slab_cone);
  assert(sample2slab_cone.size() == n_samples);

  printf("[Sample2SCS] sample2slab_cone %zu \n", sample2slab_cone.size());

  samples_dist2mat.clear();
  samples_clostprim.clear();
  load_and_compute_sample_dist2mat_gpubuffer(
      all_medial_spheres, samples, sample2slab_cone, samples_dist2mat,
      samples_clostprim, false /*is_debug*/);
}
