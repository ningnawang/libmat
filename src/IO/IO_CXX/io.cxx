
#include "io.h"

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>

#include <bitset>

#include "../extern/mshloader/MshLoader.h"
// #include "voronoi_common.h"

bool is_inverted(const Vector3& v0, const Vector3& v1, const Vector3& v2,
                 const Vector3& v3) {
  if (Predicates::orient_3d(v0, v1, v2, v3) == Predicates::ORI_POSITIVE)
    return false;
  return true;
}

void get_bbox(const std::vector<float>& vertices, float& xmin, float& ymin,
              float& zmin, float& xmax, float& ymax, float& zmax,
              float& bbox_diag_l) {
  int nb_v = vertices.size() / 3;
  xmin = xmax = vertices[0];
  ymin = ymax = vertices[1];
  zmin = zmax = vertices[2];
  for (int i = 1; i < nb_v; ++i) {
    xmin = std::min(xmin, vertices[3 * i]);
    ymin = std::min(ymin, vertices[3 * i + 1]);
    zmin = std::min(zmin, vertices[3 * i + 2]);
    xmax = std::max(xmax, vertices[3 * i]);
    ymax = std::max(ymax, vertices[3 * i + 1]);
    zmax = std::max(zmax, vertices[3 * i + 2]);
  }
  float d = xmax - xmin;
  d = std::max(d, ymax - ymin);
  d = std::max(d, zmax - zmin);
  d = 0.001f * d;
  xmin -= d;
  ymin -= d;
  zmin -= d;
  xmax += d;
  ymax += d;
  zmax += d;

  bbox_diag_l = std::sqrt(std::pow(xmax - xmin, 2) + std::pow(ymax - ymin, 2) +
                          std::pow(zmax - zmin, 2));
}

int get_tet_euler(const std::vector<float>& tet_vertices,
                  const std::vector<int>& tet_indices) {
  std::set<aint2> tet_edges;
  std::set<aint3> tet_faces;

  int faces[4][3] = {{2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
  int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

  bool is_debug = false;

  for (int tid = 0; tid < tet_indices.size() / 4; tid++) {
    // if (tid == 0)
    //   is_debug = true;
    // else
    //   is_debug = false;
    if (is_debug)
      printf("tet %d has vertices (%d,%d,%d,%d)\n", tid, tet_indices[tid * 4],
             tet_indices[tid * 4 + 1], tet_indices[tid * 4 + 2],
             tet_indices[tid * 4 + 3]);

    // store faces
    for (int i = 0; i < 4; i++) {
      aint3 one_face = {{tet_indices[tid * 4 + faces[i][0]],
                         tet_indices[tid * 4 + faces[i][1]],
                         tet_indices[tid * 4 + faces[i][2]]}};
      std::sort(one_face.begin(), one_face.end());
      tet_faces.insert(one_face);
      if (is_debug)
        printf("tet %d has face (%d,%d,%d) \n", tid, one_face[0], one_face[1],
               one_face[2]);
    }

    //  store edges
    for (int j = 0; j < 6; j++) {
      aint2 one_edge = {{tet_indices[tid * 4 + edges[j][0]],
                         tet_indices[tid * 4 + edges[j][1]]}};
      std::sort(one_edge.begin(), one_edge.end());
      tet_edges.insert(one_edge);
      if (is_debug)
        printf("tet %d has edge (%d,%d) \n", tid, one_edge[0], one_edge[1]);
    }
  }

  int euler = tet_vertices.size() / 3 - tet_edges.size() + tet_faces.size() -
              tet_indices.size() / 4;

  std::cout << "tet has v: " << tet_vertices.size() / 3
            << ", e: " << tet_edges.size() << ", f: " << tet_faces.size()
            << ", t: " << tet_indices.size() / 4 << ", and euler: " << euler
            << std::endl;
}

// vertices: tet vertices
// indices: tet 4 indices of vertices
using namespace PyMesh;
bool load_tet(const std::string& filename, std::vector<float>& vertices,
              std::vector<int>& indices, bool normalize, Parameter& params) {
  std::string s;
  int n_vertex, n_tet, temp;

  std::ifstream input(filename);
  if (input.fail()) return false;

  std::string ext = filename.substr(filename.find_last_of('.') + 1);
  if (ext == "msh") {
    MshLoader msh_loader(filename);
    vertices = msh_loader.get_nodes();
    indices = msh_loader.get_elements();
    n_vertex = vertices.size() / 3;
  } else if (ext == "tet") {
    input >> n_vertex >> n_tet;
    vertices.resize(3 * n_vertex);
    indices.resize(n_tet << 2);

    for (int i = 0; i < n_vertex; ++i)
      input >> vertices[3 * i] >> vertices[3 * i + 1] >> vertices[3 * i + 2];

    for (int i = 0; i < n_tet; ++i) {
      input >> temp >> indices[(i << 2)] >> indices[(i << 2) + 1] >>
          indices[(i << 2) + 2] >> indices[(i << 2) + 3];
      assert(temp == 4);
    }
  } else if (ext == "vtk") {
    for (int i = 0; i < 4; ++i) std::getline(input, s);  // skip first 4 lines

    input >> s >> n_vertex >> s;
    vertices.resize(3 * n_vertex);
    for (int i = 0; i < n_vertex; ++i)
      input >> vertices[3 * i] >> vertices[3 * i + 1] >> vertices[3 * i + 2];

    input >> s >> n_tet >> s;
    indices.resize(n_tet << 2);
    for (int i = 0; i < n_tet; ++i) {
      // A single left shift multiplies a binary number by 2:
      input >> temp >> indices[(i << 2)] >> indices[(i << 2) + 1] >>
          indices[(i << 2) + 2] >> indices[(i << 2) + 3];
      assert(temp == 4);
      for (uint j = 0; j < 4; ++j) --indices[(i << 2) + j];
    }
  } else {
    input.close();
    return false;
  }
  input.close();
  std::cout << "loaded tet_mesh #v: " << n_vertex
            << ", #t: " << indices.size() / 4 << std::endl;
  get_tet_euler(vertices, indices);

  // normalize vertices between [0,1000]^3
  // float xmin, ymin, zmin, xmax, ymax, zmax;
  get_bbox(vertices, params.xmin_orig, params.ymin_orig, params.zmin_orig,
           params.xmax_orig, params.ymax_orig, params.zmax_orig,
           params.bbox_diag_l_prev);
  get_bbox(vertices, params.xmin, params.ymin, params.zmin, params.xmax,
           params.ymax, params.zmax, params.bbox_diag_l);
  if (normalize) {
    float maxside =
        std::max(std::max(params.xmax - params.xmin, params.ymax - params.ymin),
                 params.zmax - params.zmin);
    params.scale_maxside_orig = maxside;
    printf("tet has maxside: %f\n", maxside);
    // #pragma omp parallel for
    for (int i = 0; i < n_vertex; i++) {
      vertices[3 * i] =
          params.scale_max * (vertices[3 * i] - params.xmin) / maxside;
      vertices[3 * i + 1] =
          params.scale_max * (vertices[3 * i + 1] - params.ymin) / maxside;
      vertices[3 * i + 2] =
          params.scale_max * (vertices[3 * i + 2] - params.zmin) / maxside;
    }
    get_bbox(vertices, params.xmin, params.ymin, params.zmin, params.xmax,
             params.ymax, params.zmax, params.bbox_diag_l);
    std::cerr << "bbox [" << params.xmin << ":" << params.xmax << "], ["
              << params.ymin << ":" << params.ymax << "], [" << params.zmin
              << ":" << params.zmax << "], bbox_diag_l: " << params.bbox_diag_l
              << std::endl;
  }

  // update 8 bbox points
  Vector3 pmin(params.xmin, params.ymin, params.zmin);
  Vector3 pmax(params.xmax, params.ymax, params.zmax);
  for (int i = 0; i < 8; i++) {
    std::array<double, 3> p;
    std::bitset<sizeof(int) * 8> flag(i);
    for (int j = 0; j < 3; j++) {
      if (flag.test(j))
        p[j] = pmax[j];
      else
        p[j] = pmin[j];
    }
    params.bb_points.push_back(p[0]);
    params.bb_points.push_back(p[1]);
    params.bb_points.push_back(p[2]);
  }

  return true;
}

void normalize_mesh_given_scale(const float scale_given, const float xmax,
                                const float xmin, const float ymax,
                                const float ymin, const float zmax,
                                const float zmin,
                                std::vector<std::array<float, 3>>& vertices) {
  printf("Rescaling & Centering, scale_given: %d\n", scale_given);
  float maxside = std::max(std::max(xmax - xmin, ymax - ymin), zmax - zmin);
  // #pragma omp parallel for
  for (int i = 0; i < vertices.size(); i++) {
    vertices[i][0] = scale_given * (vertices[i][0] - xmin) / maxside;
    vertices[i][1] = scale_given * (vertices[i][1] - ymin) / maxside;
    vertices[i][2] = scale_given * (vertices[i][2] - zmin) / maxside;
  }
}

/**
 * 1. We store the number of shared adjacent cells for furture calculation of
 * Euler.
 * 2. We also stores unique id (int2) for each face, defined by (fid, -1)
 *
 * Vertex Eulers are stored as indices
 * Euler(v) = 1. / (#adjacent cells)
 *
 * Edge Eulers are stored as a diagonal adjacency matrix
 *
 * Face Eulers (originally from tets) are always 1/2 (adjacent to 2 cells)
 * except boundary
 *
 * tet_vs2sf_fids: tet vertices -> sf_mesh adjacent facets, if exists
 * n_sf_facets: number of surface facets (on boundary of tets)
 */
void load_tet_adj_info(const std::map<int, std::set<int>>& v2tets,
                       const std::vector<int>& tet_indices,
                       const std::map<int, std::set<int>>& tet_vs2sf_fids,
                       const int n_sf_facets, std::vector<int>& v_adjs,
                       std::vector<int>& e_adjs, std::vector<int>& f_adjs,
                       std::vector<int>& f_ids) {
  if (tet_vs2sf_fids.empty() || v2tets.empty()) {
    printf("tet_vs2sf_fids size: %ld, v2tets size %ld \n",
           tet_vs2sf_fids.size(), v2tets.size());
    assert(false);
  }
  v_adjs.clear();
  e_adjs.clear();
  f_adjs.clear();
  f_ids.clear();
  int nb_p = v2tets.size();
  int nb_v = tet_indices.size() / 4;

  // load #adjacent cells of vertices
  v_adjs.resize(nb_p, UNK_INT);
  for (uint v = 0; v < nb_p; v++) {
    int nb_vcells = v2tets.at(v).size();
    v_adjs[v] = nb_vcells;
  }

  // load #adjacent cells of edges
  e_adjs.resize(nb_p * (1 + nb_p) / 2 + 1, UNK_INT);
  std::set<int> neigh_cells;
  for (uint t = 0; t < nb_v; t++) {
    // 6 edges of a tet
    for (uint lv = 0; lv < 4; lv++) {
      for (uint lv_next = lv + 1; lv_next < 4; lv_next++) {
        int v1 = tet_indices[t * 4 + lv];
        int v2 = tet_indices[t * 4 + lv_next];
        int eid = get_edge_idx_copy(v1, v2, nb_p);
        if (e_adjs[eid] != UNK_INT) continue;  // already stored
        // check common cells of this edge
        neigh_cells.clear();
        set_intersection<int>(v2tets.at(v1), v2tets.at(v2), neigh_cells);
        assert(!neigh_cells.empty());
        e_adjs[eid] = neigh_cells.size();
      }
    }
  }

  // load #adjacent faces, each tet stores 4 faces
  //
  // each face has a unique face id (even shared by two tets),
  // we store int2 as face id, {{fid, -1}} represents face from tets.
  //
  // If the fid is on boundary (sf_mesh), we store the sf_mesh facet id,
  // otherwise we increment fid
  std::array<int, 3> f_vids;
  std::map<std::array<int, 2>, int> f_2tids_visited;
  int fid = n_sf_facets;
  for (uint t = 0; t < nb_v; t++) {
    FOR(i, 4) {
      // matching halplanes inited in ConvexCell::ConvexCell()
      f_vids = {tet_indices[t * 4 + tet_faces_lvid_host[i][0]],
                tet_indices[t * 4 + tet_faces_lvid_host[i][1]],
                tet_indices[t * 4 + tet_faces_lvid_host[i][2]]};
      std::sort(f_vids.begin(), f_vids.end());
      neigh_cells.clear();
      set_intersection<int>(v2tets.at(f_vids[0]), v2tets.at(f_vids[1]),
                            v2tets.at(f_vids[2]), neigh_cells);
      assert(!neigh_cells.empty() && neigh_cells.size() < 3);
      f_adjs.push_back(neigh_cells.size());
      // store face id, each id is unique
      //
      // boundary face of tet, must be on sf_mesh
      if (neigh_cells.size() == 1) {
        std::set<int> neigh_sf_fids;
        set_intersection<int>(tet_vs2sf_fids.at(f_vids[0]),
                              tet_vs2sf_fids.at(f_vids[1]),
                              tet_vs2sf_fids.at(f_vids[2]), neigh_sf_fids);
        assert(neigh_sf_fids.size() == 1);
        f_ids.push_back(*neigh_sf_fids.begin());
        continue;
      }
      // each face is shared by {tid1, tid2}
      // fetch the stored fid if exists
      std::array<int, 2> f_2tid = {
          {*neigh_cells.begin(), *neigh_cells.rbegin()}};
      std::sort(f_2tid.begin(), f_2tid.end());
      if (f_2tids_visited.find(f_2tid) != f_2tids_visited.end())
        f_ids.push_back(f_2tids_visited.at(f_2tid));
      else {
        f_2tids_visited[f_2tid] = fid;
        f_ids.push_back(fid++);
      }
    }
  }
  assert(f_adjs.size() == tet_indices.size());
  assert(f_ids.size() == tet_indices.size());  // each tet has 4 faces

  printf("loaded #adjacent cells for v: %d, and e: %d, and f %d\n",
         v_adjs.size(), e_adjs.size(), f_adjs.size());
}

void load_spheres_to_sites_given(std::vector<MedialSphere>& all_medial_spheres,
                                 bool& site_is_transposed,
                                 std::vector<float>& site,
                                 std::vector<float>& site_weights,
                                 int& n_site) {
  std::vector<std::array<float, 4>> spheres;
  // spheres.push_back({{-2.251256, -1.251256, 1.543551, 1.543575}});
  // spheres.push_back({{0.000000, -3.000000, 1.684660, 1.684660}});
  // spheres.push_back({{0.000000, 3.000000, 2.552344, 2.552344}});

  // spheres.push_back({{-2.5, -1.5, 1.5, 1.0}});
  // spheres.push_back({{0.0, -3.0, 1.5, 1.0}});
  // spheres.push_back({{0.0, 3.0, 2.5, 1.0}});

  // spheres.push_back({{-2.251262, -1.251262, 1.543560, 1.543562}});
  // spheres.push_back({{0.000000, -3.000000, 1.684660, 1.684660}});
  // spheres.push_back({{2.251262, -1.251262, 1.543560, 1.643562}});

  spheres.push_back({{1.944161, -12.897811, -2.092603, 0.528340}});
  spheres.push_back({{-3.450441, 13.444244, -0.001805, 0.891773}});

  n_site = spheres.size();
  int dim = 4;
  std::cout << "n_site: " << n_site << ", dim: " << dim << std::endl;

  all_medial_spheres.clear();
  for (int i = 0; i < n_site; i++) {
    MedialSphere msphere(all_medial_spheres.size(),
                         Vector3(spheres[i][0], spheres[i][1], spheres[i][2]),
                         spheres[i][3], SphereType::T_UNK);
    all_medial_spheres.push_back(msphere);
  }

  site_is_transposed = true;
  site.resize(n_site * 3);
  site_weights.resize(n_site);
  for (int i = 0; i < n_site; ++i) {
    const auto& msphere = all_medial_spheres.at(i);
    if (site_is_transposed) {
      site[i] = msphere.center[0];
      site[i + n_site] = msphere.center[1];
      site[i + (n_site << 1)] = msphere.center[2];
    } else {
      site[3 * i] = msphere.center[0];
      site[3 * i + 1] = msphere.center[1];
      site[3 * i + 2] = msphere.center[2];
    }
    // add weight (sq_radius) info
    if (dim == 4) {
      site_weights[i] = std::pow(float(msphere.radius), 2);
    } else
      site_weights[i] = 1;
    // printf("site %d has sphere: (%lf %lf %lf %lf) \n", i, msphere.center[0],
    //        msphere.center[1], msphere.center[2], msphere.radius);
  }
}

void load_spheres_from_file(const char* filename,
                            std::vector<MedialSphere>& all_medial_spheres,
                            bool is_load_type, bool is_load_deleted) {
  std::ifstream file(filename);
  int dim, n_site;
  int type, is_deleted;
  file >> dim >> n_site;
  assert(dim == 3 || dim == 4);
  std::cout << "n_site: " << n_site << ", dim: " << dim << std::endl;

  all_medial_spheres.clear();
  for (int i = 0; i < n_site; i++) {
    bool is_store_msphere = true;
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         INIT_RADIUS, SphereType::T_UNK);
    file >> msphere.center[0] >> msphere.center[1] >> msphere.center[2];
    if (dim == 4)
      file >> msphere.radius;
    else
      msphere.radius = 1;
    if (is_load_type) {
      file >> type;
      msphere.type = SphereType(type);
      // do not load
      if (msphere.type == SphereType::T_UNK) is_store_msphere = false;
    }
    if (is_load_deleted) {
      file >> is_deleted;
      msphere.is_deleted = false;
      if (is_deleted) msphere.is_deleted = true;
    }
    if (is_store_msphere) {
      all_medial_spheres.push_back(msphere);
      // printf("site %d has sphere: (%lf %lf %lf %lf), type: %d\n", i,
      //        msphere.center[0], msphere.center[1], msphere.center[2],
      //        msphere.radius, msphere.type);
    }
  }
  file.close();
}

void save_spheres_file(const std::vector<MedialSphere>& all_medial_spheres,
                       const std::string filename, bool is_save_type,
                       bool is_load_deleted, int num_rpd_itr) {
  std::string folder_name = "../out/" + filename + "/sph/";
  create_dir(folder_name);
  std::string sphere_path =
      folder_name + "sph_" + filename + "itr" +
      (num_rpd_itr == -1 ? get_timestamp() : std::to_string(num_rpd_itr)) +
      ".sph";

  int n_site = all_medial_spheres.size();
  std::fstream file;
  file.open(sphere_path, std::ios_base::out);
  file << 4 << " " << n_site << std::endl;
  for (int i = 0; i < n_site; i++) {
    const auto& msphere = all_medial_spheres.at(i);
    file << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << msphere.center[0] << " " << msphere.center[1] << " "
         << msphere.center[2] << " " << msphere.radius;
    if (is_save_type) file << " " << msphere.type;
    if (is_load_deleted) file << " " << msphere.is_deleted;
    file << std::endl;
  }
  file.close();
  printf("saved .sph file %s\n", sphere_path.c_str());
}

void load_v2tets(const std::vector<float>& vertices,
                 const std::vector<int>& indices,
                 std::map<int, std::set<int>>& v2tets) {
  v2tets.clear();
  int n_tets = indices.size() / 4;
  for (int t = 0; t < n_tets; t++) {
    for (uint lv = 0; lv < 4; lv++) {
      v2tets[indices[t * 4 + lv]].insert(t);
    }
  }
  // sanity
  if (v2tets.size() != vertices.size() / 3) {
    std::cerr << "ERROR: vertices size / 3: " << vertices.size() / 3
              << " not equal to v2tets size: " << v2tets.size() << std::endl;
    std::cout << "indices size / 4: " << indices.size() / 4 << std::endl;
    exit(1);
  }
};

void load_surface_vertices(const std::vector<float>& vertices,
                           const std::vector<int>& indices,
                           std::vector<float>& surf_vertices) {
  uint dim = 4;
  assert(indices.size() % dim == 0);

  // facet -> count, surface facet only count 1, otherwise 2
  std::map<std::array<int, 3>, int> map_fcnt;  // facet count
  for (uint i = 0; i < indices.size() / dim; i++) {
    uint index = i * dim;
    for (uint j = 0; j < dim; j++) {
      std::array<int, 3> facet = {{indices.at(index + j),
                                   indices.at(index + (j + 1) % dim),
                                   indices.at(index + (j + 2) % dim)}};

      std::sort(facet.begin(), facet.end());
      if (map_fcnt.find(facet) == map_fcnt.end()) {
        map_fcnt[facet] = 0;
      }
      map_fcnt[facet]++;
    }
  }

  std::set<int> surf_vset;
  for (auto& pair : map_fcnt) {
    if (pair.second != 1) continue;  // not on surface
    auto& facet = pair.first;
    for (int i = 0; i < 3; i++) surf_vset.insert(facet[i]);
  }
  printf("[Surface] found %lu surface vertices\n", surf_vset.size());

  surf_vertices.clear();
  for (uint vid : surf_vset) {
    surf_vertices.push_back(vertices[vid]);
  }
}

bool load_surface_mesh(const std::string& path, GEO::Mesh& input) {
  std::cout << "Loading mesh at" << path << std::endl;
  input.clear(false, false);
  const bool ok = GEO::mesh_load(path, input);
  std::cout << "ok: " << ok << std::endl;
  return ok;
}

/**
 * Here we save the vertex's old id in tet_vertices as attributes
 */
bool load_surface_mesh_geogram(const std::string& path, GEO::Mesh& input) {
  std::cout << "Loading mesh at" << path << std::endl;
  if (get_file_ext(path) != "geogram") {
    printf("Please use mesh format as .geogram");
    return false;
  }
  input.clear(false, false);
  GEO::InputGeoFile geo_file(path);
  GEO::MeshIOFlags flags;
  flags.set_element(GEO::MeshElementsFlags::MESH_ALL_ELEMENTS);
  flags.set_attributes(GEO::MeshAttributesFlags::MESH_ALL_ATTRIBUTES);
  const bool ok = GEO::mesh_load(geo_file, input, flags);
  if (!ok) return false;
  return ok;
}

/**
 * sf2tet_vs_mapping: mapping to tet_vertices indices
 */
void load_sf_mesh_from_internal(const std::vector<std::array<float, 3>>& points,
                                const std::vector<std::array<int, 3>>& faces,
                                const std::vector<int>& sf2tet_vs_mapping,
                                GEO::Mesh& input) {
  std::cout << "Loading mesh from internal data ..." << std::endl;
  input.clear(false, false);
  GEO::Attribute<int> tet_vid_attr(input.vertices.attributes(), "tet_vid");

  // Setup vertices
  input.vertices.create_vertices(points.size());
  for (uint i = 0; i < input.vertices.nb(); ++i) {
    GEO::vec3& p = input.vertices.point(i);
    p[0] = points[i][0];
    p[1] = points[i][1];
    p[2] = points[i][2];
    if (!sf2tet_vs_mapping.empty()) tet_vid_attr[i] = sf2tet_vs_mapping[i];
  }

  // Setup faces
  input.facets.create_triangles(faces.size());
  for (uint c = 0; c < input.facets.nb(); ++c) {
    for (uint lv = 0; lv < 3; ++lv) {
      input.facets.set_vertex(c, lv, faces[c][lv]);
    }
  }

  // Setup edges
  std::vector<aint2> edges;
  for (uint i = 0; i < faces.size(); i++) {
    const auto& f = faces[i];
    for (uint j = 0; j < 3; j++) {
      aint2 e = {{f[j], f[(j + 1) % 3]}};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      edges.push_back(e);
    }
  }
  vector_unique(edges);
  input.edges.create_edges(edges.size());
  for (uint e = 0; e < edges.size(); e++) {
    for (uint lv = 0; lv < 2; ++lv) {
      input.edges.set_vertex(e, lv, edges[e][lv]);
    }
  }

  GEO::mesh_reorder(input, GEO::MESH_ORDER_MORTON);
  // we did not setup the adjacent info till now
  input.facets.connect();
}

// void write_convex_cells(std::vector<float3>& voro_points,
//                         std::vector<std::array<int, 4>>& voro_tets,
//                         std::vector<int>& voro_tets_sites) {
//   printf("writing cells to file ...");
//   std::ofstream file("cells.tet");
//   file << voro_points.size() << " " << voro_tets_sites.size() << std::endl;

//   for (int i = 0; i < voro_points.size(); i++)
//     file << voro_points[i].x << " " << voro_points[i].y << " "
//          << voro_points[i].z << std::endl;

//   // tet + site id
//   for (int i = 0; i < voro_tets.size(); i++)
//     file << "5 " << voro_tets[i][0] << " " << voro_tets[i][1] << " "
//          << voro_tets[i][2] << " " << voro_tets[i][3] << " "
//          << voro_tets_sites[i] << std::endl;

//   file.close();
// }

/**
 * sf2tet_vs_mapping: mapping from surface vertex id to tet_vertices id
 */
void get_surface_from_tet(const std::vector<float>& tet_vertices,
                          const std::vector<int>& tet_indices,
                          std::vector<std::array<float, 3>>& surf_vertices,
                          std::vector<std::array<int, 3>>& surf_faces,
                          std::vector<int>& sf2tet_vs_mapping) {
  // Orignal code:
  // https://github.com/wildmeshing/fTetWild/blob/master/src/MeshImprovement.cpp#L2611
  //
  // find all faces
  assert(tet_indices.size() % 4 == 0);
  std::vector<std::array<int, 5>> faces;
  for (uint i = 0; i < tet_indices.size() / 4; i++) {
    uint idx = i * 4;
    for (uint j = 0; j < 4; j++) {
      std::array<int, 3> f = {{
          tet_indices[idx + (j + 1) % 4],
          tet_indices[idx + (j + 2) % 4],
          tet_indices[idx + (j + 3) % 4],
      }};
      std::sort(f.begin(), f.end());
      faces.push_back({{f[0], f[1], f[2], int(i), int(j)}});
    }
  }
  std::sort(faces.begin(), faces.end(),
            [](const std::array<int, 5>& a, const std::array<int, 5>& b) {
              return std::make_tuple(a[0], a[1], a[2]) <
                     std::make_tuple(b[0], b[1], b[2]);
            });
  if (faces.empty()) return;

  // find boundary faces
  std::vector<std::array<int, 3>> b_faces;
  bool is_boundary = true;
  for (uint i = 0; i < faces.size() - 1; i++) {
    if (std::make_tuple(faces[i][0], faces[i][1], faces[i][2]) ==
        std::make_tuple(faces[i + 1][0], faces[i + 1][1], faces[i + 1][2])) {
      is_boundary = false;
    } else {
      if (is_boundary) {
        b_faces.push_back({{faces[i][0], faces[i][1], faces[i][2]}});
        // check inv
        std::vector<Vector3> vs_four;
        uint v1_idx = tet_indices[faces[i][3] * 4 + faces[i][4]];  // i*4+j
        vs_four.push_back(Vector3(tet_vertices[v1_idx * 3],
                                  tet_vertices[v1_idx * 3 + 1],
                                  tet_vertices[v1_idx * 3 + 2]));
        for (uint j = 0; j < 3; j++) {
          vs_four.push_back(Vector3(tet_vertices[faces[i][j] * 3],
                                    tet_vertices[faces[i][j] * 3 + 1],
                                    tet_vertices[faces[i][j] * 3 + 2]));
        }
        bool is_inv =
            is_inverted(vs_four[0], vs_four[1], vs_four[2], vs_four[3]);
        if (!is_inv) std::swap(b_faces.back()[1], b_faces.back()[2]);
      }
      is_boundary = true;
    }
  }
  if (is_boundary) {
    b_faces.push_back({{faces.back()[0], faces.back()[1], faces.back()[2]}});
    // check inv
    std::vector<Vector3> vs_four;
    int v1_idx = tet_indices[faces.back()[3] * 4 + faces.back()[4]];  // i*4+j
    vs_four.push_back(Vector3(tet_vertices[v1_idx * 3],
                              tet_vertices[v1_idx * 3 + 1],
                              tet_vertices[v1_idx * 3 + 2]));
    for (int j = 0; j < 3; j++) {
      vs_four.push_back(Vector3(tet_vertices[faces.back()[j] * 3],
                                tet_vertices[faces.back()[j] * 3 + 1],
                                tet_vertices[faces.back()[j] * 3 + 2]));
    }
    bool is_inv = is_inverted(vs_four[0], vs_four[1], vs_four[2], vs_four[3]);
    if (!is_inv) std::swap(b_faces.back()[1], b_faces.back()[2]);
  }

  // find boundary vertices
  std::vector<int> b_v_ids;
  for (uint i = 0; i < b_faces.size(); i++) {
    for (uint j = 0; j < 3; j++) {
      b_v_ids.push_back(b_faces[i][j]);
    }
  }
  vector_unique(b_v_ids);

  // output
  surf_vertices.resize(b_v_ids.size());
  surf_faces.resize(b_faces.size());
  sf2tet_vs_mapping.resize(b_v_ids.size());
  std::map<int, int> map_v_ids /*old -> new*/;
  for (uint i = 0; i < b_v_ids.size(); i++) {
    map_v_ids[b_v_ids[i]] = i;
    // update id mapping, new -> old
    sf2tet_vs_mapping[i] = b_v_ids[i];
    surf_vertices[i] = {{tet_vertices[b_v_ids[i] * 3],
                         tet_vertices[b_v_ids[i] * 3 + 1],
                         tet_vertices[b_v_ids[i] * 3 + 2]}};
  }
  for (uint i = 0; i < b_faces.size(); i++) {
    surf_faces[i] = {{map_v_ids[b_faces[i][0]], map_v_ids[b_faces[i][1]],
                      map_v_ids[b_faces[i][2]]}};
  }
}

bool save_sf_mesh(const std::string sf_path, const GEO::Mesh& sf_mesh) {
  bool ok = GEO::mesh_save(sf_mesh, sf_path);
  std::cout << "saving sf_mesh at " << sf_path << ", ok: " << ok << std::endl;
  return ok;
}

bool save_sf_mesh_with_extf(const std::string sf_path,
                            const GEO::Mesh& sf_mesh) {
  std::ofstream fout(sf_path);
  fout << sf_mesh.vertices.nb() << " " << sf_mesh.edges.nb() << " "
       << sf_mesh.facets.nb() << std::endl;

  // save vertices
  for (int i = 0; i < sf_mesh.vertices.nb(); i++) {
    Vector3 pos = sf_mesh.vertices.point(i);
    fout << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << pos[0] << " " << pos[1] << " " << pos[2] << " 1 2 0" << std::endl;
  }

  // save edges
  GEO::Attribute<int> attr_se(sf_mesh.edges.attributes(), "se");
  GEO::Attribute<int> attr_cce(sf_mesh.edges.attributes(), "cce");
  for (int e = 0; e < sf_mesh.edges.nb(); e++) {
    fout << "e " << sf_mesh.edges.vertex(e, 0) << " "
         << sf_mesh.edges.vertex(e, 1);
    if (attr_se[e] == 1) {
      // printf("input mesh edge %d is on se\n", e);
      fout << " " << 2;  // external feature
    } else if (attr_cce[e] == 1) {
      fout << " " << 1;  // concave edges
    } else
      fout << " " << 0;
    fout << std::endl;
  }

  // save faces
  for (int f = 0; f < sf_mesh.facets.nb(); f++) {
    fout << "f " << sf_mesh.facets.vertex(f, 0) << " "
         << sf_mesh.facets.vertex(f, 1) << " " << sf_mesh.facets.vertex(f, 2)
         << std::endl;
  }
  fout.close();
  printf("saved sf_mesh with extf at: %s \n", sf_path.c_str());
}

bool save_sf_mesh_scaled_01(const std::string sf_path_scaled,
                            const GEO::Mesh& sf_mesh, const Parameter& params) {
  // load vertices from sf_mesh
  std::vector<std::array<float, 3>> vertices(sf_mesh.vertices.nb());
  for (uint v = 0; v < sf_mesh.vertices.nb(); ++v) {
    Vector3 pos = sf_mesh.vertices.point(v);
    for (uint lv = 0; lv < 3; ++lv) {
      vertices[v][lv] = pos[lv];
    }
  }
  // load faces from sf_mesh
  std::vector<std::array<int, 3>> faces(sf_mesh.facets.nb());
  for (uint c = 0; c < sf_mesh.facets.nb(); ++c) {
    for (uint lv = 0; lv < 3; ++lv) {
      faces[c][lv] = sf_mesh.facets.vertex(c, lv);
    }
  }

  // re-scale the sf_mesh vertices back to [0,1]
  normalize_mesh_given_scale(1 / params.scale_max, params.xmax, params.xmin,
                             params.ymax, params.ymin, params.zmax, params.zmin,
                             vertices);

  // save scaled mesh
  std::vector<int> _;
  GEO::Mesh sf_mesh_scaled;
  load_sf_mesh_from_internal(vertices, faces, _, sf_mesh_scaled);
  save_sf_mesh(sf_path_scaled, sf_mesh_scaled);
}

bool save_sf_mesh_geogram(const std::string sf_path, GEO::Mesh& sf_mesh) {
  GEO::OutputGeoFile geo_file(sf_path);
  GEO::MeshIOFlags flags;
  flags.set_element(GEO::MeshElementsFlags::MESH_ALL_ELEMENTS);
  flags.set_attributes(GEO::MeshAttributesFlags::MESH_ALL_ATTRIBUTES);
  if (!GEO::mesh_save(sf_mesh, geo_file, flags)) {
    std::cout << "Unable to save file at " << sf_path << std::endl;
    return false;
  }
  std::cout << "Done saving file  " << sf_path << std::endl;
  return true;
}

bool is_slice_by_plane(const Vector3& bary, const Parameter& params) {
  Vector3 slice_pl_pos(params.xmax / 2.2, 0, 0);
  Vector3 slice_pl_normal(-1, 0, 0);
  // Vector3 slice_pl_pos(0, params.ymax / 2.2, 0);
  // Vector3 slice_pl_normal(0, -1, 0);
  // Vector3 slice_pl_pos(0, 0, params.zmax / 2.5);
  // Vector3 slice_pl_normal(0, 0, 1);
  double dist = GEO::dot(bary - slice_pl_pos, slice_pl_normal);
  if (dist >= 0) return true;  // clip and do not save
  return false;
}

void unnormalize_matfp(const Parameter& params, MedialMesh& mmesh_matfp) {
  assert(!mmesh_matfp.vertices->empty());
  const Vector3 p0 = mmesh_matfp.vertices->at(0).center;
  Vector3 min(p0[0], p0[1], p0[2]);
  Vector3 max(p0[0], p0[1], p0[2]);
  int size_max = 10;

  for (int i = 0; i < mmesh_matfp.vertices->size(); ++i) {
    const Vector3& p = mmesh_matfp.vertices->at(i).center;
    for (int j = 0; j < 3; j++) {
      min[j] = std::min(min[j], p[j]);
      max[j] = std::max(max[j], p[j]);
    }
  }

  // same as original shape
  // load_tet() will update params.scale_maxside_orig
  assert(params.scale_maxside_orig != 0.f);
  double size_orig = params.scale_maxside_orig / 10.f;
  double xcenter = (params.xmax_orig + params.xmin_orig) * 0.5;
  double ycenter = (params.ymax_orig + params.ymin_orig) * 0.5;
  double zcenter = (params.zmax_orig + params.zmin_orig) * 0.5;

  // un-normalize matfp to original
  for (int i = 0; i < mmesh_matfp.vertices->size(); ++i) {
    Vector3& p = mmesh_matfp.vertices->at(i).center;
    p[0] = (p[0] * size_orig) + xcenter;
    p[1] = (p[1] * size_orig) + ycenter;
    p[2] = (p[2] * size_orig) + zcenter;
  }
}

void renormalize_matfp(const Parameter& params, MedialMesh& mmesh_matfp) {
  float maxside = std::max(std::max(params.xmax_orig - params.xmin_orig,
                                    params.ymax_orig - params.ymin_orig),
                           params.zmax_orig - params.zmin_orig);
  for (int i = 0; i < mmesh_matfp.vertices->size(); ++i) {
    Vector3& p = mmesh_matfp.vertices->at(i).center;
    p[0] = params.scale_max * (p[0] - params.xmin_orig) / maxside;
    p[1] = params.scale_max * (p[1] - params.ymin_orig) / maxside;
    p[2] = params.scale_max * (p[2] - params.zmin_orig) / maxside;
  }
}

// this export is matching MATFP (contains flag_type and flag_delete)
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r flag_type flag_delete
 * e v1 v2 flag_delete
 * f v1 v2 v3
 */
void load_matfp(const std::string& ma_path,
                std::vector<MedialSphere>& all_medial_spheres,
                MedialMesh& mat) {
  std::ifstream ma_file(ma_path);
  int num_vs, num_es, num_fs;
  int type, is_deleted;
  char ch;
  ma_file >> num_vs >> num_es >> num_fs;
  printf("loading .ma with num_vs: %d, num_es: %d, num_fs: %d\n", num_vs,
         num_es, num_fs);

  // load spheres
  all_medial_spheres.clear();
  for (int i = 0; i < num_vs; i++) {
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         Vector3(0, 0, 0), SphereType::T_2);
    ma_file >> ch >> msphere.center[0] >> msphere.center[1] >>
        msphere.center[2];
    ma_file >> msphere.radius;
    ma_file >> type;
    if (type != SphereType::T_UNK)
      msphere.type = SphereType(type);  // else as T2 sphere
    ma_file >> is_deleted;
    msphere.is_deleted = false;
    if (is_deleted) msphere.is_deleted = true;
    all_medial_spheres.push_back(msphere);
    // printf("ma v %d has sphere: (%lf %lf %lf %lf), type: %d\n", i,
    //        msphere.center[0], msphere.center[1], msphere.center[2],
    //        msphere.radius, msphere.type);
  }

  // load mat
  mat.clear();
  mat.vertices = &all_medial_spheres;
  int e1, e2;
  for (int e = 0; e < num_es; e++) {
    ma_file >> ch >> e1 >> e2 >> is_deleted;
    if (is_deleted) continue;
    // printf("ma e %d has (%d,%d)\n", e, e1, e2);
    mat.create_edge(e1, e2);
  }

  int f1, f2, f3;
  for (int f = 0; f < num_fs; f++) {
    ma_file >> ch >> f1 >> f2 >> f3;
    // printf("ma f %d has (%d,%d,%d)\n", f, f1, f2, f3);
    aint3 fvs = {{f1, f2, f3}};
    mat.create_face(fvs);
  }

  ma_file.close();
}

// this export is matching MATFP clean format
// see MATFP commit fdd2d19dd73e10262fed885f414c338d0d8c6151
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r
 * e v1 v2
 * f v1 v2 v3
 */
void load_mat_clean(const std::string& ma_path,
                    std::vector<MedialSphere>& all_medial_spheres,
                    MedialMesh& mat) {
  std::ifstream ma_file(ma_path);
  int num_vs, num_es, num_fs;
  int type, is_deleted;
  char ch;
  ma_file >> num_vs >> num_es >> num_fs;
  printf("loading .ma with num_vs: %d, num_es: %d, num_fs: %d\n", num_vs,
         num_es, num_fs);

  // load spheres
  all_medial_spheres.clear();
  for (int i = 0; i < num_vs; i++) {
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         Vector3(0, 0, 0), SphereType::T_2);
    ma_file >> ch >> msphere.center[0] >> msphere.center[1] >>
        msphere.center[2];
    ma_file >> msphere.radius;
    all_medial_spheres.push_back(msphere);
    // printf("ma v %d has sphere: (%lf %lf %lf %lf), type: %d\n", i,
    //        msphere.center[0], msphere.center[1], msphere.center[2],
    //        msphere.radius, msphere.type);
  }

  // load mat
  mat.clear();
  mat.vertices = &all_medial_spheres;
  int e1, e2;
  for (int e = 0; e < num_es; e++) {
    ma_file >> ch >> e1 >> e2;
    // printf("ma e %d has (%d,%d)\n", e, e1, e2);
    mat.create_edge(e1, e2);
  }

  int f1, f2, f3;
  for (int f = 0; f < num_fs; f++) {
    ma_file >> ch >> f1 >> f2 >> f3;
    // printf("ma f %d has (%d,%d,%d)\n", f, f1, f2, f3);
    aint3 fvs = {{f1, f2, f3}};
    mat.create_face(fvs);
  }

  ma_file.close();
}

// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
/* # .ma format
 * # all flag_delete might not exist
 * numVertices numEdges numFaces
 * v x y z r flag_type flag_delete
 * e v1 v2 flag_type flag_delete
 * f v1 v2 v3 flag_delete
 */
void export_ma(const std::string& maname, const MedialMesh& mat) {
  std::string folder_name = "../out/" + maname + "/mat/";
  create_dir(folder_name);
  std::string ma_name_full =
      folder_name + "mat_" + maname + "_" + get_timestamp() + ".ma";

  std::ofstream fout(ma_name_full);
  fout << mat.vertices->size() << " " << mat.numEdges_no_dup << " "
       << mat.numFaces_no_dup << std::endl;

  // we might have redundant vertices
  // but we need those indices, so store it as delete flag = true
  for (int i = 0; i < mat.vertices->size(); i++) {
    Vector3 pos = mat.vertices->at(i).center;
    fout << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << pos[0] << " " << pos[1] << " " << pos[2] << " "
         << mat.vertices->at(i).radius << " " << int(mat.vertices->at(i).type);

    // store this for not showing redundant vertices
    if (mat.vertices->at(i).is_deleted) {
      fout << " 1";
    } else {
      fout << " 0";
    }
    fout << std::endl;
  }

  // do not save edge if deleted
  for (int i = 0; i < mat.edges.size(); i++) {
    const auto& mat_e = mat.edges[i];
    if (mat_e.is_deleted) continue;
    aint2 edge = {{mat_e.vertices_[0], mat_e.vertices_[1]}};
    std::sort(edge.begin(), edge.end());
    fout << "e " << edge[0] << " " << edge[1];
    // is this edge a feature edge
    if (mat_e.is_intf) {
      fout << " " << 1;  // internal
    } else if (mat_e.is_extf) {
      fout << " " << 2;  // external
    } else {
      fout << " " << 0;
    }
    fout << std::endl;
  }

  // do not save face if deleted
  for (int i = 0; i < mat.faces.size(); i++) {
    const auto& mat_f = mat.faces[i];
    if (mat_f.is_deleted) continue;
    fout << "f";
    for (uint v = 0; v < 3; v++) fout << " " << mat_f.vertices_[v];
    fout << std::endl;
  }
  fout.close();

  printf("saved mat at: %s \n", ma_name_full.c_str());
}

// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
//
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r
 * e v1 v2
 * f v1 v2 v3
 */
void export_ma_given(const std::string& maname,
                     const std::vector<Vector4>& mat_vertices,
                     const std::vector<aint2>& mat_edges,
                     const std::vector<std::array<int, 3>>& mat_faces,
                     bool is_use_given_name) {
  std::string ma_name_full = maname;
  if (!is_use_given_name) {
    std::string folder_name = "../out/" + maname + "/mat/";
    create_dir(folder_name);
    ma_name_full =
        folder_name + "mat_" + maname + "_" + get_timestamp() + ".ma";
  }

  std::ofstream fout;
  fout.open(ma_name_full, std::ofstream::out | std::ofstream::app);  //   append
  fout << mat_vertices.size() << " " << mat_edges.size() << " "
       << mat_faces.size() << std::endl;

  // save vertices
  for (int i = 0; i < mat_vertices.size(); i++) {
    const auto& mat_v = mat_vertices.at(i);
    fout << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << mat_v[0] << " " << mat_v[1] << " " << mat_v[2] << " " << mat_v[3];
    fout << std::endl;
  }

  //  save edges
  for (int i = 0; i < mat_edges.size(); i++) {
    const auto& mat_e = mat_edges[i];
    fout << "e " << mat_e[0] << " " << mat_e[1];
    fout << std::endl;
  }

  // save faces
  for (int i = 0; i < mat_faces.size(); i++) {
    const auto& mat_f = mat_faces[i];
    fout << "f";
    for (uint v = 0; v < 3; v++) fout << " " << mat_f[v];
    fout << std::endl;
  }
  fout.close();

  printf("saved mat at: %s \n", ma_name_full.c_str());
}

// helper function for export_ma_clean() and export_ma_ply()
void get_mat_clean(const MedialMesh& mat, std::vector<Vector4>& vertices,
                   std::vector<aint2>& edges,
                   std::vector<std::array<int, 3>>& faces) {
  std::map<int, int> map_vertices;  // mat vertex tag to new
  vertices.clear();
  edges.clear();
  faces.clear();

  auto get_vertex_mapped_id = [&](const int vid) {
    if (map_vertices.find(vid) == map_vertices.end()) {
      // add a new vertex
      map_vertices[vid] = map_vertices.size();
    }
    return map_vertices.at(vid);
  };

  // faces
  for (int f = 0; f < mat.faces.size(); f++) {
    const auto& face = mat.faces[f];
    if (face.is_deleted) continue;
    std::array<int, 3> one_f;
    for (uint j = 0; j < 3; j++) {
      int vid = get_vertex_mapped_id(face.vertices_[j]);
      one_f[j] = vid;
    }
    faces.push_back(one_f);
  }
  printf("faces: %d \n", faces.size());

  // edges
  for (int e = 0; e < mat.edges.size(); e++) {
    const auto& edge = mat.edges[e];
    if (edge.is_deleted) continue;
    // if (edge.faces_.empty()) continue;
    int vid1 = get_vertex_mapped_id(edge.vertices_[0]);
    int vid2 = get_vertex_mapped_id(edge.vertices_[1]);
    edges.push_back({{vid1, vid2}});
  }
  printf("edges: %d \n", edges.size());

  // vertices
  // save from map_vertices, to avoid (0,0,0,0) in .ma file
  vertices.resize(map_vertices.size());
  for (const auto& v_pair : map_vertices) {
    int old_vid = v_pair.first;
    int new_vid = v_pair.second;
    const auto& mat_v = mat.vertices->at(old_vid);
    vertices[new_vid] = Vector4(mat_v.center[0], mat_v.center[1],
                                mat_v.center[2], mat_v.radius);
  }
  // for (int v = 0; v < mat.vertices->size(); v++) {
  //   const auto& mat_v = mat.vertices->at(v);
  //   if (mat_v.is_deleted) continue;
  //   if (mat_v.edges_.empty() && mat_v.faces_.empty()) continue;  //
  //   isolated？ int vid = get_vertex_mapped_id(mat_v.id); if (vid ==
  //   map_vertices.size()) {
  //     // vertices.push_back(Vector4(vertex.center[0], vertex.center[1],
  //     //                            vertex.center[2], vertex.radius));
  //     printf("mat sphere: %d not connect to any edge/face, wrong?\n", v);
  //     assert(false);
  //   } else {
  //     vertices[vid] = Vector4(mat_v.center[0], mat_v.center[1],
  //     mat_v.center[2],
  //                             mat_v.radius);
  //   }
  // }
  printf("vertcies: %d \n", vertices.size());
}

// remove deleted edges/faces
// should save the same result as export_ma_ply but with different format
//
// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r
 * e v1 v2
 * f v1 v2 v3
 */
void export_ma_clean(const std::string& maname, const MedialMesh& mat,
                     bool is_use_given_name) {
  std::string ma_name_full = maname;
  if (!is_use_given_name) {
    std::string folder_name = "../out/" + maname + "/mat/";
    create_dir(folder_name);
    ma_name_full =
        folder_name + "mat_" + maname + "_" + get_timestamp() + ".ma";
  }
  printf("start saving mat .m file: %s \n", ma_name_full.c_str());

  std::vector<Vector4> vertices;
  std::vector<aint2> edges;
  std::vector<std::array<int, 3>> faces;
  get_mat_clean(mat, vertices, edges, faces);

  // export
  export_ma_given(maname, vertices, edges, faces, true);
}

// [no use]
// some vertices/faces are deleted
void write_ma_ply(const std::string& maname, const MedialMesh& mat,
                  bool is_use_given_name) {
  std::string ma_name_full = maname;
  if (!is_use_given_name) {
    std::string folder_name = "../out/" + maname + "/mat/";
    create_dir(folder_name);
    ma_name_full =
        folder_name + "mat_" + maname + "_" + get_timestamp() + ".ply";
  }

  printf("start saving mat .ply file: %s \n", ma_name_full.c_str());

  std::map<int, int> map_vertices;  // mat vertex tag to new ply
  int num_vertices = 0;
  int num_edges = 0;
  int num_faces = 0;
  std::vector<Vector3> vertices;
  std::vector<aint2> edges;
  std::vector<std::array<int, 3>> faces;

  auto get_vertex_mapped_id = [&](const int vid) {
    if (map_vertices.find(vid) == map_vertices.end()) {
      map_vertices[vid] = num_vertices;
      num_vertices++;
    }
    return map_vertices.at(vid);
  };

  for (int f = 0; f < mat.faces.size(); f++) {
    const auto& face = mat.faces[f];
    if (face.is_deleted) continue;
    num_faces++;
    std::array<int, 3> one_f;
    for (uint j = 0; j < 3; j++) {
      int vid = get_vertex_mapped_id(face.vertices_[j]);
      one_f[j] = vid;
    }
    faces.push_back(one_f);
  }
  printf("faces: %d \n", num_faces);

  for (int e = 0; e < mat.edges.size(); e++) {
    const auto& edge = mat.edges[e];
    if (edge.is_deleted) continue;
    num_edges++;
    int vid1 = get_vertex_mapped_id(edge.vertices_[0]);
    int vid2 = get_vertex_mapped_id(edge.vertices_[1]);
    edges.push_back({{vid1, vid2}});
  }
  printf("edges: %d \n", num_edges);

  vertices.resize(map_vertices.size());
  for (int v = 0; v < mat.vertices->size(); v++) {
    auto& vertex = mat.vertices->at(v);
    if (vertex.is_deleted) continue;
    int vid = get_vertex_mapped_id(vertex.id);
    if (vid == vertices.size()) {
      vertices.push_back(vertex.center);
    } else {
      vertices[vid] = vertex.center;
    }
  }
  printf("vertcies: %d \n", num_vertices);

  // TODO: use geogram to save? libigl is not included
  //
  // MatrixXs V(num_vertices, 3);
  // Eigen::MatrixXi E(num_edges, 2);
  // Eigen::MatrixXi F(num_faces, 3);
  // for (int v = 0; v < num_vertices; v++) {
  //   V(v, 0) = vertices[v][0];
  //   V(v, 1) = vertices[v][1];
  //   V(v, 2) = vertices[v][2];
  // }
  // for (int e = 0; e < num_edges; e++) {
  //   E(e, 0) = edges[e][0];
  //   E(e, 1) = edges[e][1];
  // }
  // for (int f = 0; f < num_faces; f++) {
  //   F(f, 0) = faces[f][0];
  //   F(f, 1) = faces[f][1];
  //   F(f, 2) = faces[f][2];
  // }
  // igl::writePLY(ma_name_full, V, F, E);
}