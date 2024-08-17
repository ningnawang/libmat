
#include "main_gui_cxx.h"

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>

#include "../src/IO/IO_CXX/io.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// MainGuiWindow functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
MainGuiWindow* MainGuiWindow::instance_ = nullptr;
MainGuiWindow::MainGuiWindow() {
  if (instance_ != nullptr) {
    printf("ERROR: GuiWindow instance is not nullptr!!");
    exit(1);
  }
  instance_ = this;
}

MainGuiWindow::~MainGuiWindow() {
  // we delete here since we allocate memory
  // (called new) for each of these variables
  delete sf_mesh;
  delete tet_mesh;
  delete params;
  delete all_medial_spheres;
  delete mmesh;
  // delete tet_vertices;
  // delete tet_indices;
  // delete sf_vertices;
  // delete sf_faces;
  this->instance_ = nullptr;
}
void MainGuiWindow::set_file_name(const std::string& _file_name) {
  file_name = _file_name;
}
void MainGuiWindow::set_tet_mesh(TetMesh& _tet_mesh) { tet_mesh = &_tet_mesh; }
void MainGuiWindow::set_sf_mesh(SurfaceMesh& _sf_mesh) { sf_mesh = &_sf_mesh; }
void MainGuiWindow::set_params(Parameter& _params) { params = &_params; }
void MainGuiWindow::set_all_medial_spheres(
    std::vector<MedialSphere>& _all_medial_spheres) {
  all_medial_spheres = &_all_medial_spheres;
}
void MainGuiWindow::set_medial_mesh(MedialMesh& _mmesh) { mmesh = &_mmesh; }
// void MainGuiWindow::set_tet_mesh(std::vector<float>& _tet_vertices,
//                                  std::vector<int>& _tet_indices) {
//   tet_vertices = &_tet_vertices;
//   tet_indices = &_tet_indices;
// }
// void MainGuiWindow::set_sf_mesh_internal(
//     std::vector<std::array<float, 3>>& _sf_vertices,
//     std::vector<std::array<int, 3>>& _sf_faces) {
//   sf_vertices = &_sf_vertices;
//   sf_faces = &_sf_faces;
// }

void MainGuiWindow::show() {
  // a few camera options
  polyscope::view::upDir = polyscope::UpDir::ZUp;

  // Initialize Polyscope
  polyscope::init();

  // Add the callback
  polyscope::state::userCallback = MainGuiWindow::callbacks;

  // For debug
  // auto sf_mesh_internal = polyscope::registerSurfaceMesh(
  //     "Surface mesh internal", *(this->sf_vertices), *(this->sf_faces));
  // Surface
  // instance_->show_surface_mesh(*(instance_->sf_mesh), -1
  // /*given_sf_face_id*/);

  instance_->show_tet_mesh(*(instance_->tet_mesh), false /*is_shown_meta*/);
  instance_->show_medial_mesh(*this->mmesh);
  instance_->show_medial_edges(*this->mmesh);

  // Show the GUI
  polyscope::show();
}

// Your callback functions
void MainGuiWindow::callbacks() {
  ImGui::PushItemWidth(100);
  if (ImGui::Button("Scale to MATTopo")) {
    unnormalize_matfp(*instance_->params, *instance_->mmesh);
    renormalize_matfp(*instance_->params, *instance_->mmesh);
    printf("Done unnormalize_matfp & renormalize_matfp.\n");
    std::string ma_scale_name = "../out/mat/scaled_" + instance_->file_name +
                                "_" + get_timestamp() + ".ma";
    export_ma_clean(ma_scale_name, *(instance_->mmesh),
                    true /*is_use_given_name*/);
    instance_->show_medial_mesh(*instance_->mmesh);
    instance_->show_medial_edges(*instance_->mmesh);
  }

  if (ImGui::Button("Show MAT")) {
    instance_->show_medial_mesh(*instance_->mmesh);
    instance_->show_medial_edges(*instance_->mmesh);
  }

  if (ImGui::Button("Show TetMesh")) {
    instance_->show_tet_mesh(*(instance_->tet_mesh), false /*is_shown_meta*/);
  }

  // Surface
  // if (ImGui::Button("show surface mesh")) {
  //   instance_->show_surface_mesh(*(instance_->sf_mesh),
  //                                -1 /*given_sf_face_id*/);
  // }
  // ImGui::SameLine();
  if (ImGui::Button("show pins")) {
    instance_->show_pin_points(*(instance_->all_medial_spheres),
                               -1 /*given_sphere_id*/);
  }

  ImGui::InputInt("#face_id", &(instance_->given_sf_face_id));
  ImGui::SameLine();
  if (ImGui::Button("show")) {
    instance_->show_surface_mesh(*(instance_->sf_mesh),
                                 instance_->given_sf_face_id);
  }

  // Medial spheres
  if (ImGui::SmallButton("show all mat centers")) {
    instance_->show_all_medial_spheres(*(instance_->all_medial_spheres));
  }

  ImGui::InputInt("#sphere_id", &(instance_->given_sphere_id));
  ImGui::SameLine();
  if (ImGui::SmallButton("show/clear")) {
    instance_->show_one_sphere(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id,
                               true, /*is_show_radius*/
                               true  /*is_clear_all*/
    );
    instance_->show_pin_points(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id);
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("add sphere no R")) {
    instance_->show_one_sphere(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id,
                               false, /*is_show_radius*/
                               false  /*is_clear_all*/
    );
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("with R")) {
    instance_->show_one_sphere(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id,
                               true, /*is_show_radius*/
                               false /*is_clear_all*/
    );
  }

  //   ImGui::SameLine();
  ImGui::PopItemWidth();
}

void MainGuiWindow::show_tet_mesh(const TetMesh& tet_mesh, bool is_shown_meta) {
  std::vector<std::array<float, 3>> tet_vertices_new;
  std::vector<std::array<int, 4>> tet_indices_new;
  std::vector<int> tet_ids;
  for (int i = 0; i < tet_mesh.tet_vertices.size() / 3; i++) {
    tet_vertices_new.push_back(
        {{tet_mesh.tet_vertices[i * 3], tet_mesh.tet_vertices[i * 3 + 1],
          tet_mesh.tet_vertices[i * 3 + 2]}});
  }
  for (int i = 0; i < tet_mesh.tet_indices.size() / 4; i++) {
    tet_indices_new.push_back(
        {{tet_mesh.tet_indices[i * 4], tet_mesh.tet_indices[i * 4 + 1],
          tet_mesh.tet_indices[i * 4 + 2], tet_mesh.tet_indices[i * 4 + 3]}});
    tet_ids.push_back(i);
  }
  auto my_tet =
      polyscope::registerTetMesh("my tet", tet_vertices_new, tet_indices_new);
  my_tet->addCellScalarQuantity("tet_id", tet_ids);
  my_tet->setEnabled(true);

  if (!is_shown_meta) return;
  // show concave corners
  std::vector<int> is_vertex_cc_corner(tet_vertices_new.size(), false);
  for (const auto& cc_corner : tet_mesh.cc_corners) {
    is_vertex_cc_corner[cc_corner.tvid] = true;
  }
  my_tet->addVertexScalarQuantity("is_cc_corner", is_vertex_cc_corner);
  // show concave feature lines
  std::vector<Vector3> ce_edge_vs;
  std::vector<aint2> ce_edges;
  std::vector<int> ce_fl_ids;
  for (const auto& one_ce_line : tet_mesh.ce_lines) {
    for (const auto& one_ce_id : one_ce_line.fe_ids) {  // FeatureEdge::id
      const auto one_ce = tet_mesh.feature_edges.at(one_ce_id);
      ce_edge_vs.push_back(one_ce.t2vs_pos[0]);  // v0
      ce_edge_vs.push_back(one_ce.t2vs_pos[1]);  // v1
      ce_edges.push_back({{ce_edge_vs.size() - 2, ce_edge_vs.size() - 1}});
      ce_fl_ids.push_back(one_ce_line.id);  // FeatureLine::id
    }
  }
  auto my_ce_edges =
      polyscope::registerCurveNetwork("ce edges", ce_edge_vs, ce_edges);
  my_ce_edges->addEdgeScalarQuantity("fl_id", ce_fl_ids)->setEnabled(true);
  my_ce_edges->setEnabled(false);
}

void MainGuiWindow::show_surface_mesh(const GEO::Mesh& sf_mesh,
                                      const int& given_sf_face_id) {
  std::vector<std::array<double, 3>> sf_vertices;
  std::vector<std::array<int, 3>> sf_faces;  // vs size might be > 3
  std::vector<bool> is_selected(sf_mesh.facets.nb(), false);

  for (uint v = 0; v < sf_mesh.vertices.nb(); v++) {
    const Vector3& p = sf_mesh.vertices.point(v);
    sf_vertices.push_back({{p[0], p[1], p[2]}});
  }

  for (uint f = 0; f < sf_mesh.facets.nb(); f++) {
    if (given_sf_face_id == f) {
      is_selected[f] = true;
    }
    int f_nb_v = sf_mesh.facets.nb_vertices(f);
    assert(f_nb_v == 3);
    std::array<int, 3> one_f;
    for (uint lv = 0; lv < f_nb_v; lv++) {
      int v = sf_mesh.facets.vertex(f, lv);
      one_f[lv] = v;
    }
    sf_faces.push_back(one_f);
  }
  auto poly_sf_mesh =
      polyscope::registerSurfaceMesh("Surface mesh", sf_vertices, sf_faces);
  poly_sf_mesh->addFaceScalarQuantity("is_selected", is_selected)
      ->setEnabled(true);
}

void MainGuiWindow::show_pin_points(
    const std::vector<MedialSphere>& all_medial_spheres,
    const int given_sphere_id) {
  std::vector<Vector3> all_pins;
  std::vector<Vector3> all_pin_normals;

  for (uint i = 0; i < all_medial_spheres.size(); i++) {
    if (given_sphere_id != -1 && i != given_sphere_id) continue;
    const auto& msphere = all_medial_spheres.at(i);
    all_pins.push_back(msphere.ss.p);
    all_pin_normals.push_back(msphere.ss.p_normal);
  }

  auto all_pin = polyscope::registerPointCloud("All pins", all_pins);
  all_pin->addVectorQuantity("normal", all_pin_normals)->setEnabled(false);
}

void MainGuiWindow::show_all_medial_spheres(
    const std::vector<MedialSphere>& all_medial_spheres) {
  std::vector<Vector3> all_spheres;
  std::vector<int> sphere_ids;
  for (uint i = 0; i < all_medial_spheres.size(); i++) {
    const auto& msphere = all_medial_spheres.at(i);
    if (msphere.is_deleted) continue;
    all_spheres.push_back(msphere.center);
    sphere_ids.push_back(msphere.id);
  }
  auto all_pc = polyscope::registerPointCloud("All spheres", all_spheres);
  all_pc->setPointRadius(MainGuiWindow::point_radius_rel)->setEnabled(true);
  all_pc->addScalarQuantity("id", sphere_ids);
}

void MainGuiWindow::show_one_sphere(
    const std::vector<MedialSphere>& all_medial_spheres,
    const int given_sphere_id, const bool is_show_radius,
    const bool is_clear_all) {
  if (is_clear_all == false &&
      (given_sphere_id < 0 || given_sphere_id >= all_medial_spheres.size())) {
    polyscope::warning("msphere tag range is [0, " +
                       std::to_string(all_medial_spheres.size() - 1) + "]");
    return;
  }
  std::vector<Vector3> positions;
  std::vector<double> radii;
  const auto& msphere = all_medial_spheres.at(given_sphere_id);
  positions.push_back(msphere.center);
  radii.push_back(msphere.radius);
  msphere.print_info();

  std::string name_sphere = "Sphere ";
  // if is_clear_all, then shows nothing
  if (is_clear_all) {
    for (uint i = 1; i < 20; i++) {
      std::string new_sphere =
          name_sphere + "(with radius)" + std::to_string(i);
      polyscope::removePointCloud(new_sphere, false);
      new_sphere = name_sphere + "(no radius)" + std::to_string(i);
      polyscope::removePointCloud(new_sphere, false);
    }
  }

  if (is_show_radius)
    name_sphere += "(with radius)";
  else
    name_sphere += "(no radius)";
  for (uint i = 1; i < 20; i++) {
    std::string new_sphere = name_sphere + std::to_string(i);
    if (polyscope::hasPointCloud(new_sphere) == false) {
      auto single_msphere =
          polyscope::registerPointCloud(new_sphere, positions);
      if (is_show_radius) {
        single_msphere->addScalarQuantity("radius", radii);
        single_msphere->setPointRadiusQuantity("radius", false);
      }
      break;
    }
  }

  // Tangent planes
  std::vector<Vector3> tan_pl_normals, tan_pins;
  for (const auto& tan_pl : msphere.tan_planes) {
    tan_pl_normals.push_back(tan_pl.normal);
    tan_pins.push_back(tan_pl.tan_point);
  }
  auto mat_p_tan_elements =
      polyscope::registerPointCloud("Tangent Elements", tan_pins);
  mat_p_tan_elements->addVectorQuantity("normal", tan_pl_normals)
      ->setEnabled(true);
}

void MainGuiWindow::show_medial_mesh(const MedialMesh& mmesh,
                                     int given_mat_face_id,
                                     std::string mmname) {
  if (mmesh.vertices == nullptr) return;
  const auto& mspheres = *(mmesh.vertices);
  const auto& mfaces = mmesh.faces;

  std::vector<Vector3> mat_pos(mspheres.size());
  std::vector<double> mat_radius(mspheres.size());
  std::vector<aint3> mat_faces(mfaces.size(), {{0, 0, 0}});
  std::vector<double> mat_faces_importance(mfaces.size(), 2.0);
  std::vector<int> mat_faces_cnt(mfaces.size(), 0);

  // set true when given_mat_face_id != -1
  std::vector<bool> is_given_mat_face(mmesh.faces.size(), false);
  std::vector<int> mat_faces_mstruc_ids(mfaces.size(), -1);
  std::vector<Vector3> dual_segment_vs;
  std::vector<aint2> dual_segment_edge;
  // // for simple triangles
  // std::vector<Vector3> st_points;
  // std::vector<aint3> st_faces;
  // std::vector<Vector3> st_centroids;
  // std::vector<Vector3> st_normals;
  // std::vector<Vector3> st_nearest_p;

  // store mat vertices
  for (uint i = 0; i < mspheres.size(); i++) {
    mat_pos[i] = mspheres[i].center;
    mat_radius[i] = mspheres[i].radius;
  }

  // store mat faces
  // need to enbale culling for drawing ma faces twice
  for (uint f = 0; f < mfaces.size(); f++) {
    if (mfaces[f].is_deleted) continue;
    // if (mfaces[f].dup_cnt <= 1) continue;
    mat_faces_cnt[f] = mfaces[f].dup_cnt;
    mat_faces_importance[f] = mfaces[f].importance;
    mat_faces_mstruc_ids[f] = mfaces[f].mstruc_id;
    // draw facets
    uint lv = 0;
    for (const auto& v : mfaces[f].vertices_) {
      mat_faces[f][lv++] = v;
    }
  }

  // if this mat face is what we are looking for
  if (given_mat_face_id > 0 && given_mat_face_id < mmesh.faces.size()) {
    int f = given_mat_face_id;
    const auto& matf = mmesh.faces[f];
    is_given_mat_face[f] = true;
    for (int i = 0; i < matf.dual_edge_endpoints.size(); i++) {
      dual_segment_vs.push_back(matf.dual_edge_endpoints[i][0].first);
      dual_segment_vs.push_back(matf.dual_edge_endpoints[i][1].first);
      dual_segment_edge.push_back({{i * 2, i * 2 + 1}});
    }
    // // show simple triangles
    // for (int i = 0; i < 2; i++) {
    //   st_centroids.push_back(matf.st[i].centroid);
    //   st_normals.push_back(matf.st[i].normal);
    //   st_nearest_p.push_back(matf.st[i].nearest_point);
    //   for (int j = 0; j < 3; j++) st_points.push_back(matf.st[i].v[j]);
    //   st_faces.push_back({{i * 3, i * 3 + 1, i * 3 + 2}});
    // }

    matf.print_medial_face();
    for (const auto tid : matf.tets_) {
      mmesh.tets.at(tid).print_medial_tet();
    }
    for (int i = 0; i < 2; i++) {
      printf("mface %d has nearest surface fid %d, dist_to_sf %f\n", matf.fid,
             matf.st[i].nearest_sf_fid, matf.st[i].dist_to_sf);
    }
  }

  // Register medial mesh
  auto medial_mesh = polyscope::registerSurfaceMesh(mmname, mat_pos, mat_faces);
  medial_mesh->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
  medial_mesh->addVertexScalarQuantity("radius", mat_radius,
                                       polyscope::DataType::MAGNITUDE);
  medial_mesh->addFaceScalarQuantity("importance", mat_faces_importance)
      ->setEnabled(false);
  medial_mesh->addFaceScalarQuantity("dup_cnt", mat_faces_cnt)
      ->setEnabled(false);
  medial_mesh->addFaceScalarQuantity("mstruc_id", mat_faces_mstruc_ids)
      ->setEnabled(true);

  if (given_mat_face_id != -1) {
    medial_mesh->addFaceScalarQuantity("is_given_mat_face", is_given_mat_face)
        ->setEnabled(true);
    // polyscope::registerCurveNetwork("dual_segment", dual_segment_vs,
    //                                 dual_segment_edge);
    // polyscope::registerSurfaceMesh("mface_st", st_points, st_faces);
    // auto tmp = polyscope::registerPointCloud("st centroids", st_centroids);
    // tmp->addVectorQuantity("st normals", st_normals)->setEnabled(true);
    // polyscope::registerPointCloud("mface_st nearest_p", st_nearest_p);
  }
}

void MainGuiWindow::show_medial_edges(const MedialMesh& mmesh) {
  const auto& mspheres = *(mmesh.vertices);
  const auto& medges = mmesh.edges;
  std::vector<Vector3> mat_pos(mspheres.size());
  std::vector<aint2> mat_edges;
  std::vector<int> mat_edges_extf;
  // store mat vertices
  for (uint i = 0; i < mspheres.size(); i++) {
    mat_pos[i] = mspheres[i].center;
  }

  // store mat edges
  for (uint e = 0; e < medges.size(); e++) {
    if (medges[e].is_deleted) continue;
    mat_edges.push_back(medges[e].vertices_);
    mat_edges_extf.push_back(medges[e].is_extf);
  }

  auto me = polyscope::registerCurveNetwork("medial edges", mat_pos, mat_edges);
  me->setRadius(MainGuiWindow::edge_radius_rel, true /*isRelative*/);
  // me->addEdgeScalarQuantity("extf", mat_edges_extf)->setEnabled(true);
}