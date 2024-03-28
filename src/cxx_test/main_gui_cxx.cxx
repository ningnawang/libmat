
#include "main_gui_cxx.h"

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>

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
  delete all_medial_spheres;
  // delete tet_vertices;
  // delete tet_indices;
  // delete sf_vertices;
  // delete sf_faces;
  this->instance_ = nullptr;
}
void MainGuiWindow::set_sf_mesh(GEO::Mesh& _sf_mesh) { sf_mesh = &_sf_mesh; }
void MainGuiWindow::set_all_medial_spheres(
    std::vector<MedialSphere>& _all_medial_spheres) {
  all_medial_spheres = &_all_medial_spheres;
}
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
  instance_->show_surface_mesh(*(instance_->sf_mesh), -1 /*given_sf_face_id*/);

  // Show the GUI
  polyscope::show();
}

// Your callback functions
void MainGuiWindow::callbacks() {
  ImGui::PushItemWidth(100);

  // // Surface
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
    tan_pins.push_back(tan_pl.points[0]);
  }
  auto mat_p_tan_elements =
      polyscope::registerPointCloud("Tangent Elements", tan_pins);
  mat_p_tan_elements->addVectorQuantity("normal", tan_pl_normals)
      ->setEnabled(true);
}