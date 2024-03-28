#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "../io.h"
#include "../params.h"
#include "../shrinking.h"
#include "../updating.h"
#include "main_gui_cxx.h"

int main(int argc, char** argv) {
  if (argc < 1) {
    std::cerr << "Usage: " << argv[0]
              << " <tet_mesh.tet/vtk> (e.g.: " << argv[0]
              << " ../data/joint.tet)" << std::endl;
    return 1;
  }

  std::srand(RAN_SEED);  // set random seed
  std::string tet_mesh_path = argv[1];
  // read surface file from tetwild/ftetwild
  // this will make sure all geogram algoriths can be used
  GEO::initialize();
  GEO::CmdLine::import_arg_group("algo");

  // ftetwild/tetwild will return "_surf.obj" for surface triangle mesh
  std::string surface_path = get_other_file_path(tet_mesh_path);
  std::cout << "surface path: " << surface_path << std::endl;
  GEO::Mesh sf_mesh;
  if (!is_file_exist(surface_path)) {
    std::cout << "surface path: " << surface_path
              << " not exist, load from tet_mesh: " << tet_mesh_path
              << std::endl;
    std::vector<float> tet_vertices;  // tet vertices
    std::vector<int> tet_indices;     // tet 4 indices of vertices
    // please set normalize to true for model joint.tet/.xyz
    if (!load_tet(tet_mesh_path, tet_vertices, tet_indices,
                  false /*normalize*/)) {
      std::cerr << tet_mesh_path << ": could not load file" << std::endl;
      return 1;
    }
    std::vector<std::array<float, 3>> sf_vertices;
    std::vector<std::array<int, 3>> sf_faces;
    get_surface_from_tet(tet_vertices, tet_indices, sf_vertices, sf_faces);
    load_sf_mesh_from_internal(sf_vertices, sf_faces, sf_mesh);
    if (!save_sf_mesh(surface_path, sf_mesh)) {
      std::cerr << surface_path << ": could not save surface file" << std::endl;
      return 1;
    }
  } else {
    if (!load_surface_mesh(surface_path, sf_mesh)) {
      std::cerr << surface_path << ": could not load surface file" << std::endl;
      return 1;
    }
  }
  std::cout << "Done loading surface mesh" << std::endl;

  // std::vector<MedialSphere> all_medial_spheres;
  // AABBWrapper aabb_wrapper;
  // pre_and_init_aabb(sf_mesh, aabb_wrapper);
  // init_and_shrink(sf_mesh, aabb_wrapper, all_medial_spheres, 10,
  //                 -1 /*itr_limit*/, true);
  // init_shrink_and_update(sf_mesh, aabb_wrapper, all_medial_spheres, 10,
  //                        true /*is_debug*/);

  // show Gui
  MainGuiWindow main_gui;
  main_gui.set_sf_mesh(sf_mesh);
  main_gui.set_all_medial_spheres(all_medial_spheres);
  // main_gui.set_tet_mesh(tet_vertices, tet_indices);
  // main_gui.set_sf_mesh_internal(sf_vertices, sf_faces);
  main_gui.show();

  return 0;
}
