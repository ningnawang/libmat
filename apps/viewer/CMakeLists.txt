# ###############################################################################
# Polyscope viewer: load a tet mesh + a MATFP .ma file and rescale it to the
# MATTopo convention (was tests_cxx/). Enable with -DLIBMAT_BUILD_VIEWER=ON.
#
#   ./build/apps/viewer/libmat_viewer <tet_mesh.msh> <medial_mesh.ma>
# ###############################################################################
include(FetchContent)
if(NOT TARGET polyscope)
    FetchContent_Declare(polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG eb07f8acaf5c8dba30de9587f84b14ba0c411ca1
    )
    FetchContent_MakeAvailable(polyscope)
endif()

add_executable(libmat_viewer main_cxx.cxx main_gui_cxx.cxx)
target_link_libraries(libmat_viewer PRIVATE libmat::libmat polyscope)
if(NOT MSVC)
    # same flag as the former tests_cxx build
    target_compile_options(libmat_viewer PRIVATE -Wno-narrowing)
endif()
