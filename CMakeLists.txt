# ###############################################################################
# LibMAT: a C++ (+ CUDA) library for the 3D medial axis transform,
# with Python bindings.
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release                # CPU only
#   cmake -B build -DLIBMAT_WITH_CUDA=ON                     # + CUDA modules
#   cmake -B build -DLIBMAT_BUILD_PYTHON=ON                  # + python module
#
# Targets:
#   libmat::core    CPU library (always built)
#   libmat::cuda    CUDA library (LIBMAT_WITH_CUDA=ON)
#   libmat::libmat  interface target linking everything that was built
#   _libmat         python extension module (LIBMAT_BUILD_PYTHON=ON)
# ###############################################################################
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(libmat VERSION 0.1.0 LANGUAGES C CXX)

# ###############################################################################
# options
option(LIBMAT_WITH_CUDA "Build the CUDA modules (RPD, dist2mat, matfun_fix)" OFF)
option(LIBMAT_BUILD_PYTHON "Build the python bindings" OFF)
option(LIBMAT_BUILD_VIEWER "Build the polyscope viewer app in apps/viewer" OFF)
option(LIBMAT_AUTO_DOWNLOAD "Download missing dependencies" ON)

# backward compatibility with the old option names
if(DEFINED LIBMAT_IS_USE_GPU)
    set(LIBMAT_WITH_CUDA ${LIBMAT_IS_USE_GPU} CACHE BOOL "" FORCE)
endif()
if(DEFINED LIBMAT_IS_BUILD_TEST)
    set(LIBMAT_BUILD_VIEWER ${LIBMAT_IS_BUILD_TEST} CACHE BOOL "" FORCE)
endif()

message(STATUS "LIBMAT_WITH_CUDA:     ${LIBMAT_WITH_CUDA}")
message(STATUS "LIBMAT_BUILD_PYTHON:  ${LIBMAT_BUILD_PYTHON}")
message(STATUS "LIBMAT_BUILD_VIEWER:  ${LIBMAT_BUILD_VIEWER}")

# ###############################################################################
# global settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(LIBMAT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(LIBMAT_EXTERNAL_DIR "${LIBMAT_ROOT}/extern")
list(APPEND CMAKE_MODULE_PATH "${LIBMAT_ROOT}/cmake")

if(LIBMAT_WITH_CUDA)
    # must be set before enable_language(CUDA), which picks a default otherwise
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
            set(CMAKE_CUDA_ARCHITECTURES native)
        else()
            set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86)
        endif()
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# geogram, CGAL, Eigen (+ nlohmann_json for CUDA)
include(libmat_dependencies)

# generated config header
set(LIBMAT_WITH_CUDA_VALUE 0)
if(LIBMAT_WITH_CUDA)
    set(LIBMAT_WITH_CUDA_VALUE 1)
endif()
configure_file(cmake/config.h.in "${CMAKE_CURRENT_BINARY_DIR}/include/libmat/config.h" @ONLY)

# ###############################################################################
# libmat::core (CPU)
add_subdirectory(${LIBMAT_EXTERNAL_DIR}/mshloader ${CMAKE_CURRENT_BINARY_DIR}/mshloader)

set(LIBMAT_CORE_SOURCES
    # legacy modules
    src/common/predicates.cxx
    src/inputs/input_types.cxx
    src/inputs/mesh_AABB.cxx
    src/inputs/sharp_feature_detection.cxx
    src/matbase/medial_mesh.cxx
    src/matbase/medial_primitives.cxx
    src/matbase/medial_sphere.cxx
    src/matfun/shrinking.cxx
    src/matfun/thinning.cxx
    src/matfun/updating.cxx
    src/io/io.cxx
    src/rpd3d/rpd_update.cxx
    src/rpd3d/triangulation.cxx
    src/rpd3d/voronoi_defs.cxx
    # functional API (one function per file)
    src/api/internal/surface_context.cpp
    src/api/detect_sharp_features.cpp
    src/api/distance_to_medial_mesh.cpp
    src/api/extract_surface.cpp
    src/api/normalize_vertices.cpp
    src/api/read_ma.cpp
    src/api/read_surface_mesh.cpp
    src/api/read_tet_mesh.cpp
    src/api/regular_triangulation.cpp
    src/api/shrink_spheres.cpp
    src/api/slab_tangent_triangles.cpp
    src/api/thin_medial_mesh.cpp
    src/api/write_ma.cpp
)
if(NOT LIBMAT_WITH_CUDA)
    # CPU stubs that raise a clear error for CUDA-only functions
    list(APPEND LIBMAT_CORE_SOURCES src/api/cuda_stubs.cpp)
endif()

add_library(libmat_core STATIC ${LIBMAT_CORE_SOURCES})
add_library(libmat::core ALIAS libmat_core)
target_include_directories(libmat_core
    PUBLIC
        $<BUILD_INTERFACE:${LIBMAT_ROOT}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    PRIVATE
        ${LIBMAT_ROOT}/src
        ${LIBMAT_EXTERNAL_DIR}/mshloader
)
target_link_libraries(libmat_core
    PUBLIC geogram CGAL::CGAL Eigen3::Eigen
    PRIVATE mshloader
)
if(MSVC)
    target_compile_options(libmat_core PRIVATE /bigobj)
else()
    # the legacy code base is not warning-clean, keep the build log readable,
    # but never let a missing return through: gcc turns it into a crash
    target_compile_options(libmat_core PRIVATE
        -Wno-narrowing -Wno-format -Wno-format-extra-args -Werror=return-type)
endif()

# ###############################################################################
# libmat::cuda
if(LIBMAT_WITH_CUDA)
    set(LIBMAT_CUDA_SOURCES
        src/cuda/rpd3d/convex_cell.cu
        src/cuda/rpd3d/voronoi.cu
        src/cuda/rpd3d/knncuda.cu
        src/cuda/rpd3d/rpd_api.cxx
        src/cuda/dist2mat/dist2mat.cu
        src/cuda/matfun_fix/fix_common.cxx
        src/cuda/matfun_fix/fix_extf.cxx
        src/cuda/matfun_fix/fix_geo.cxx
        src/cuda/matfun_fix/fix_geo_error.cxx
        src/cuda/matfun_fix/fix_intf.cxx
        src/cuda/matfun_fix/fix_topo.cxx
        src/cuda/io/io_cuda.cxx
        src/cuda/io/io_utils.cpp
        # functional API
        src/cuda/api/distance_to_medial_mesh_cuda.cpp
        src/cuda/api/medial_mesh_from_spheres.cpp
    )
    add_library(libmat_cuda STATIC ${LIBMAT_CUDA_SOURCES})
    add_library(libmat::cuda ALIAS libmat_cuda)
    # same flags as the original rpd3d target
    set_source_files_properties(
        src/cuda/rpd3d/convex_cell.cu src/cuda/rpd3d/voronoi.cu
        PROPERTIES COMPILE_OPTIONS "--generate-line-info;--use_fast_math")
    target_include_directories(libmat_cuda PRIVATE ${LIBMAT_ROOT}/src)
    target_link_libraries(libmat_cuda
        PUBLIC libmat_core CUDA::cudart CUDA::cublas CUDA::cuda_driver
        PRIVATE nlohmann_json::nlohmann_json
    )
endif()

# ###############################################################################
# libmat::libmat (everything that was built)
add_library(libmat INTERFACE)
add_library(libmat::libmat ALIAS libmat)
target_link_libraries(libmat INTERFACE libmat_core)
if(LIBMAT_WITH_CUDA)
    target_link_libraries(libmat INTERFACE libmat_cuda)
endif()

# ###############################################################################
# python bindings
if(LIBMAT_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# ###############################################################################
# polyscope viewer (C++)
if(LIBMAT_BUILD_VIEWER)
    add_subdirectory(apps/viewer)
endif()
