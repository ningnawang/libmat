# ###############################################################################
# Third-party dependencies of libmat
#   - geogram  (v1.7.9 downloaded, built statically; set LIBMAT_GEOGRAM_DIR
#               to use an existing checkout)
#   - CGAL     (system, header-only; needs Boost + GMP/MPFR)
#   - Eigen3   (system, downloaded if missing)
#   - nlohmann_json (downloaded, CUDA build only)
# ###############################################################################
include(FetchContent)

# ###############################################################################
# geogram
if(NOT TARGET geogram)
    # 1.7.9 is API compatible with 1.7.5 (used by MATFP/MATTopo) and also
    # builds on Apple silicon
    set(LIBMAT_GEOGRAM_DIR "" CACHE PATH
        "Path to the geogram sources (downloaded if empty)")
    if(NOT EXISTS "${LIBMAT_GEOGRAM_DIR}/src/lib/geogram/basic/common.h")
        if(NOT LIBMAT_AUTO_DOWNLOAD)
            message(FATAL_ERROR "geogram not found in ${LIBMAT_GEOGRAM_DIR}")
        endif()
        FetchContent_Declare(geogram
            GIT_REPOSITORY https://github.com/alicevision/geogram.git
            GIT_TAG v1.7.9
            GIT_SHALLOW TRUE
            SOURCE_SUBDIR do-not-add-subdirectory
        )
        FetchContent_MakeAvailable(geogram)
        set(LIBMAT_GEOGRAM_DIR "${geogram_SOURCE_DIR}" CACHE PATH "" FORCE)
    endif()
    include(geogram)
endif()

# ###############################################################################
# CGAL (header-only since 5.0)
if(NOT TARGET CGAL::CGAL)
    find_package(CGAL REQUIRED)
endif()

# ###############################################################################
# Eigen
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.3 CONFIG QUIET)
    if(NOT TARGET Eigen3::Eigen)
        if(NOT LIBMAT_AUTO_DOWNLOAD)
            message(FATAL_ERROR "Eigen3 not found")
        endif()
        FetchContent_Declare(eigen
            URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
            SOURCE_SUBDIR do-not-add-subdirectory
        )
        FetchContent_MakeAvailable(eigen)
        add_library(Eigen3_Eigen INTERFACE)
        add_library(Eigen3::Eigen ALIAS Eigen3_Eigen)
        target_include_directories(Eigen3_Eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
    endif()
endif()

# ###############################################################################
# nlohmann_json (only used by the CUDA IO module)
if(LIBMAT_WITH_CUDA AND NOT TARGET nlohmann_json::nlohmann_json)
    FetchContent_Declare(json
        URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
    )
    FetchContent_MakeAvailable(json)
endif()
