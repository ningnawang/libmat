# ###############################################################################
# Third-party dependencies of libmat
#   - geogram  (v1.7.9 downloaded, built statically; set LIBMAT_GEOGRAM_DIR
#               to use an existing checkout)
#   - CGAL     (system, header-only, downloaded if missing; needs Boost
#               headers, and uses GMP/MPFR unless CGAL_DISABLE_GMP=ON)
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
    find_package(CGAL QUIET)
    if(NOT CGAL_FOUND)
        if(NOT LIBMAT_AUTO_DOWNLOAD)
            message(FATAL_ERROR "CGAL not found")
        endif()
        # 5.6 is the last release that accepts Boost 1.66 (manylinux_2_28)
        FetchContent_Declare(cgal
            URL https://github.com/CGAL/cgal/releases/download/v5.6.2/CGAL-5.6.2-library.tar.xz
            URL_HASH SHA256=2d93d96c8c73c73c9999556d512775648d26e1e04035ceb958208d25e6e2d10f
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            SOURCE_SUBDIR do-not-add-subdirectory
        )
        FetchContent_MakeAvailable(cgal)
        # without GMP/MPFR CGAL falls back to Boost.Multiprecision
        if(NOT DEFINED CGAL_DISABLE_GMP)
            find_path(LIBMAT_GMP_INCLUDE_DIR gmp.h)
            find_library(LIBMAT_GMP_LIBRARY gmp)
            find_library(LIBMAT_MPFR_LIBRARY mpfr)
            if(NOT LIBMAT_GMP_INCLUDE_DIR OR NOT LIBMAT_GMP_LIBRARY OR NOT LIBMAT_MPFR_LIBRARY)
                message(STATUS "GMP/MPFR not found, building CGAL without them")
                set(CGAL_DISABLE_GMP ON)
            endif()
        endif()
        find_package(CGAL REQUIRED PATHS ${cgal_SOURCE_DIR} NO_DEFAULT_PATH)
    endif()
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
        # 3.12 is the first release whose CMakeLists accepts CMake >= 4
        URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
        URL_HASH SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(json)
endif()
