# ###############################################################################
project(MAT_MODULES_CXX LANGUAGES CXX)
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# ###############################################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wno-narrowing")

# ###############################################################################
if(IS_USE_GPU)
    find_package(CUDA REQUIRED)
    link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)
    include_directories("${CUDA_INCLUDE_DIRS}")
endif(IS_USE_GPU)

# find_package(OpenMP)
# if(OPENMP_FOUND)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

find_package(ZLIB REQUIRED)

find_package(CGAL)

# Get Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    # message(STATUS "${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# add dependencies
list(APPEND CMAKE_MODULE_PATH ${MAT_MODULE_CMAKE})
option(AUTO_DOWNLOAD "Auto download dependencies" ON)
include(rpdDependencies)
include(geogram) # check geogram.cmake
add_subdirectory(${MAT_MODULE_EXTERNAL}/polyscope ${CMAKE_CURRENT_BINARY_DIR}/polyscope)
add_subdirectory(${MAT_MODULE_EXTERNAL}/mshloader ${CMAKE_CURRENT_BINARY_DIR}/mshloader)
if(IS_USE_GPU)
    add_subdirectory(${MAT_MODULE_SRC}/dist2mat ${CMAKE_CURRENT_BINARY_DIR}/dist2mat)
    add_subdirectory(${MAT_MODULE_SRC}/rpd3d ${CMAKE_CURRENT_BINARY_DIR}/rpd3d)
    add_subdirectory(${MAT_MODULE_SRC}/rpd3d_api ${CMAKE_CURRENT_BINARY_DIR}/rpd3d_api)
    add_subdirectory(${MAT_MODULE_SRC}/matfun_fix ${CMAKE_CURRENT_BINARY_DIR}/matfun_fix)
    add_subdirectory(${MAT_MODULE_SRC}/IO/IO_CUDA ${CMAKE_CURRENT_BINARY_DIR}/IO/IO_CUDA)
endif(IS_USE_GPU)
add_subdirectory(${MAT_MODULE_SRC}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs)
add_subdirectory(${MAT_MODULE_SRC}/matbase ${CMAKE_CURRENT_BINARY_DIR}/matbase)
add_subdirectory(${MAT_MODULE_SRC}/matfun ${CMAKE_CURRENT_BINARY_DIR}/matfun)
add_subdirectory(${MAT_MODULE_SRC}/IO/IO_CXX ${CMAKE_CURRENT_BINARY_DIR}/IO/IO_CXX)
add_subdirectory(${MAT_MODULE_SRC}/rpd3d_base ${CMAKE_CURRENT_BINARY_DIR}/rpd3d_base)


# ###############################################################################
set(SOURCE_LIST
	main_cxx.cxx
	main_gui_cxx.cxx
)
add_executable(${PROJECT_NAME} ${SOURCE_LIST})
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${MAT_MODULE_INCLUDE}>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(${PROJECT_NAME} PUBLIC
    ### mat_modules lib
    optimized mat_modules_cxx
    ### other libs
    optimized mshloader
    optimized polyscope
    optimized geogram
    optimized nlohmann_json::nlohmann_json
    optimized ZLIB::ZLIB
    optimized CGAL::CGAL
    optimized Eigen3::Eigen
    optimized igl::core
)

if(IS_USE_GPU)
    target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        --generate-line-info
        --use_fast_math
        >)
endif(IS_USE_GPU)