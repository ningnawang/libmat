# ###############################################################################
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(libmat)

# ###############################################################################
option(LIBMAT_IS_USE_GPU "Enable GPU for RPD" OFF)
option(LIBMAT_IS_BUILD_TEST "Build tests, either CXX or CUDA" ON)
message("LIBMAT_IS_USE_GPU: ${LIBMAT_IS_USE_GPU}")
message("LIBMAT_IS_BUILD_TEST: ${LIBMAT_IS_BUILD_TEST}")

# ###############################################################################
set(LIBMAT_MODULE_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(LIBMAT_MODULE_CMAKE "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(LIBMAT_MODULE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(LIBMAT_MODULE_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(LIBMAT_PROJECT_NAME "libmat")

# ###############################################################################
# add dependencies if NOT exists
# 1. geogram
# 2. libigl
# 3. polyscope
# 4. nlohmann_json
list(APPEND CMAKE_MODULE_PATH ${LIBMAT_MODULE_CMAKE})
option(AUTO_DOWNLOAD "Auto download dependencies" ON)
include(rpdDependencies)

# ###############################################################################
# add subdirectories
add_subdirectory(${LIBMAT_MODULE_EXTERNAL}/mshloader)
add_subdirectory(src/inputs)
add_subdirectory(src/IO/IO_CXX)
add_subdirectory(src/rpd3d_base)
add_subdirectory(src/matbase)
add_subdirectory(src/matfun)
if(LIBMAT_IS_USE_GPU)
    add_subdirectory(src/rpd3d)
    add_subdirectory(src/rpd3d_api)
    add_subdirectory(src/dist2mat)
    add_subdirectory(src/matfun_fix)
    add_subdirectory(src/IO/IO_CUDA)
endif(LIBMAT_IS_USE_GPU)

# ###############################################################################
# create the interface library
add_library(${LIBMAT_PROJECT_NAME} INTERFACE)
target_compile_features(${LIBMAT_PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(${LIBMAT_PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${LIBMAT_MODULE_INCLUDE}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_directories(${LIBMAT_PROJECT_NAME} INTERFACE "${CMAKE_BINARY_DIR}/lib")
target_include_directories(${LIBMAT_PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${LIBMAT_PROJECT_NAME} INTERFACE
    optimized inputs
    optimized IO_CXX
    optimized rpd3d_base
    optimized matbase
    optimized matfun
)
if(LIBMAT_IS_USE_GPU)
    target_link_libraries(${LIBMAT_PROJECT_NAME} INTERFACE
        general cuda
        general cublas
        ###
        optimized rpd3d
        optimized rpd3d_api
        optimized dist2mat
        optimized matfun_fix
        optimized IO_CUDA
    )
    target_link_libraries(${LIBMAT_PROJECT_NAME} INTERFACE ${CUDA_LIBRARIES})
endif(LIBMAT_IS_USE_GPU)


# ###############################################################################
# build test
if(LIBMAT_IS_BUILD_TEST)
    message("Building CXX tests ...")
    add_subdirectory(tests_cxx)
endif(LIBMAT_IS_BUILD_TEST)