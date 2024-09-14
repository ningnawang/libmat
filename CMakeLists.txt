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
# build test, before add subdirectories
if(LIBMAT_IS_BUILD_TEST)
    message("Building CXX tests ...")
    add_subdirectory(tests_cxx)
endif(LIBMAT_IS_BUILD_TEST)

# ###############################################################################
# add subdirectories
#
# we need if here because add_subdirectory(tests_cxx) will add the targets
if (NOT TARGET mshloader)
    add_subdirectory(${LIBMAT_MODULE_EXTERNAL}/mshloader)
endif()
if (NOT TARGET inputs)
    add_subdirectory(src/inputs)
endif()
if (NOT TARGET IO_CXX)
    add_subdirectory(src/IO/IO_CXX)
endif()
if (NOT TARGET rpd3d_base)
    add_subdirectory(src/rpd3d_base)
endif()
if (NOT TARGET matbase)
    add_subdirectory(src/matbase)
endif()
if (NOT TARGET matfun)
    add_subdirectory(src/matfun)
endif()
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
