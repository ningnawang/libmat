# ###############################################################################
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(mat_modules)

# ###############################################################################
option(IS_USE_GPU "Enable GPU for RPD" OFF)
option(IS_BUILD_TEST "Build tests, either CXX or CUDA" ON)

# ###############################################################################
set(MAT_MODULE_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(MAT_MODULE_CMAKE "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(MAT_MODULE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(MAT_MODULE_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(MAT_PROJECT_NAME "mat_modules_cxx")
if(IS_USE_GPU)
    set(MAT_PROJECT_NAME "mat_modules_cuda")
endif()

# ###############################################################################
if(IS_BUILD_TEST)
    if(IS_USE_GPU)
        message("Building CUDA tests ...")
        # add_subdirectory(tests_cuda)
    else()
        message("Building CXX tests ...")
        add_subdirectory(tests_cxx)
    endif(IS_USE_GPU)
endif(IS_BUILD_TEST)

# ###############################################################################
add_library(${MAT_PROJECT_NAME} INTERFACE)
target_include_directories(${MAT_PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${MAT_MODULE_INCLUDE}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_directories(${MAT_PROJECT_NAME} INTERFACE "${CMAKE_BINARY_DIR}/lib")
target_link_libraries(${MAT_PROJECT_NAME} INTERFACE
    optimized IO_CXX
    optimized rpd3d_base
    optimized inputs
    optimized matbase
    optimized matfun
)
if(IS_USE_GPU)
    target_link_libraries(${MAT_PROJECT_NAME} INTERFACE
        general cuda
        general cublas
        optimized rpd3d
        optimized rpd3d_api
        optimized dist2mat
        optimized matfun_fix
        optimized IO_CUDA
    )
    target_link_libraries(${MAT_PROJECT_NAME} general ${CUDA_LIBRARIES})
endif(IS_USE_GPU)
