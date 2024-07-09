cmake_minimum_required(VERSION 3.16)
add_library(rpd3d STATIC   
    convex_cell.h
    convex_cell.cu
    voronoi.h
    voronoi.cu
    stopwatch.h
)
target_compile_options(rpd3d PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
    --generate-line-info
    --use_fast_math
    >)

# add lib
add_subdirectory(kNN-CUDA)
# add common headers
target_include_directories(rpd3d PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/rpd3d_base
)

target_link_libraries(rpd3d
    optimized rpd3d_base    
    debug kNN-CUDAd
    optimized kNN-CUDA
)