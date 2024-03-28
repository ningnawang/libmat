cmake_minimum_required(VERSION 3.16)
add_library(rpd3d STATIC  
    voronoi_common.h    
    convex_cell.h
    convex_cell.cu
    voronoi.h
    voronoi.cu
    stopwatch.h
    voronoi_defs.h
    voronoi_defs.cxx
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
)

target_link_libraries(rpd3d
    debug kNN-CUDAd
    optimized kNN-CUDA
)