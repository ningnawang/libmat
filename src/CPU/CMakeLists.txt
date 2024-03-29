cmake_minimum_required(VERSION 3.17)

add_library(cpu_fallback STATIC 
    convex_cell.h
    convex_cell.cpp
    Delaunay_psm.h
    Delaunay_psm.cpp
    voronoi_fallback.h
    voronoi_fallback.cpp
)

# add lib
target_link_libraries(cpu_fallback
    optimized rpd3d
)

# add common headers
target_include_directories(cpu_fallback PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/rpd3d
)

