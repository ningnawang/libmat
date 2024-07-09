cmake_minimum_required(VERSION 3.16)

find_package(CGAL)

add_library(rpd3d_base STATIC  
    voronoi_common.h   
    voronoi_defs.h
    voronoi_defs.cxx
    triangulation.h
    triangulation.cxx
)

# add lib
target_link_libraries(rpd3d_base
    optimized inputs
    optimized CGAL::CGAL
    optimized matbase
)

# add common headers
target_include_directories(rpd3d_base PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/inputs
    ${PROJECT_SOURCE_DIR}/src/matbase
)
