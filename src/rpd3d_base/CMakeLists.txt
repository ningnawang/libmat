cmake_minimum_required(VERSION 3.16)

find_package(CGAL)

add_library(rpd3d_base STATIC  
    voronoi_common.h   
    voronoi_defs.h
    voronoi_defs.cxx
    triangulation.h
    triangulation.cxx
    rpd_update.h
    rpd_update.cxx
)

# add lib
target_link_libraries(rpd3d_base PUBLIC
    optimized inputs
    optimized CGAL::CGAL
    optimized matbase
)

# add common headers
target_include_directories(rpd3d_base PUBLIC
    ${LIBMAT_MODULE_INCLUDE}
    ${LIBMAT_MODULE_SRC}/inputs
    ${LIBMAT_MODULE_SRC}/matbase
)
