cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)

# geogram
if(NOT TARGET geogram)
    message(FATAL_ERROR "geogram not found in inputs")
endif()
#CGAL
if(NOT TARGET CGAL)
    message(FATAL_ERROR "CGAL not found in inputs")
endif()
find_package(CGAL REQUIRED)
message("CGAL Version: ${CGAL_VERSION}")


add_library(inputs STATIC 
    input_types.cxx
    input_types.h
    params.h
    sharp_feature_detection.h
    sharp_feature_detection.cxx
    mesh_AABB.h
    mesh_AABB.cxx
    CGAL_poly.h
)


# add common headers
target_include_directories(inputs PUBLIC 
    ${LIBMAT_MODULE_INCLUDE}
    ${LIBMAT_MODULE_SRC}/inputs
)

target_link_libraries(inputs
    optimized geogram
    optimized CGAL::CGAL
)