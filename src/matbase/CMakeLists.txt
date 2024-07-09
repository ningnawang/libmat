cmake_minimum_required(VERSION 3.17)

add_library(matbase STATIC 
    medial_mesh.h
    medial_mesh.cxx
    medial_primitives.h
    medial_primitives.cxx
    medial_sphere.h
    medial_sphere.cxx
)

# add lib
target_link_libraries(matbase
    optimized inputs
)

# add common headers
target_include_directories(matbase PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
)

