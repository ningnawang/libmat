cmake_minimum_required(VERSION 3.17)

add_library(IO_CXX STATIC 
    io.h
    io.cxx
    Predicates.cxx
    Predicates.hpp
)

# add lib
target_link_libraries(IO_CXX
    optimized inputs
    optimized matbase
    optimized nlohmann_json::nlohmann_json
    optimized mshloader
)

# add common headers
target_include_directories(IO_CXX PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/extern/mshloader
    ${PROJECT_SOURCE_DIR}/src/external
    ${PROJECT_SOURCE_DIR}/src/matbase
    ${PROJECT_SOURCE_DIR}/src/rpd3d
)
