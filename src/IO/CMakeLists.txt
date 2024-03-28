cmake_minimum_required(VERSION 3.17)

add_library(IO STATIC 
    io_utils.hpp
    io_utils.cpp
    io.h
    io.cxx
)

# add lib
target_link_libraries(IO
    optimized inputs
    optimized rpd3d
    optimized matbase
    optimized nlohmann_json::nlohmann_json
    optimized mshloader
)

# add common headers
target_include_directories(IO PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/extern/mshloader
    ${PROJECT_SOURCE_DIR}/src/external
    ${PROJECT_SOURCE_DIR}/src/matbase
)

