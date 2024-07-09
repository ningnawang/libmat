cmake_minimum_required(VERSION 3.17)

add_library(IO_CUDA STATIC 
    io_utils.hpp
    io_utils.cpp
    io_cuda.h
    io_cuda.cxx
)

# add lib
target_link_libraries(IO_CUDA  
    optimized inputs
    optimized matbase
    optimized rpd3d
    optimized IO_CXX
    optimized nlohmann_json::nlohmann_json
)

# add common headers
target_include_directories(IO_CUDA PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/external
    ${PROJECT_SOURCE_DIR}/src/matbase
    ${PROJECT_SOURCE_DIR}/src/rpd3d
    ${PROJECT_SOURCE_DIR}/src/IO/IO_CXX
)

