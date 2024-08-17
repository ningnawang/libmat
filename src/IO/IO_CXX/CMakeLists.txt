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
    optimized mshloader
)

# add common headers
target_include_directories(IO_CXX PUBLIC 
    ${MAT_MODULE_INCLUDE}
    ${MAT_MODULE_EXTERNAL}/mshloader
    ${MAT_MODULE_SRC}/external
    ${MAT_MODULE_SRC}/matbase
)
