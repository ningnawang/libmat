cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 17)

add_library(matfun STATIC 
    shrinking.h
    shrinking.cxx
    updating.h
    updating.cxx
    thinning.h
    thinning.cxx
)

# add lib
target_link_libraries(matfun
    optimized inputs
    optimized matbase
)

# add common headers
target_include_directories(matfun PUBLIC 
    ${LIBMAT_MODULE_INCLUDE}
    ${LIBMAT_MODULE_SRC}/matbase
)

