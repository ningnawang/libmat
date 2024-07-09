cmake_minimum_required(VERSION 3.16)
add_library(rpd3d_api STATIC  
    rpd_api.h
    rpd_api.cxx
    rpd_update.h
    rpd_update.cxx
)

# add lib
target_link_libraries(rpd3d_api
    optimized rpd3d_base
    optimized inputs
    optimized matbase
)

# add common headers
target_include_directories(rpd3d_api PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/inputs
    ${PROJECT_SOURCE_DIR}/src/matbase
    ${PROJECT_SOURCE_DIR}/src/rpd3d_base
)
