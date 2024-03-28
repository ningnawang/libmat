cmake_minimum_required(VERSION 3.17)

# geogram
if(NOT TARGET geogram)
    message(FATAL_ERROR "geogram not found in inputs" )
endif()
include(geogram) # check geogram.cmake

add_library(inputs STATIC 
    input_types.cxx
    input_types.h
    params.h
    sharp_feature_detection.h
    sharp_feature_detection.cxx
)

# add common headers
target_include_directories(inputs PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/inputs
)

target_link_libraries(inputs
    optimized geogram
)