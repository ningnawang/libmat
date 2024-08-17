cmake_minimum_required(VERSION 3.17)

# geogram
if(NOT TARGET geogram)
    message(FATAL_ERROR "geogram not found in inputs" )
endif()

add_library(inputs STATIC 
    input_types.cxx
    input_types.h
    params.h
    sharp_feature_detection.h
    sharp_feature_detection.cxx
)


# add common headers
target_include_directories(inputs PUBLIC 
    ${MAT_MODULE_INCLUDE}
    ${MAT_MODULE_SRC}/inputs
)

target_link_libraries(inputs
    optimized geogram
)