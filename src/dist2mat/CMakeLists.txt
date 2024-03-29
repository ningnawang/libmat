cmake_minimum_required(VERSION 2.6)
add_library(dist2mat STATIC dist2mat.h dist2mat.cu)

# add common headers
target_include_directories(dist2mat PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
)

