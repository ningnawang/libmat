# ninwang: we do not use Eigen
cmake_minimum_required(VERSION 2.6)

add_library(mshloader STATIC MshLoader.cxx MshLoader.h) 

set_target_properties(mshloader PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)