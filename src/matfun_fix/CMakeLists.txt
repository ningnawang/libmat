cmake_minimum_required(VERSION 3.17)

add_library(matfun_fix STATIC 
    fix_common.cxx
    fix_common.h
    fix_topo.cxx
    fix_topo.h
    fix_geo.cxx
    fix_geo.h
    fix_geo_error.cxx
    fix_geo_error.h
    fix_extf.cxx
    fix_extf.h
    fix_intf.cxx
    fix_intf.h
)

# add lib
target_link_libraries(matfun_fix
    optimized dist2mat
    optimized inputs
    optimized matbase
    optimized matfun
)

# add common headers
target_include_directories(matfun_fix PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/matbase
    ${PROJECT_SOURCE_DIR}/src/matfun
    ${PROJECT_SOURCE_DIR}/src/dist2mat
)




