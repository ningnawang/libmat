# ###############################################################################
project(MAT_MODULES_CXX LANGUAGES CXX)
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# ###############################################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wno-narrowing")

# ###############################################################################
find_package(ZLIB REQUIRED)

find_package(CGAL REQUIRED)

# Get Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    # message(STATUS "${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# add dependencies
list(APPEND CMAKE_MODULE_PATH ${MAT_MODULE_CMAKE})
option(AUTO_DOWNLOAD "Auto download dependencies" ON)
include(rpdDependencies)
add_subdirectory(${MAT_MODULE_EXTERNAL}/polyscope ${CMAKE_CURRENT_BINARY_DIR}/polyscope)
add_subdirectory(${MAT_MODULE_EXTERNAL}/mshloader ${CMAKE_CURRENT_BINARY_DIR}/mshloader)
add_subdirectory(${MAT_MODULE_SRC}/inputs ${CMAKE_CURRENT_BINARY_DIR}/inputs)
add_subdirectory(${MAT_MODULE_SRC}/matbase ${CMAKE_CURRENT_BINARY_DIR}/matbase)
add_subdirectory(${MAT_MODULE_SRC}/matfun ${CMAKE_CURRENT_BINARY_DIR}/matfun)
add_subdirectory(${MAT_MODULE_SRC}/IO/IO_CXX ${CMAKE_CURRENT_BINARY_DIR}/IO/IO_CXX)
add_subdirectory(${MAT_MODULE_SRC}/rpd3d_base ${CMAKE_CURRENT_BINARY_DIR}/rpd3d_base)

# ###############################################################################
set(SOURCE_LIST
	main_cxx.cxx
	main_gui_cxx.cxx
)
add_executable(${PROJECT_NAME} ${SOURCE_LIST})
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${MAT_MODULE_INCLUDE}>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(${PROJECT_NAME} PUBLIC
    ### mat_modules lib
    optimized mat_modules
    ### other libs
    optimized polyscope
    optimized geogram
    optimized ZLIB::ZLIB
    optimized CGAL::CGAL
    optimized Eigen3::Eigen
    optimized igl::core
)