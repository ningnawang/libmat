# ###############################################################################
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# ###############################################################################
option(USE_GPU "Enable GPU for RPD" ON) #ON by default

# ###############################################################################
if(USE_GPU)
    set(CMAKE_CUDA_ARCHITECTURES 61 75)
    project(MAT_MODULES LANGUAGES CXX CUDA)
else()
    project(MAT_MODULES LANGUAGES CXX)
endif(USE_GPU)

# ###############################################################################
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# ###############################################################################
if(USE_GPU)
    find_package(CUDA REQUIRED)
    link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)
    include_directories("${CUDA_INCLUDE_DIRS}")
endif(USE_GPU)

# find_package(OpenMP)
# if(OPENMP_FOUND)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

find_package(ZLIB REQUIRED)

find_package(CGAL)

# Get Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    # message(STATUS "${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# add dependencies
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(AUTO_DOWNLOAD "Auto download dependencies" ON)
include(rpdDependencies)
add_subdirectory(extern/mshloader)
if(USE_GPU)
    add_subdirectory(src/dist2mat)
    add_subdirectory(src/rpd3d)
    add_subdirectory(src/rpd3d_api)
    add_subdirectory(src/matfun_fix)
    add_subdirectory(src/IO/IO_CUDA)
endif(USE_GPU)
add_subdirectory(src/inputs)
add_subdirectory(src/matbase)
add_subdirectory(src/matfun)
add_subdirectory(src/IO/IO_CXX)

# ###############################################################################
# For VolumeVoronoiGPU binary
# ###############################################################################
set(SOURCE_LIST
    "src/main.cpp"
)

add_executable(${PROJECT_NAME} ${SOURCE_LIST})
if(USE_GPU)
    target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        --generate-line-info
        --use_fast_math
        >)
endif(USE_GPU)
target_compile_definitions(${PROJECT_NAME} PUBLIC -Dgeogram_EXPORTS)

# add common_headers include
target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

if(${MSVC})
    if(${CMAKE_VERSION} VERSION_GREATER 3.6.0)
        set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "${PROJECT_NAME}")
    endif()
endif()

target_link_directories(${PROJECT_NAME} PUBLIC "${CMAKE_BINARY_DIR}/lib")
if(USE_GPU)
    target_link_libraries(${PROJECT_NAME}
        general cuda
        general cublas
        optimized rpd3d
        optimized rpd3d_api
        optimized dist2mat
        optimized matfun_fix
        optimized IO_CUDA
    )
    target_link_libraries(${PROJECT_NAME} general ${CUDA_LIBRARIES})
endif(USE_GPU)
target_link_libraries(${PROJECT_NAME}
    optimized IO_CXX
    optimized inputs
    optimized matbase
    optimized matfun
    # optimized cpu_fallback
    optimized mshloader
    optimized polyscope
    optimized geogram
    optimized nlohmann_json::nlohmann_json
    optimized ZLIB::ZLIB
    optimized CGAL::CGAL
    optimized Eigen3::Eigen
    optimized igl::core
)

# ###############################################################################
# For VORO_GUI binary
# ###############################################################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pthread -D_DEBUG")

add_executable(VORO_GUI src/test_gui.cpp)
target_link_libraries(VORO_GUI PUBLIC polyscope mshloader)