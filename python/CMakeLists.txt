# ###############################################################################
# Python bindings: builds the extension module libmat._libmat
#   - pip install .                      (scikit-build-core, see pyproject.toml)
#   - cmake -DLIBMAT_BUILD_PYTHON=ON     (in-place: python/libmat/_libmat.*.so,
#                                          then `export PYTHONPATH=<repo>/python`)
# ###############################################################################
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 2.12 CONFIG QUIET)
if(NOT pybind11_FOUND)
    include(FetchContent)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.13.6
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(pybind11)
endif()

file(GLOB LIBMAT_BINDING_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/*.cpp)
pybind11_add_module(_libmat MODULE ${LIBMAT_BINDING_SOURCES})
target_link_libraries(_libmat PRIVATE libmat::libmat)
target_compile_definitions(_libmat PRIVATE LIBMAT_VERSION_INFO="${PROJECT_VERSION}")

if(SKBUILD)
    install(TARGETS _libmat LIBRARY DESTINATION libmat)
else()
    # in-place build for development
    set_target_properties(_libmat PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/libmat"
        LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/libmat"
        LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/libmat")
endif()
