# ###############################################################################
# Build geogram (from LIBMAT_GEOGRAM_DIR) as a static library, so that the
# python extension module is self-contained.
# ###############################################################################

if(TARGET geogram)
	return()
endif()

set(GEOGRAM_SOURCE_INCLUDE_DIR "${LIBMAT_GEOGRAM_DIR}/src/lib")

# ###############################################################################
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set(VORPALINE_ARCH_64 TRUE CACHE BOOL "" FORCE)
	set(VORPALINE_PLATFORM Win-vs-generic CACHE STRING "" FORCE)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
		# geogram 1.7.x uses x86 assembly for its atomics on any Linux
		message(FATAL_ERROR "geogram 1.7.9 does not build on Linux ARM: "
			"libmat supports Linux x86_64 and macOS (Intel and Apple silicon)")
	endif()
	set(VORPALINE_PLATFORM Linux64-gcc CACHE STRING "" FORCE)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(VORPALINE_PLATFORM Darwin-clang CACHE STRING "" FORCE)
endif()

option(GEOGRAM_WITH_GRAPHICS "Viewers and geogram_gfx library" OFF)
option(GEOGRAM_WITH_LEGACY_NUMERICS "Legacy numerical libraries" OFF)
option(GEOGRAM_WITH_HLBFGS "Non-linear solver (Yang Liu's HLBFGS)" OFF)
option(GEOGRAM_WITH_TETGEN "Tetrahedral mesher (Hang Si's TetGen)" OFF)
option(GEOGRAM_WITH_TRIANGLE "Triangle mesher (Jonathan Shewchuk's triangle)" OFF)
option(GEOGRAM_WITH_EXPLORAGRAM "Experimental code (hexahedral meshing vpipeline and optimal transport)" OFF)
option(GEOGRAM_WITH_LUA "Built-in LUA interpreter" OFF)
option(GEOGRAM_LIB_ONLY "Libraries only (no example programs/no viewer)" ON)
option(GEOGRAM_WITH_FPG "Predicate generator (Sylvain Pion's FPG)" OFF)
option(GEOGRAM_USE_SYSTEM_GLFW3 "Use the version of GLFW3 installed in the system if found" OFF)

# ###############################################################################
# geogram asks for CMake 2.8, which CMake >= 4 refuses unless told otherwise
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
	set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
	set(LIBMAT_UNSET_POLICY_VERSION_MINIMUM TRUE)
endif()
add_subdirectory(${LIBMAT_GEOGRAM_DIR} ${CMAKE_BINARY_DIR}/geogram EXCLUDE_FROM_ALL)
if(LIBMAT_UNSET_POLICY_VERSION_MINIMUM)
	unset(CMAKE_POLICY_VERSION_MINIMUM)
endif()
target_include_directories(geogram SYSTEM PUBLIC ${GEOGRAM_SOURCE_INCLUDE_DIR})
add_library(geogram::geogram ALIAS geogram)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	# geogram is compiled with -fopenmp on Linux, so whatever links the static
	# library also needs the OpenMP runtime (plain signature, as geogram uses)
	find_package(OpenMP REQUIRED COMPONENTS CXX)
	target_link_libraries(geogram OpenMP::OpenMP_CXX)
endif()

# ###############################################################################
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	# remove warning for multiply defined symbols (caused by multiple
	# instanciations of STL templates)
	target_compile_options(geogram INTERFACE /wd4251)

	# remove all unused stuff from windows.h
	target_compile_definitions(geogram INTERFACE -DWIN32_LEAN_AND_MEAN)
	target_compile_definitions(geogram INTERFACE -DVC_EXTRALEAN)

	# do not define a min() and a max() macro, breaks
	# std::min() and std::max() !!
	target_compile_definitions(geogram INTERFACE -DNOMINMAX)

	# we want M_PI etc...
	target_compile_definitions(geogram INTERFACE -D_USE_MATH_DEFINES)
endif()
