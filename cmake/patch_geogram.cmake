# ###############################################################################
# Patches the downloaded geogram 1.7.9 sources for current Apple toolchains
# (FetchContent PATCH_COMMAND, run with -DGEOGRAM_DIR=<geogram sources>).
#   - Darwin: drop -msse3, which clang >= 17 rejects when targeting arm64
#     (x86_64 macOS targets have SSE3 anyway)
#   - zlib: do not redefine fdopen() on macOS, which breaks the stdio.h of
#     the macOS 15 SDK (the same change as zlib 1.2.12)
# Patching twice is a no-op.
# ###############################################################################

function(libmat_patch_file file from to)
	file(READ "${GEOGRAM_DIR}/${file}" text)
	string(REPLACE "${from}" "${to}" patched "${text}")
	if(NOT patched STREQUAL text)
		file(WRITE "${GEOGRAM_DIR}/${file}" "${patched}")
	endif()
endfunction()

libmat_patch_file(cmake/platforms/Darwin-clang.cmake
	"# Enable SSE3 instruction set\nadd_flags(CMAKE_CXX_FLAGS -msse3)\nadd_flags(CMAKE_C_FLAGS -msse3)\n" "")
libmat_patch_file(src/lib/geogram/third_party/zlib/zutil.h
	"#if defined(MACOS) || defined(TARGET_OS_MAC)" "#if defined(MACOS)")
