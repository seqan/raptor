cmake_minimum_required (VERSION 3.8)

project (fasta_to_fastq VERSION 1.0.0)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE
    )
endif ()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Dependency: SeqAn3.
set (SEQAN3_SUBMODULES_DIR "../../../../lib")
find_package (SeqAn3 QUIET REQUIRED HINTS ../../../../lib/seqan3/build_system)

include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
seqan3_require_ccache ()

add_subdirectory (src)
