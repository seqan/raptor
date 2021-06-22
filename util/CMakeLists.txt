cmake_minimum_required (VERSION 3.8)

# Define the application name.
project (raptor_utility)

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS ../lib/seqan3/build_system)

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Use ccache.
include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
seqan3_require_ccache ()

# Optionally build mason and DREAM-Yara.
option (RAPTOR_UTILITY_BUILD_MASON "Build mason. Needed to generate artificial data set." OFF)
option (RAPTOR_UTILITY_BUILD_DREAM_YARA "Build DREAM-Yara. Needed to map reads to the NCBI RefSeq." OFF)

if (RAPTOR_UTILITY_BUILD_MASON)
    message ("-- Configuring Mason")
    include (cmake/raptor_require_mason.cmake)
    raptor_require_mason()
    message ("-- Configuring Mason - done")
else ()
    message ("-- Mason will not be built. Enable via 'cmake -DRAPTOR_UTILITY_BUILD_MASON=1 .'")
endif ()

if (RAPTOR_UTILITY_BUILD_DREAM_YARA)
    message ("-- Configuring DREAM-Yara")
    include (cmake/raptor_require_dream_yara.cmake)
    raptor_require_dream_yara()
    message ("-- Configuring DREAM-Yara - done")
else ()
    message ("-- DREAM-Yara will not be built. Enable via 'cmake -DRAPTOR_UTILITY_BUILD_DREAM_YARA=1 .'")
endif()

add_library (common INTERFACE)
target_link_libraries ("common" INTERFACE seqan3::seqan3)
target_include_directories ("common" INTERFACE ../include)

add_executable ("generate_reads" src/applications/generate_reads.cpp)
target_link_libraries ("generate_reads" "common")
install (TARGETS generate_reads DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("split_sequence" src/applications/split_sequence.cpp)
target_link_libraries ("split_sequence" "common")
install (TARGETS split_sequence DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("count_minimiser" src/applications/count_minimiser.cpp)
target_link_libraries ("count_minimiser" "common")
install (TARGETS count_minimiser DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("apply_taxsbp" src/applications/apply_taxsbp.cpp)
target_link_libraries ("apply_taxsbp" "common")
install (TARGETS apply_taxsbp DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_executable ("generate_reads_refseq" src/applications/generate_reads_refseq.cpp)
target_link_libraries ("generate_reads_refseq" "common")
install (TARGETS generate_reads_refseq DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

message (STATUS "${FontBold}You can run `make install` to build the application.${FontReset}")
