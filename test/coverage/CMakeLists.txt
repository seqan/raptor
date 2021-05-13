cmake_minimum_required (VERSION 3.8)

# Add a custom build type: Coverage

set (CMAKE_CXX_FLAGS_COVERAGE
     "${CMAKE_CXX_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage" CACHE STRING
     "Flags used by the C++ compiler during coverage builds."
     FORCE)
set (CMAKE_C_FLAGS_COVERAGE
     "${CMAKE_C_FLAGS_DEBUG} --coverage -fprofile-arcs -ftest-coverage" CACHE STRING
     "Flags used by the C compiler during coverage builds."
     FORCE)
set (CMAKE_EXE_LINKER_FLAGS_COVERAGE
     "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -Wl,-lgcov" CACHE STRING
     "Flags used for linking binaries during coverage builds."
     FORCE)
set (CMAKE_SHARED_LINKER_FLAGS_COVERAGE
     "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -Wl,-lgcov" CACHE STRING
     "Flags used by the shared libraries linker during coverage builds."
     FORCE)

mark_as_advanced (
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE)
