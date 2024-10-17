# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Minimum cmake version
cmake_minimum_required (VERSION 3.17...3.30)

set (RAPTOR_DOXYGEN_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component (RAPTOR_CLONE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)

include (${RAPTOR_CLONE_DIR}/cmake/version.cmake)

set (CPM_INDENT "CMake Package Manager CPM: ")
include (${RAPTOR_CLONE_DIR}/cmake/CPM.cmake)
CPMUsePackageLock (${RAPTOR_CLONE_DIR}/cmake/package-lock.cmake)
CPMGetPackage (doxygen_awesome)

project (docs
         LANGUAGES NONE
         VERSION "${RAPTOR_PROJECT_VERSION}"
)

include (raptor-doxygen.cmake)

enable_testing ()

message (STATUS "Configuring doc.")

set (RAPTOR_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (RAPTOR_DOXYGEN_EXCLUDE_SYMBOLS "")
set (RAPTOR_DOXYGEN_PREDEFINED_NDEBUG "")
set (RAPTOR_DOXYGEN_ENABLED_SECTIONS "DEV")
set (RAPTOR_SHOW_DEV_GROUP "YES")
set (RAPTOR_DOXYGEN_EXTRACT_PRIVATE "YES")

configure_file (${RAPTOR_DOXYFILE_IN} ${RAPTOR_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${RAPTOR_FOOTER_HTML_IN} ${RAPTOR_DOXYGEN_OUTPUT_DIR}/raptor_footer.html)
include (raptor-doxygen-layout.cmake)

add_custom_target (doc ALL
                   COMMAND ${DOXYGEN_EXECUTABLE}
                   WORKING_DIRECTORY ${RAPTOR_DOXYGEN_OUTPUT_DIR}
                   BYPRODUCTS html/
                   DEPENDS download-cppreference-doxygen-web-tag
                   COMMENT "Generating developer API documentation with Doxygen"
                   VERBATIM
)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (DIRECTORY "${RAPTOR_DOXYGEN_OUTPUT_DIR}/html"
         DESTINATION "${CMAKE_INSTALL_DOCDIR}"
         COMPONENT doc
         EXCLUDE_FROM_ALL
         PATTERN "*.md5" EXCLUDE
         PATTERN "*.map" EXCLUDE
         PATTERN "formula.repository" EXCLUDE
)

### Enable testing

enable_testing ()

message (STATUS "Add doc test.")
add_test (NAME doc_no_warnings_test
          COMMAND bash -c "${RAPTOR_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties (doc_no_warnings_test PROPERTIES DEPENDS cppreference-doxygen-web-tag)
