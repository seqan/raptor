# --------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

project (docs LANGUAGES NONE VERSION "2.0.1")

set (RAPTOR_DOXYGEN_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component (RAPTOR_CLONE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)

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

add_custom_target(doc ALL
                  COMMAND ${DOXYGEN_EXECUTABLE}
                  WORKING_DIRECTORY ${RAPTOR_DOXYGEN_OUTPUT_DIR}
                  BYPRODUCTS html/
                  DEPENDS download-cppreference-doxygen-web-tag
                  COMMENT "Generating developer API documentation with Doxygen"
                  VERBATIM)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (
    DIRECTORY "${RAPTOR_DOXYGEN_OUTPUT_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_DOCDIR}"
    COMPONENT doc
    EXCLUDE_FROM_ALL
    PATTERN "*.md5" EXCLUDE
    PATTERN "*.map" EXCLUDE
    PATTERN "formula.repository" EXCLUDE
)

### Enable testing

enable_testing()

message (STATUS "Add doc test.")
add_test(NAME doc_no_warnings_test
         COMMAND bash -c "${RAPTOR_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
