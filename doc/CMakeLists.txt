cmake_minimum_required (VERSION 3.8)

# Find doxygen.
find_package (Doxygen QUIET)

if (${DOXYGEN_FOUND})
    message (STATUS "Found Doxygen: ${DOXYGEN_EXECUTABLE}")

    # Configure doxygen options.
    set (APP_TEMPLATE_DOXYFILE_IN ${CMAKE_SOURCE_DIR}/doc/doxygen_cfg)
    set (APP_TEMPLATE_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/doc")
    set (APP_TEMPLATE_DOXYGEN_EXCLUDE_SYMBOLS "")
    set (APP_TEMPLATE_DOXYGEN_PREDEFINED_NDEBUG "")
    set (APP_TEMPLATE_DOXYGEN_ENABLED_SECTIONS "DEV")
    set (APP_TEMPLATE_DOXYGEN_EXTRACT_PRIVATE "YES")

    configure_file (${APP_TEMPLATE_DOXYFILE_IN} ${APP_TEMPLATE_DOXYGEN_OUTPUT_DIR}/Doxyfile)
    add_custom_target (doc
                       COMMAND ${DOXYGEN_EXECUTABLE}
                       WORKING_DIRECTORY ${APP_TEMPLATE_DOXYGEN_OUTPUT_DIR}
                       COMMENT "Generating (developer) API documentation with Doxygen."
                       VERBATIM)
    message (STATUS "${FontBold}You can run `make doc` to build api documentation.${FontReset}")
else ()
    message (STATUS "Doxygen not found.")
endif ()
