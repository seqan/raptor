# --------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.17)

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_DOT_FOUND})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (RAPTOR_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (RAPTOR_DOXYGEN_HAVE_DOT "YES")
endif ()

### Use mathjax instead of latex to render formulas.
set (RAPTOR_DOXYGEN_USE_MATHJAX "NO")

### Number of threads to use for dot. Doxygen's default is 0 (all threads).
set (RAPTOR_DOXYGEN_DOT_NUM_THREADS "0")

### Configure doc/developer targets.
set (RAPTOR_DOXYGEN_SOURCE_DIR "${RAPTOR_CLONE_DIR}")
set (RAPTOR_DOXYFILE_IN ${RAPTOR_DOXYGEN_INPUT_DIR}/raptor_doxygen_cfg.in)
set (RAPTOR_FOOTER_HTML_IN ${RAPTOR_DOXYGEN_INPUT_DIR}/raptor_footer.html.in)
set (RAPTOR_LAYOUT_IN ${RAPTOR_DOXYGEN_INPUT_DIR}/DoxygenLayout.xml)

option (RAPTOR_VERCEL_PREVIEW_DOC "Is this a preview build by vercel.com?" OFF)

set (RAPTOR_DOXYGEN_EXTERNAL_SEARCH "NO")

if (RAPTOR_VERCEL_PREVIEW_DOC)
    set (RAPTOR_DOXYGEN_EXTERNAL_SEARCH "YES")
    set (RAPTOR_DOXYGEN_USE_MATHJAX "YES")
    set (RAPTOR_DOXYGEN_DOT_NUM_THREADS "2")
    set (RAPTOR_DOXYFILE_OPTION_POWERED_BY_VERCEL
         "HTML_EXTRA_FILES       += ${RAPTOR_DOXYGEN_SOURCE_DIR}/test/documentation/.vercel/powered-by-vercel.svg"
    )
    set (RAPTOR_FOOTER_HTML_OPTION_POWERED_BY_VERCEL
         "<li class='footer'><a href='https://vercel.com/?utm_source=seqan&utm_campaign=oss'><img src='$relpath^powered-by-vercel.svg' height='29px' alt='Powered by Vercel'/></a></li>"
    )
endif ()

### Download and extract cppreference-doxygen-web.tag.xml for std:: documentation links
set (RAPTOR_DOXYGEN_STD_TAGFILE "${PROJECT_BINARY_DIR}/cppreference-doxygen-web.tag.xml")
include (ExternalProject)
ExternalProject_Add (download-cppreference-doxygen-web-tag
                     URL "https://github.com/PeterFeicht/cppreference-doc/releases/download/v20220730/html-book-20220730.tar.xz"
                     URL_HASH SHA256=71f15003c168b8dc5a00cbaf19b6480a9b3e87ab7e462aa39edb63d7511c028b
                     TLS_VERIFY ON
                     DOWNLOAD_DIR "${PROJECT_BINARY_DIR}"
                     DOWNLOAD_NAME "html-book.tar.xz"
                     DOWNLOAD_NO_EXTRACT YES
                     BINARY_DIR "${PROJECT_BINARY_DIR}"
                     CONFIGURE_COMMAND /bin/sh -c "xzcat html-book.tar.xz | tar -xf - cppreference-doxygen-web.tag.xml"
                     BUILD_COMMAND rm "html-book.tar.xz"
                     INSTALL_COMMAND ""
)

### TEST HELPER

# doxygen does not show any warnings (doxygen prints warnings / errors to cerr)
set (RAPTOR_TEST_DOXYGEN_FAIL_ON_WARNINGS
     "
    ${DOXYGEN_EXECUTABLE} > doxygen.cout 2> doxygen.cerr;
    cat \"doxygen.cerr\";
    test ! -s \"doxygen.cerr\""
)

### install helper

# make sure that prefix path is /usr/local/share/doc/raptor/
if (NOT DEFINED CMAKE_SIZEOF_VOID_P)
    # we need this to suppress GNUInstallDirs AUTHOR_WARNING:
    #   CMake Warning (dev) at /usr/share/cmake-3.19/Modules/GNUInstallDirs.cmake:223 (message):
    #     Unable to determine default CMAKE_INSTALL_LIBDIR directory because no
    #     target architecture is known.  Please enable at least one language before
    #     including GNUInstallDirs.
    set (CMAKE_SIZEOF_VOID_P 8)
endif ()
include (GNUInstallDirs) # this is needed to prefix the install paths
