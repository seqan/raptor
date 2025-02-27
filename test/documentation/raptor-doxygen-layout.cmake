# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.17...3.30)

include (${CMAKE_CURRENT_LIST_DIR}/../cmake/raptor_test_files.cmake)

# Replaces documentation entries in variable `DOXYGEN_LAYOUT`
#
# Important: the variable of the parent scope must be named `DOXYGEN_LAYOUT` for this function to work properly.
#
# Example:
# replace_in_doxygen_layout ("${RAPTOR_INCLUDE_DIR}/../doc/tutorial/" "Tutorial")
#
# (1) Find all 'index.md' files in ${RAPTOR_INCLUDE_DIR}/../doc/tutorial/
# (2) For each 'index.md'
# (a) Use the very first line, e.g. `# FOOOMOO {#fooooo_bar_ref}`,
# for the title, e.g. "FOOOMOO" and doxygen reference, e.g. "fooooo_bar_ref".
# (b) Generate doxygen html layout entry, e.g.
# "      <tab type=\"user\" visible=\"yes\" title=\"FOOOMOO" url=\"\\\\ref fooooo_bar_ref\" intro=\"\"/>\n"
# and append it to a list
# (3) Replace the doxygen html layout entry list from (2) in the ${DOXYGEN_LAYOUT} input variable
#
function (replace_in_doxygen_layout doc_path doxygen_layout_tag)
    set (DOXYGEN_LAYOUT_TAG_LINE
         "<tab type=\"usergroup\" visible=\"yes\" title=\"${doxygen_layout_tag}\" intro=\"\">\n"
    )
    set (DOXYGEN_LAYOUT_DOC_PAGES ${DOXYGEN_LAYOUT_TAG_LINE}) # append header line

    # iterate over all index.md
    raptor_test_files (doc_how_to_filenames "${doc_path}" "index.md")

    foreach (doc_how_to_filename ${doc_how_to_filenames})
        set (doc_howto_filepath "${doc_path}/${doc_how_to_filename}")
        execute_process (COMMAND head -n 1 ${doc_howto_filepath} OUTPUT_VARIABLE DOC_HEADER_LINE)
        string (REGEX MATCH "^# \(.*\) {#\(.*\)}" DUMMY ${DOC_HEADER_LINE})
        set (doc_title ${CMAKE_MATCH_1})
        set (doc_ref_name ${CMAKE_MATCH_2})
        string (APPEND
                DOXYGEN_LAYOUT_DOC_PAGES
                "      <tab type=\"user\" visible=\"yes\" title=\"${doc_title}\" url=\"\\\\ref ${doc_ref_name}\" intro=\"\"/>\n"
        )

        unset (doc_howto_filepath)
        unset (doc_title)
        unset (doc_ref_name)
    endforeach ()

    # Replace header line and appended list of doc entries with header line
    string (REGEX REPLACE "${DOXYGEN_LAYOUT_TAG_LINE}" "${DOXYGEN_LAYOUT_DOC_PAGES}" NEW_DOXYGEN_LAYOUT
                          ${DOXYGEN_LAYOUT}
    )

    set (DOXYGEN_LAYOUT
         ${NEW_DOXYGEN_LAYOUT}
         PARENT_SCOPE
    ) # replace new Doxygen layout

    unset (DOXYGEN_LAYOUT_TAG_LINE)
    unset (DOXYGEN_LAYOUT_DOC_PAGES)
endfunction ()

# ## Add all documentation pages to DoxygenLayout.xml.in
# ## ---------------------------------------------------
# Note: variable name DOXYGEN_LAYOUT must not be changed because it is directly used within `replace_in_doxygen_layout`
file (READ "${RAPTOR_LAYOUT_IN}" DOXYGEN_LAYOUT)

replace_in_doxygen_layout ("${RAPTOR_DOXYGEN_SOURCE_DIR}/doc/usage/" "Usage")
replace_in_doxygen_layout ("${RAPTOR_DOXYGEN_SOURCE_DIR}/doc/methods/" "Methods")

file (WRITE "${RAPTOR_DOXYGEN_OUTPUT_DIR}/DoxygenLayout.xml" ${DOXYGEN_LAYOUT})
