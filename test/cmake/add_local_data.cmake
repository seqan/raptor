# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

set (DATASOURCES_DATA_DIR "${Raptor_SOURCE_DIR}/test/data")

file (GLOB_RECURSE datasources
      LIST_DIRECTORIES false
      RELATIVE ${DATASOURCES_DATA_DIR}
      CONFIGURE_DEPENDS ${DATASOURCES_DATA_DIR}/*
)
list (REMOVE_ITEM datasources datasources.cmake README.md REUSE.toml)

foreach (datasource IN LISTS datasources)
    get_filename_component (datasource_name "${datasource}" NAME)
    set (data_dir "${CMAKE_CURRENT_BINARY_DIR}/data")
    file (MAKE_DIRECTORY "${data_dir}")

    if (datasource_name MATCHES ".*\.layout$")
        configure_file ("${DATASOURCES_DATA_DIR}/${datasource}" "${data_dir}/${datasource_name}")
    else ()
        configure_file ("${DATASOURCES_DATA_DIR}/${datasource}" "${data_dir}/${datasource_name}" COPYONLY)
    endif ()
endforeach ()
