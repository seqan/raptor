# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# SYSTEM property is only implemented with CMake 3.25
macro (CPMGetSystemPackage package)
    CPMGetPackage (${package})

    if (CMAKE_VERSION VERSION_LESS 3.25)
        set (target_name "${package}")
        if ("${package}" STREQUAL "seqan3")
            set (target_name "seqan3_seqan3")
        elseif ("${package}" STREQUAL "sharg")
            set (target_name "sharg_sharg")
        endif ()

        if (${package}_ADDED)
            set (interface_include "$<TARGET_PROPERTY:${target_name},INTERFACE_INCLUDE_DIRECTORIES>")
            set (include "$<TARGET_PROPERTY:${target_name},INCLUDE_DIRECTORIES>")
            set_target_properties (${target_name}
                                   PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                                              "$<$<BOOL:${interface_include}>:${interface_include}>;$<$<BOOL:${include}>:${include}>"
            )
            unset (interface_include)
            unset (include)
        endif ()

        unset (target_name)
    endif ()
endmacro ()
