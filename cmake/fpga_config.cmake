# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

if (TARGET raptor::fpga::interface)
    return ()
endif ()

# Effect: -DRAPTOR_FPGA=ON/-DRAPTOR_FPGA_HARDWARE=ON is only honored if the IntelLLVM compiler is used
cmake_dependent_option (RAPTOR_FPGA "Enable FPGA support." OFF "${CMAKE_CXX_COMPILER_ID} STREQUAL IntelLLVM" OFF)
cmake_dependent_option (RAPTOR_FPGA_HARDWARE "Compile for FPGA hardware." OFF
                        "${CMAKE_CXX_COMPILER_ID} STREQUAL IntelLLVM" OFF
)

if (RAPTOR_FPGA_HARDWARE)
    set (RAPTOR_FPGA ON)
endif ()

if (NOT RAPTOR_FPGA)
    target_compile_definitions (raptor_interface INTERFACE RAPTOR_FPGA=0)
    raptor_config_print ("FPGA support:               disabled")
    return ()
endif ()

if (RAPTOR_FPGA_HARDWARE)
    raptor_config_print ("FPGA support:               Hardware")
else ()
    raptor_config_print ("FPGA support:               Emulation")
endif ()

target_compile_definitions (raptor_interface INTERFACE RAPTOR_FPGA=1)

# CPM ignores DOWNLOAD_ONLY when package overrides are used
if (DEFINED CPM_ibf-fpga_SOURCE)
    set (ibf-fpga_SOURCE_DIR "${CPM_ibf-fpga_SOURCE}")
else ()
    CPMGetPackage (ibf-fpga)
endif ()

# Shared Includes between all targets
add_library (raptor_fpga_interface INTERFACE)
target_link_libraries (raptor_fpga_interface INTERFACE "raptor::interface")
target_include_directories (raptor_fpga_interface SYSTEM INTERFACE "${ibf-fpga_SOURCE_DIR}/include")
target_compile_options (raptor_fpga_interface INTERFACE "-fsycl" "-fintelfpga" "-Xshyper-optimized-handshaking=off"
                                                        "-qactypes"
)
target_link_options (raptor_fpga_interface INTERFACE "-fsycl" "-fintelfpga" "-Xshyper-optimized-handshaking=off"
                     "-qactypes"
)
add_library (raptor::fpga::interface ALIAS raptor_fpga_interface)

if (RAPTOR_FPGA_HARDWARE)
    # Try to detect FPGA board by searching for known boards in the output of "aoc -list-boards"
    if (NOT DEFINED FPGA_DEVICE)
        set (known_boards "ofs_ia840f_usm" "ofs_d5005_usm")

        execute_process (COMMAND aoc -list-boards
                         OUTPUT_VARIABLE output
                         OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
        )

        foreach (search_string IN LISTS known_boards)
            string (FIND "${output}" "${search_string}" string_pos)

            # If string is found
            if (NOT string_pos EQUAL -1)
                raptor_config_print ("Found FPGA board '${search_string}' in the output of 'aoc -list-boards'")
                set (FPGA_DEVICE "${search_string}")
                break ()
            endif ()
        endforeach ()
    endif ()

    # Fall back to pre defined device for FPGA board selection
    if (NOT DEFINED FPGA_DEVICE)
        set (FPGA_DEVICE "intel_s10sx_pac:pac_s10_usm")
        message (STATUS "    FPGA_DEVICE was not specified.\n"
                        "         Configuring the design to run on the default FPGA board ${FPGA_DEVICE}.\n"
                        "         Please refer to the README for information on board selection."
        )
    else ()
        raptor_config_print ("Configuring the design to run on FPGA board ${FPGA_DEVICE}")
    endif ()
    # use cmake -D USER_HARDWARE_FLAGS=<flags> to set extra flags for FPGA backend compilation
    target_link_options (raptor_fpga_interface INTERFACE -Xshardware -Xstarget=${FPGA_DEVICE} ${USER_HARDWARE_FLAGS})
endif ()

if (NOT DEFINED WINDOW_SIZE_LIST)
    set (WINDOW_SIZE_LIST 23)
    raptor_config_print ("  No WINDOW_SIZE_LIST supplied. Defaulting to '${WINDOW_SIZE_LIST}'.")
endif ()
list (JOIN WINDOW_SIZE_LIST "," WINDOW_SIZE_STRING)

if (NOT DEFINED MIN_IBF_K_LIST)
    set (MIN_IBF_K_LIST 19)
    raptor_config_print ("  No MIN_IBF_K_LIST supplied. Defaulting to '${MIN_IBF_K_LIST}'.")
endif ()
list (JOIN MIN_IBF_K_LIST "," MIN_IBF_K_STRING)

if (NOT DEFINED BIN_COUNT_LIST)
    set (BIN_COUNT_LIST 64 128)
    raptor_config_print ("  No BIN_COUNT_LIST supplied. Defaulting to '${BIN_COUNT_LIST}'.")
endif ()
list (JOIN BIN_COUNT_LIST "," BIN_COUNT_STRING)

if (NOT DEFINED KERNEL_COPYS_LIST)
    set (KERNEL_COPYS_LIST 1 2)
    raptor_config_print ("  No KERNEL_COPYS_LIST supplied. Defaulting to '${KERNEL_COPYS_LIST}'.")
endif ()
list (JOIN KERNEL_COPYS_LIST "," KERNEL_COPYS_STRING)
