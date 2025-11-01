# Search with Raptor on FPGAs {#tutorial_search_fpga}

<!--
SPDX-FileCopyrightText: 2025 Tobias Baumann & Zuse Institute Berlin
SPDX-License-Identifier: CC-BY-4.0
-->

You will learn how to use an FPGA to search a Raptor index.

\tutorial_head{Easy, 30 min, \ref tutorial_search, }

[TOC]

# Use a FPGA for searching an index

\note
The support for FPGAs in Raptor is considered EXPERIMENTAL.

## Perquisites

- a oneAPI compatible FPGA supporting Unified Shared Memory (USM)
- oneAPI <= 2025.0 (Intel dropped support for FPGAs in the oneAPI 2025.1 release)
- Quartus (recommended version depends on your FPGA and used oneAPI version)
- oneAPI compatible Board Support Package (BSP) and/or Card Support Package (CSP) for your FPGA supporting Unified Shared Memory (USM)

## Building

FPGA support can run in either emulation (one the CPU) or hardware mode. While emulation is good for verifying functional correctness the FPGA will be only be used in hardware mode. Compiling for the FPGA (hardware mode) will take multiple hours while compiling for the emulator is usually done in minutes.

To enable FPGA support in Raptor use `icpx` as compiler and supply `-DRAPTOR_FPGA=ON` (emulation) or `-DRAPTOR_FPGA_HARDWARE=ON` (hardware) to the CMake command line when building Raptor from source:

```bash
CXX=icpx cmake -DRAPTOR_FPGA=ON ..
CXX=icpx cmake -DRAPTOR_FPGA_HARDWARE=ON ..
```
CMake will try to detect which FPGA is installed in your system by parsing the output of `aoc -list-boards`. This only works if your FPGA is in a list of known boards. You can specify your board manually by supplying `-DFPGA_DEVICE=<your_fpga>` to the CMake command line. Note that Unified Shared Memory (USM) is required and that the USM base image of your board needs to be specified.

Additionally, you need to specify at compile time which exact parameters (`w`, `k`, `number of bins`) you like to use with the FPGA, otherwise defaults will be used.
You can specify these parameters by supplying `WINDOW_SIZE_LIST`, `MIN_IBF_K_LIST` and `BIN_COUNT_LIST`. During the build process, all possible combinations will be build.

Specific to FPGAs, you can supply a list of factors by which the logic on the FPGA will be replicated via `KERNEL_COPYS_LIST`. During execution, you will use `--kernels` to select on of the factors supplied via `KERNEL_COPYS_LIST` to select the right hardware bitstream generated during the build process

A full CMake command can look like this:

```bash
CXX=icpx cmake -DRAPTOR_FPGA_HARDWARE=ON -DFPGA_DEVICE=ia840f:ofs_ia840f_usm -DWINDOW_SIZE_LIST="23;24" -DMIN_IBF_K_LIST="19;20" -DBIN_COUNT_LIST=8192 -DKERNEL_COPYS_LIST="3;4" ..
```

When compiling for FPGA (hardware), sometimes the build process will fail or produce unusable bitstreams. In this case try lower factors for `KERNEL_COPYS_LIST` or supply another seed via `-DUSER_HARDWARE_FLAGS=-Xsseed=<number>`, e.g.:

```bash
CXX=icpx cmake -DRAPTOR_FPGA_HARDWARE=ON -DUSER_HARDWARE_FLAGS=-Xsseed=42 ..
```

## Execution

Make sure your hardware and software setup is correct and that your FPGA is initialized with the USM variant of your oneAPI-compatible base image. This is usually done by running `aocl initialize acl0 ofs_<fpga model>_usm` but can vary depending on your FPGA.

To enable execution on the FPGA, supply `--fpga` to `raptor search`. Additionally, you can provide a buffer size (in this example 10 MiB) via `--buffer 10` and select a bitstream with replicated logic (in this example 4) via `--kernels 4`.

Note that a pre-compiled hardware bitstream has to exist for each parameter (`w`, `k`, `number of bins`, `number of kernel copys`) combination you try to run the search with. If you try to run a parameter combination without an existing pre-compiled bitstream, execution will fail. See the above section of this document on how to pre-compile hardware bitstreams.
