# raptor search FPGA {#usage_search_fpga}

<!--
SPDX-FileCopyrightText: 2025 Tobias Baumann & Zuse Institute Berlin
SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

## FPGA-only Parameters

\note
This Document only lists additional parameters for the FPGA variant of `raptor search`.

### --fpga
Use the FPGA.

### --buffer
The size (in MiB) of the host side double buffer to use.

### -​-kernels
The number of kernel copys to use. Note that this parameter has to correspond to the available bitstreams created during compilation.
