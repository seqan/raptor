# Setup {#setup}

<!--
SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

# Download and Installation
There may be performance benefits when compiling from source as the build can be optimized for the host system.

## Install with conda and bioconda (Linux)
```bash
conda install -c bioconda -c conda-forge raptor
```

## Compile from source

### Prerequisites
* CMake >= 3.21
* GCC 12, 13 or 14 (most recent minor version)
* Clang 17 or 18 (most recent minor version)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in depth
information.

### Download current main branch
```bash
git clone https://github.com/seqan/raptor
cd raptor
git submodule update --init
```

### Download specific version
E.g., for version `1.1.0`:
```bash
git clone --branch raptor-v1.1.0 --recurse-submodules https://github.com/seqan/raptor
```
Or from within an existing repository
```bash
git checkout raptor-v1.1.0
```

### Building
```bash
cd raptor
mkdir -p build
cd build
cmake ..
make
```

The binary can be found in `bin`.

You may want to add the Raptor executable to your PATH:
```
export PATH=$(pwd)/bin:$PATH
raptor --version
```

By default, Raptor will be built with host specific optimizations (`-march=native`). This behavior can be disabled by
passing `-DHIBF_NATIVE_BUILD=OFF` to CMake.
