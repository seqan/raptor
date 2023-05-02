# Quick Setup {#setup}

<b>Learning Objective:</b><br>
In this short guide you will learn how to set up Raptor or how to compile a specific version.

\tutorial_head{Easy, 30 Minutes, No prerequisites, }

[TOC]

# Download and Installation
There may be performance benefits when compiling from source as the build can be optimized for the host system.

## Install with conda and bioconda (Linux)

```bash
conda install -c bioconda -c conda-forge raptor
```

## Install with Homebrew (Linux, macOS)

```bash
brew install brewsci/bio/raptor
```

## Compile from source

\startcollapsible{Prerequisites (click to expand)}

* CMake >= 3.18
* GCC 11, 12 or 13 (most recent minor version)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in depth
information.

\endcollapsible

\startcollapsible{Download current main branch (click to expand)}

```bash
git clone https://github.com/seqan/raptor
git submodule update --init
```

\endcollapsible

\startcollapsible{Download specific version (click to expand)}

E.g., for version `1.1.0`:
```bash
git clone --branch raptor-v1.1.0 --recurse-submodules https://github.com/seqan/raptor
```
Or from within an existing repository
```bash
git checkout raptor-v1.1.0
```

\endcollapsible

\startcollapsible{Building (click to expand)}

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
passing `-DRAPTOR_NATIVE_BUILD=OFF` to CMake.

\endcollapsible
