# Setup {#setup}

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

### Prerequisites
* CMake >= 3.18
* GCC 11, 12 or 13 (most recent minor version)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in depth
information.

### Download current main branch
```bash
git clone https://github.com/seqan/raptor
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
passing `-DRAPTOR_NATIVE_BUILD=OFF` to CMake.
