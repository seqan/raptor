# Raptor [![build status](https://github.com/seqan/raptor/workflows/Raptor%20CI/badge.svg?branch=master)](https://github.com/seqan/raptor/actions)
### A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences

## Download and Installation
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/raptor/README.html)

```bash
conda install -c bioconda -c conda-forge raptor
```

* There may be performance benefits when compiling from source, especially when using `-march=native` as compiler directive.
* Raptor is compatible with macOS but, due to compiler restrictions, Raptor is not available for macOS on conda.

### Quick Start without conda
<details><summary>Prerequisites</summary>

* CMake >= 3.8
* GCC 7, 8, 9 or 10 (most recent minor version)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in depth information.
</details>

<details><summary>Download current master branch</summary>

```bash
git clone --recurse-submodules https://github.com/seqan/raptor
```

</details>

<details><summary>Download specific version</summary>

E.g., for version `1.0.0`:
```bash
git clone --branch raptor-v1.0.0 --recurse-submodules https://github.com/seqan/raptor
```
Or from within an existing repository
```bash
git checkout raptor-v1.0.0
```
</details>

<details><summary>Building</summary>

```bash
cd raptor
mkdir -p build
cd build
cmake ..
make
```

The binary can be found in `bin`.

You may want to add the raptor executable yo your PATH:
```
export PATH=$(pwd)/bin:$PATH
raptor --version
```

</details>

### Example Data and Usage
A toy data set can be found [here](https://ftp.imp.fu-berlin.de/pub/seiler/raptor/).

```bash
wget https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
tar xfz example_data.tar.gz
```

After extraction, the `example_data` will look like:

```console
$ tree -L 2 example_data
example_data
├── 1024
│   ├── bins
│   └── reads
└── 64
    ├── bins
    └── reads
```

The `bins` folder contains a FASTA file for each bin and the `reads` directory contains a FASTQ file for each bin containing reads from the respective bin (with 2 errors).
Additionally, `mini.fastq` (5 reads of all bins), `all.fastq` (concatenation of all FASTQ files) and `all10.fastq` (`all.fastq` repeated 10 times) are provided in the `reads` folder.

In the following, we will use the `64` data set.
We can now build an index over all the bins:

```
raptor build --kmer 19 --window 23 --size 8m --output index.raptor $(seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63)
# You can replace `$(seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63)` by `example_data/64/bins/bin_{00..63}.fasta` if your shell supports this syntax.
# The equivalent command for 1,024 bins is `$(seq -f "example_data/1024/bins/bin_%04g.fasta" 0 1 1023)`
```
You may be prompted to enable or disable automatic update notifications. For questions, please consult [the SeqAn documentation](https://github.com/seqan/seqan3/wiki/Update-Notifications).

Afterwards, we can search for all reads from bin 1:

```
raptor search --kmer 19 --window 23 --error 2 --index index.raptor --query example_data/64/reads/mini.fastq --output search.output
```

Each line of the output consists of the read ID (in the toy example these are numbers) and the corresponding bins in which they were found:
```text
0       0,
1       0,
2       0,
3       0,
4       0,
16384   1,
...
1015812 62,
1032192 63,
1032193 63,
1032194 63,
1032195 63,
1032196 63,
```

For a list of options, see the help pages:
```console
raptor --help
raptor build --help
raptor search --help
```

## Authorship and Copyright
Raptor is being developed by [Enrico Seiler](mailto:enrico.seiler@fu-berlin.de), but also incorporates much work from
other members of [SeqAn](https://www.seqan.de).

### Citation
In your academic works (also comparisons and pipelines) please cite:
  * Seiler, E. et al. (2020). Raptor: A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences. bioRxiv 2020.10.08.330985. doi: https://doi.org/10.1101/2020.10.08.330985

### License
Raptor is open source software. However, certain conditions apply when you (re-)distribute and/or modify Raptor, please see the [license](https://github.com/seqan/raptor/blob/master/LICENSE.md).
