# Raptor [![build status][1]][2] [![codecov][3]][4] [![install with bioconda][5]][6] [![install with brew][7]][8]

[1]: https://img.shields.io/github/workflow/status/seqan/raptor/CI%20on%20Linux/main?style=flat&logo=github&label=Raptor%20CI
[2]: https://github.com/seqan/raptor/actions?query=branch%3Amain
[3]: https://codecov.io/gh/seqan/raptor/branch/main/graph/badge.svg?token=SJVMYRUKW2
[4]: https://codecov.io/gh/seqan/raptor
[5]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[6]: #install-with-bioconda-linux
[7]: https://img.shields.io/badge/install%20with-brew-brightgreen.svg?style=flat
[8]: #install-with-brew-linux-macos

### A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences

## Download and Installation
There may be performance benefits when compiling from source as the build can be optimized for the host system.

### Install with [bioconda](https://bioconda.github.io/recipes/raptor/README.html) (Linux)

```bash
conda install -c bioconda -c conda-forge raptor
```

### Install with [brew](https://brew.sh/) (Linux, macOS)

```bash
brew install brewsci/bio/raptor
```

### Compile from source
<details><summary>Prerequisites (click to expand)</summary>

* CMake >= 3.18
* GCC 10, 11 or 12 (most recent minor version)
* git

Refer to the [Seqan3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in depth
information.
</details>

<details><summary>Download current main branch (click to expand)</summary>

```bash
git clone https://github.com/seqan/raptor
git submodule update --init
```

</details>

<details><summary>Download specific version (click to expand)</summary>

E.g., for version `1.1.0`:
```bash
git clone --branch raptor-v1.1.0 --recurse-submodules https://github.com/seqan/raptor
```
Or from within an existing repository
```bash
git checkout raptor-v1.1.0
```
</details>

<details><summary>Building (click to expand)</summary>

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

</details>

## Example Data and Usage
A toy data set (124 MiB compressed, 983 MiB decompressed) can be found
[here](https://ftp.imp.fu-berlin.de/pub/seiler/raptor/).

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

The `bins` folder contains a FASTA file for each bin and the `reads` directory contains a FASTQ file for each bin
containing reads from the respective bin (with 2 errors).
Additionally, `mini.fastq` (5 reads of all bins), `all.fastq` (concatenation of all FASTQ files) and `all10.fastq`
(`all.fastq` repeated 10 times) are provided in the `reads` folder.

In the following, we will use the `64` data set.
To build an index over all bins, we first prepare a file that contains one file path per line
(a line corresponds to a bin) and use this file as input:
```
seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63 > all_bin_paths.txt
raptor build --kmer 19 --window 23 --size 8m --output raptor.index all_bin_paths.txt
```

You may be prompted to enable or disable automatic update notifications. For questions, please consult
[the SeqAn documentation](https://github.com/seqan/seqan3/wiki/Update-Notifications).

Afterwards, we can search for some reads:

```
raptor search --error 2 --index raptor.index --query example_data/64/reads/mini.fastq --output search.output
```

The output starts with a header section (lines starting with `#`). The header maps a number to each input file.
After the header section, each line of the output consists of the read ID (in the toy example these are numbers) and
the corresponding bins in which they were found:
```text
#0	example_data/64/bins/bin_00.fasta
#1	example_data/64/bins/bin_01.fasta
...
#62	example_data/64/bins/bin_62.fasta
#63	example_data/64/bins/bin_63.fasta
#QUERY_NAME	USER_BINS
0	0
1	0
2	0
3	0
4	0
16384	1
...
1015812	62
1032192	63
1032193	63
1032194	63
1032195	63
1032196	63
```

For a list of options, see the help pages:
```console
raptor --help
raptor build --help
raptor search --help
raptor upgrade --help
```

### Preprocessing the input
We offer the option to precompute the minimisers of the input files. This is useful to build indices of big datasets
(in the range of several TiB) and also allows an estimation of the needed index size since the amount of minimisers is
known.
Following above example, we would change the build step as follows:

First we precompute the minimisers and store them in a directory:
```
seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63 > all_bin_paths.txt
raptor build --kmer 19 --window 23 --compute-minimiser --output precomputed_minimisers all_bin_paths.txt
```

Then we run the build step again and use the computed minimisers as input:
```
seq -f "precomputed_minimisers/bin_%02g.minimiser" 0 1 63 > all_minimiser_paths.txt
raptor build --size 8m --output minimiser_raptor.index all_minimiser_paths.txt
```

The preprocessing applies the same cutoffs as used in Mantis
([Pandey et al., 2018](https://doi.org/10.1016/j.cels.2018.05.021)).
This means that only minimisers that occur more often than the cutoff specifies are included in the output.
If you wish to process all minimisers, you can use `--disable-cutoffs`.

### Partitioned indices
To reduce the overall memory consumption, the index can be divided into multiple (a power of two) parts.
This can be done by passing `--parts n` to `raptor build`, where `n` is the number of parts you want to create.
This will create `n` files, each representing one part of the index. The `--size` parameter describes the overall size
of the index. For example, `--size 8g --parts 4` will create four 2 GiB indices. This will reduce the memory consumption
of `raptor build` and `raptor search` by approximately 6 GiB, since there will only be one part in memory at any given
time. `raptor search` will automatically detect the parts, and does not need any special parameters.

### Upgrading the index (v1.1.0 to v2.0.0)
An old index can be upgraded by running `raptor upgrade` and providing some information about how the index was
constructed.

### SOCKS interface
We implement the core interface of [SOCKS](https://gitlab.ub.uni-bielefeld.de/gi/socks).
For a list of options, see the help pages:
```console
raptor socks build --help
raptor socks lookup-kmer --help
```

## Authorship and Copyright
Raptor is being developed by [Enrico Seiler](mailto:enrico.seiler@fu-berlin.de), but also incorporates much work from
other members of [SeqAn](https://www.seqan.de).

### Citation
In your academic works (also comparisons and pipelines) please cite:
  * *Raptor: A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences*;
    Enrico Seiler, Svenja Mehringer, Mitra Darvish, Etienne Turc, and Knut Reinert; iScience 2021 24 (7): 102782.
    doi: https://doi.org/10.1016/j.isci.2021.102782

### RECOMB 2021
Raptor was presented at the [25th International Conference on Research in Computational Molecular Biology][recomb_url]:
  * [Slides][recomb_slides]
  * [Pre-recorded version of the talk][recomb_talk]

Please see the [License](#license) section for information on allowed use.

[recomb_url]: https://www.recomb2021.org/
[recomb_slides]: https://box.fu-berlin.de/s/TtM3Raxixm35Syy
[recomb_talk]: https://box.fu-berlin.de/s/YJFQnwqdE5q2Tym

### Supplementary
The subdirectory `util` contains applications and scripts related to the paper.

### License
Raptor is open source software. However, certain conditions apply when you (re-)distribute and/or modify Raptor,
please see the [license](https://github.com/seqan/raptor/blob/main/LICENSE.md).

## Sponsorships

[![Vercel][vercel_badge]][vercel_website]

[vercel_badge]: https://raw.githubusercontent.com/seqan/raptor/main/test/documentation/.vercel/powered-by-vercel.svg "Powered by Vercel"
[vercel_website]: https://vercel.com/?utm_source=seqan&utm_campaign=oss

Vercel is kind enough to build and host our documentation and even provide preview-builds within our pull requests.
Check them out!
