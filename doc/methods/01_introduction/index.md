# First steps with Raptor {#tutorial_first_steps}

***Learning Objective:***

This tutorial walks you through small Raptor applications. It is intended to give you a short overview of what to expect
in the other tutorials and how to use this documentation.

\tutorial_head{Easy, 30 min, \ref setup, }

*Every page in the tutorials begins with this section. It is recommended that you do the "prerequisite tutorials"
before the current one. You should also have a look at the links provided in "recommended reading" and maybe keep
them open in separate tabs/windows for reference.*

[TOC]

# Example Data and Usage

Raptor is a fast and space-efficient pre-filter for querying very large collections of nucleotide sequences. Before we
go into the various applications and possibilities of Raptor in the next tutorials, we want to run Raptor in this
tutorial with a first minimal example and give a first overview of some possibilities.

A toy data set (124 MiB compressed, 983 MiB decompressed) can be found
[here](https://ftp.imp.fu-berlin.de/pub/seiler/raptor/).

\snippet script.sh 01_introduction_snippet_1

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

\snippet script.sh 01_introduction_snippet_2

You may be prompted to enable or disable automatic update notifications. For questions, please consult
[the SeqAn documentation](https://github.com/seqan/seqan3/wiki/Update-Notifications).

Afterwards, we can search for some reads:

\snippet script.sh 01_introduction_snippet_3

The output starts with metainformation on how the output file was produced (lines starting with `##`).
Following this, there is a header section (lines starting with `#`) that maps a number to each input file.
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

\note
For a list of options, see the help pages:
```console
raptor --help
raptor layout --help
raptor build --help
raptor prepare --help
raptor search --help
raptor upgrade --help
```

## Preprocessing the input
We offer the option to precompute the
[minimisers](https://docs.seqan.de/seqan/3.0.3/group__views.html#ga191fcd1360fc430441567f3ed0f371d1) of the input files.
This is useful to build indices of big datasets (in the range of several TiB).
Following above example, we would change the build step as follows:

First we precompute the minimisers and store them in a directory:

\snippet script.sh 01_introduction_snippet_4

Then we run the build step again and use the computed minimisers as input:

\snippet script.sh 01_introduction_snippet_5

The preprocessing may apply the same cutoffs as used in Mantis
([Pandey et al., 2018](https://doi.org/10.1016/j.cels.2018.05.021)).
This means that only minimisers that occur more often than the cutoff specifies are included in the output.
If you wish to do so, you can use `--use-filesize-dependent-cutoff`.

By default, `--kmer-count-cutoff 1` is used, which will include all minimisers in the output.

\note
If you want to learn more about Minimiser, we recommend the corresponding SeqAn3 tutorial:
[Minimiser](https://docs.seqan.de/seqan/3.0.3/tutorial_minimiser.html).

## Partitioned indices
To reduce the overall memory consumption, the index can be divided into multiple (a power of two) parts.
This can be done by passing `--parts n` to `raptor build`, where `n` is the number of parts you want to create.
This will create `n` files, each representing one part of the index. This will reduce the memory consumption
of `raptor build` and `raptor search` by roughly `n`, since there will only be one part in memory at any given
time. `raptor search` will automatically detect the parts, and does not need any special parameters.

## Upgrading the index
An old index can be upgraded by running `raptor upgrade` and providing some information about how the index was
constructed.
