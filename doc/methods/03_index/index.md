# Indexing with Raptor {#tutorial_index}

<!--
SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

You will learn how to construct a Raptor index of large collections of nucleotide sequences.

\tutorial_head{Easy, 30 min, \ref tutorial_first_steps \, for using HIBF: \ref tutorial_layout,
[<b>Interleaved Bloom Filter (IBF)</b>](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1interleaved__bloom__filter.html#details)\,
[<b>Hierarchical Interleaved Bloom Filter (HIBF)</b>](https://hibf.vercel.app/dev/html/classseqan_1_1hibf_1_1hierarchical__interleaved__bloom__filter.html#details)}

[TOC]

# Indexing large collections of nucleotide sequences

## A first index

Raptor consists of two methods, `build` and `search`. The former creates an index over a given dataset so that it can be
searched efficiently with `raptor search`.
In this tutorial, we will look at `raptor build` in more detail.

Raptor can be used as a pre-filter for applications where searching the complete dataset is not feasible. For example,
in read mapping you might want to compare your genome to the genome of 100'000 other people.
It can also be used for metagenomic classification, i.e., which microbes are in a single sample.

\see Related tools:
- Metagenomics: [Ganon/Ganon2](https://github.com/pirovc/ganon), a tool for comparative metagenomics analysis in short
time using the whole of the quickly growing number of assembled sequences openly available.
- RNA-seq expression analysis: [Needle](https://www.seqan.de/apps/needle.html), a tool for storing sequencing
experiments in such a way that approximate quantification of large sequencing data sets is possible.

Regardless of the application, we start with a data set of different nucleotide sequences over which we want to build an
index.

These sequences are typically available as
[FASTA](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1format__fasta.html#a86790afd92e0229cccbc20be20d5758a) or
[FASTQ](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1format__fastq.html#details) files, but can be in any
format supported by
[`seqan3::sequence_file_input`](https://docs.seqan.de/seqan/3-master-user/tutorial_sequence_file.html).

We summarise these in a list of paths and can then create a first index using the default values of raptor:
\snippet script.sh 03_index_snippet_1

\assignment{Assignment 1: Create example files and index them}
Use this script to create FASTA files.
\snippet script.sh 03_index_snippet_2
Then create an all_paths.txt and run raptor build.

\note
Check the help page to read what the expected content of `all_paths.txt` looks like, which can be accessed as usual by
typing `raptor build -h` or `raptor build --help`.

\endassignment

\solution
Your `all_paths.txt` might look like:
```txt
mini_1.fasta
mini_2.fasta
mini_3.fasta
```

/note
Sometimes it would be better to use the absolute paths instead.

And you should have run:
\snippet script.sh 03_index_snippet_3
Your directory should look like this:
```bash
tmp$ ls
all_paths.txt   mini_1.fasta    mini_2.fasta    mini_3.fasta    raptor.index
```

\note
We will use this mini-example in the following, both with further parameters and then for `raptor search`. Therefore, we
recommend not deleting the files including the built indexes.

\endsolution

## General idea & main parameters

\image html ibf.svg width=90%

Before explaining parameters, we would like to briefly explain the general idea of the Raptor index.

If we want to check whether a query is contained in a sample, we can use a
[Bloom Filter (BF)](https://en.wikipedia.org/wiki/Bloom_filter) to create an index for the sample. Although this only
gives us a probability, it saves us a time-consuming complete mapping.
If our data set consists of many samples, there is a BF for each sample. By default, Raptor uses an
[Interleaved Bloom Filter (IBF)](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1interleaved__bloom__filter.html#details),
which is an efficient way to store these many Bloom Filters, called raptor index. Another possibility is to use the
Hierarchical Interleaved Bloom Filter (HIBF), more about this later in \ref hibf.

To create the index, the individual samples of the data set are chopped up into k-mers and determine in their so-called
bin the specific bit setting of the Bloom Filter by passing them through the hash functions. This means that a k-mer
from sample `i` marks in bin `i` with `j` hash functions `j` bits with a `1`.
If a query is then searched, its k-mers are thrown into the hash functions and looked at in which bins it only points
to ones. This can also result in false positives. Thus, the result only indicates that the query is probably part of a
sample.

With `--kmer` you can specify the length of the k-mers, which should be long enough to avoid random hits.
By using multiple hash functions, you can sometimes further reduce the possibility of false positives (`--hash`). We
found a useful [Bloom Filter Calculator](https://hur.st/bloomfilter/) to get a calculation if it could help. As it is
not ours, we do not guarantee its accuracy.
To use this calculator the number of inserted elements is the number of kmers in a single bin and you should use the
biggest bin to be sure.

Each Bloom Filter has a bit vector length, which over all Bloom Filters gives the size of the Interleaved Bloom Filter,
which is automatically inferred. The lower the false positive rate, the bigger the index.

\snippet script.sh 03_index_snippet_4

\assignment{Assignment 2: Using parameters}
As our example is really tiny, lets run this with a k-mer size 4, three hash functions and a false positive rate of
0.01. Name the new index `raptor2.index`.
\endassignment

\solution
You should have run
\snippet script.sh 03_index_snippet_5

Your directory should look like this:
```bash
tmp$ ls -la
...  39B 10 Nov 16:42 all_paths.txt
... 107B 14 Nov 16:11 mini_1.fasta
... 107B 14 Nov 16:11 mini_2.fasta
... 107B 14 Nov 16:11 mini_3.fasta
... 1,2K 14 Nov 16:15 raptor.index
... 1,0M 14 Nov 16:18 raptor2.index
```
\endsolution

## Using minimisers

The k-mers can also be saved as [minimisers](https://docs.seqan.de/seqan/3-master-user/group__views.html#ga191fcd1360fc430441567f3ed0f371d1),
which saves space. To do this, first use `raptor prepare`. A minimiser works with windows, which means that you also
have to define their size `--window`. A window is always larger than a k-mer because it combines several k-mers.

\snippet script.sh 03_index_snippet_6

\assignment{Assignment 3: Minimisers}
Lets use a minimiser for our small example. Use for this a window size of 6 and a k-mer size of 4 and name the new
index minimiser.index.
\hint
For the all_minimiser_paths.txt you just need to add the `*.minimiser` files.
\endhint
\endassignment

\solution
You should have run:
\snippet script.sh 03_index_snippet_7

Your directory should look like this:
```bash
tmp$ ls -la
...  39B 10 Nov 16:42 all_paths.txt
... 107B 14 Nov 16:11 mini_1.fasta
... 107B 14 Nov 16:11 mini_2.fasta
... 107B 14 Nov 16:11 mini_3.fasta
... 8,0M 14 Nov 16:21 minimiser.index
... 256B 14 Nov 16:20 precomputed_minimisers/
...  19B 10 Nov 16:40 query.fasta
... 1,2K 14 Nov 16:15 raptor.index
... 1,0M 14 Nov 16:18 raptor2.index
```
\endsolution

\note
If you want to learn more about minimisers, take a look at the SeqAn3 tutorial for
[minimisers](https://docs.seqan.de/seqan/3-master-user/tutorial_minimiser.html).

### Advanced

When hashing a sequence, there may be positions that do not count towards the final hash value. A
[shape](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1shape.html) offers an easy way to define such patterns:
`--shape`.

\note
Raptor also has an advanced help page for the experimental content, which can be accessed as usual by typing
`raptor build -hh` or `raptor build --advanced-help`.

## Others

### Parallelization

Raptor supports parallelization. By specifying `--threads`, for example, the fastq-records are processed simultaneously.

### Partitioned indices

To reduce the overall memory consumption of the search, the index can be divided into multiple (a power of two) parts.
This can be done by passing `--parts n`, where n is the number of parts you want to create. This will create `n` files,
each representing one part of the index.
This will reduce the memory consumption of `raptor build` and `raptor search` by roughly `n`, since there will
only be one part in memory at any given time.
`raptor search` will automatically detect the parts, and does not need any special parameters.

## IBF vs HIBF {#hibf}

Raptor works with the Interleaved Bloom Filter by default. A new feature is the Hierarchical Interleaved Bloom Filter
(HIBF) (raptor::hierarchical_interleaved_bloom_filter), which is enabled when a layout file from `raptor layout` is used
as input. This uses a more
space-saving method of storing the bins. It distinguishes between the user bins, which reflect the individual samples as
before, and the so-called technical bins, which throw some bins together. This is especially useful when there are
samples of very different sizes.

To use the HIBF, a layout must be created before creating an index. We have written an extra tutorial for this
\ref tutorial_layout.

### HIBF indexing with the use of the layout

The layout replaces the `--input all_bin_path.txt` and is given instead with the layout: `--input layout.txt`.

We can set the desired false positive rate with `--fpr`. Thus, for example, a call looks like this:

\snippet script.sh 03_index_snippet_8

\assignment{Assignment 4: A default HIBF}
Since we cannot see the advantages of the hibf with our small example. And certainly not the differences when we change
the parameters. Let's not go back to our small example from above, but to the one from the introduction:

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
And use the data of the `1024` Folder and the two layouts we've created in the \ref tutorial_layout tutorial :
`layout.txt` and `binning2.layout`.

Lets use the HIBF with the default parameters and call the new indexes `hibf.index` and `hibf2.index`.
\note
Your kmer size, number of hash functions and the false positive rate must be the same as in the layout.

\hint
For the second layout we took a kmer size of `16`, `3` hash functions and a false positive rate of `0.1`.
\endhint
\endassignment

\solution
You should have run:
\snippet script.sh 03_index_snippet_9

Your directory should look like this:
```bash
$ ls -la
... 384B 12 Dez 13:44 ./
... 544B 12 Dez 12:14 ../
... 128B 23 Sep  2020 1024/
... 128B 23 Sep  2020 64/
...  25K 12 Dez 12:52 all_bin_paths.txt
...  37K 12 Dez 12:56 layout.txt
...  37K 12 Dez 13:26 binning2.layout
...  55K 12 Dez 13:26 chopper_sketch.count
...  32K 12 Dez 12:52 chopper_sketch_sketches/
...  13M 12 Dez 13:44 hibf.index
... 7,2M 12 Dez 13:44 hibf2.index
...  64B  8 Dez 16:24 mini/
```
\endsolution

\note
For a more detailed explanation of the Hierarchical Interleaved Bloom Filter (HIBF), please refer to the
`raptor::hierarchical_interleaved_bloom_filter` API.
