# Create a layout with Raptor {#tutorial_layout}

<!--
SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

You will learn how to construct a Raptor layout of large collections of nucleotide sequences.

\attention
This is a tutorial only if you want to use the advanced option `--hibf` for the `raptor index`.
You can skip this chapter if you want to use raptor with the default IBF.

\tutorial_head{Easy, 30 min, \ref tutorial_first_steps,
[<b>Interleaved Bloom Filter (IBF)</b>](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1interleaved__bloom__filter.html#details)\,
[<b>Hierarchical Interleaved Bloom Filter (HIBF)</b>](https://hibf.vercel.app/dev/html/classseqan_1_1hibf_1_1hierarchical__interleaved__bloom__filter.html#details)}


[TOC]

# IBF vs HIBF

Raptor works with the Interleaved Bloom Filter by default. A new feature is the Hierarchical Interleaved Bloom Filter
(HIBF). This uses an almost always more space-saving method of storing
the bins (except if the input samples are all of the same size). It distinguishes between the *user bins*, which reflect
the individual input samples, and the *technical bins*, which are physical storage units within the HIBF.
*Technical bins* may store a single user bin, a split part of a user bin, or several (merged) user bins. This is
especially useful when samples vary dramatically in size.

To use the HIBF, a layout must be created.

# Create a Layout of the HIBF

To realise this distinction between user bins and technical bins, a layout must be calculated before creating an index.
For this purpose we have developed our own tool [Chopper](https://www.seqan.de/apps/chopper.html) and integrated it into
raptor. So you can simply call it up with `raptor layout` without having to install Chopper separately.

## General Idea & main parameters

\image html hibf.svg width=40%

The figure above shows the storage of the user bins in the technical bins. The resulting tree represents the layout.
The first step is to estimate the number of (representative) k-mers per user bin by computing
[HyperLogLog (HLL) sketches](https://dx.doi.org/10.46298/dmtcs.3545) of the input data. These HLL
sketches are used for computing an HIBF layout. We will go into more detail later \ref HLL.
The HIBF layout tries to minimize the disk space consumption of the resulting index. The space is estimated
using a k-mer count per user bin which represents the potential density in a technical bin in an Interleaved Bloom
filter.

\note
The term representative indicates that the k-mer content could be transformed by a function which reduces its size and
distribution, e.g. by using minimizers. Raptor currently offers canonical k-mers (e.g., k=20, w=20) and
minimisers (e.g., k=20, w=24).

Using all default values a first call will look like:

\snippet script.sh 02_layout_snippet_1

The `input` looks exactly as in our previous calls of `raptor index`; it contains all the paths of our database
files.

`raptor layout` chooses sensible defaults that showed to work well for the average use case.

We then get the resulting layout (default: `layout.txt`) as an output file, which we then pass to Raptor to create the
index. You can change this default with `--output`.

\note
Raptor also has a help page, which can be accessed as usual by typing `raptor layout -h` or `raptor layout --help`.

\assignment{Assignment 1: Create a first layout}
Lets take the bigger example form the introduction \ref tutorial_first_steps and create a layout for it.
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
And use the data of the `1024` Folder.

\hint
First we need a file with all paths to the fasta files. For this use the command:

\snippet script.sh 02_layout_snippet_2

\endhint

Then first determine the best tmax value and calculate a layout with default values and this tmax.
\endassignment

\solution
Your `all_bin_paths.txt` should look like:
```txt
1024/bins/bin_0001.fasta
1024/bins/bin_0002.fasta
1024/bins/bin_0003.fasta
...
1024/bins/bin_1021.fasta
1024/bins/bin_1022.fasta
1024/bins/bin_1023.fasta
```

/note
Sometimes it would be better to use the absolute paths instead.

And you should have run:
\snippet script.sh 02_layout_snippet_3
With the resulting `layout.txt`:
```bash
## ### Parameters ###
## number of user bins = 1023
## number of hash functions = 2
## false positive rate = 0.05
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	2.00	1.00	2.00	12.8MiB
```

Your directory should look like this:
```bash
$ ls
1024/                    all_bin_paths.txt        chopper_sketch.count     mini/
64/                      layout.txt              chopper_sketch_sketches/
```

\note
We will use this mini-example in the following.
Therefore, we recommend not deleting the files including the built indexes.

\endsolution

## Additional parameters

To create an index and thus a layout, the individual samples of the data set are chopped up into k-mers and determine in
their so-called bin the specific bit setting of the Bloom Filter by passing them through hash functions. This means that
a k-mer from sample `i` marks in bin `i` with `j` hash functions `j` bits with a `1`.

Example:
Query ACGT with kmers ACG, CGT.
Sample 1: AATGT, Sample 2: ACCGT, Sample 3: ACGTA
2 Hash funktions
-> Bins 1 to 3 for ACG could look like |0000|0000|0101|
   Bins 1 to 3 for CGT could look like |0000|0110|1100|
-> The query seems to match Sample 3

If a query is then searched, its k-mers are thrown into the hash functions and looked at in which bins it only points
to ones. This can also result in false positives. Thus, the result only indicates that the query is probably part of a
sample.

This principle also applies to the Hierarchical Interleaved Bloom Filter, except that the bins are then stored even more
efficiently as described above and this is described by the layout. This means that you already have to know some
parameters for the layout, which you would otherwise specify in the index:

With `--kmer` you can specify the length of the k-mers, which should be long enough to avoid random hits.
By using multiple hash functions, you can sometimes further reduce the possibility of false positives
(`--hash`). We found a useful [Bloom Filter Calculator](https://hur.st/bloomfilter/) to get a calculation
if it could help. As it is not ours, we do not guarantee its accuracy.
To use this calculator the number of inserted elements is the number of kmers in a single bin and you should use the
biggest bin to be sure.

Each Bloom Filter has a bit vector length that, across all Bloom Filters, gives the size of the Interleaved Bloom
Filter, which we can specify in the IBF case. Since the HIBF calculates the size of the index itself, it is no longer
possible to specify a size here. But we can offer the option to name the desired false positive rate with
`--fpr`.

\note These parameters must be set identically for `raptor index`.

A call could then look like this:
\snippet script.sh 02_layout_snippet_4

### Parallelization

Raptor supports parallelization. By specifying `--threads`, for example, the k-mer hashes are processed simultaneously.

\assignment{Assignment 2: Create a more specific layout}
Now lets run the above example with more parameters.

Use the same `all_bin_paths.txt` and create a `binning2.out`. Take a kmer size of `16`, `3` hash functions, a false
positive rate of `0.1` and use `2` threads.
\endassignment

\solution
And you should have run:
\snippet script.sh 02_layout_snippet_5
Your directory should look like this:
```bash
$ ls
1024/                    all_bin_paths.txt        binning2.layout          chopper_sketch_sketches/
64/                      layout.txt              chopper_sketch.count     mini/
```
\endsolution

### HyperLogLog sketch {#HLL}

The first step is to estimate the number of (representative) k-mers per user bin by computing
[HyperLogLog (HLL) sketches](http://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf) of the input data. These HLL
sketches are stored in a directory and will be used in computing an HIBF layout.

We will also give a short explanation of the HLL sketches here to explain the possible parameters, whereby each bin is
sketched individually.

\note
Most parameters are advanced and only need to be changed if the calculation takes significantly too long or the memory
usage is too high.

So the question is how many elements of our multiset are identical?
With exact calculation, we need more storage space and runtime than with the HLL estimate.
So, to find this out, we form (binary) 64 bit hash values of the data. These are equally distributed over all possible
hash values. If you go through this hashed data, you can then estimate how many different elements you have seen so far
only by reading leading zeros. (For the i'th element with `p` leading zeros, it is estimated that you have seen
\f$2^p\f$ different elements). You then simply save the maximum of these (\f$2^{p_{max}}\f$ different elements).

However, if we are unlucky and come across a hash value that consists of only `0`'s, then \f$p_{max}\f$ is of course
maximum of all possible hash values, no matter how many different elements are actually present.
To avoid this, we cut the stream of hash values into `m` substreams and use the first `b` bits of each hash value to
determine into which substream it belongs. From these, we calculate the *harmonic mean* as the total \f$p_{max}\f$.

We can influence this m with `--sketch-bits`. `m` must be a power of two so that we can divide the `64` bit evenly, so
we use `--sketch-bits` to set a `b` with \f$m = 2^b\f$.

If we choose our `b` (`m`) to be very large, then we need more memory but get higher accuracy. (Storage consumption is
growing exponentially.) In addition, calculating the layout can take longer with a high `b` (`m`).  If we have many user
bins and observe a long layout computation time, then it is worth choosing a somewhat smaller `b` (`m`). Furthermore,
the relative error of the HLL estimate increases with a decreasing `b` (`m`). Based on our benchmarks, we believe that
anything above `m = 512` should be fine.

#### Advanced options for HLL sketches

The following options should only be touched if the calculation takes a long time.

We have implemented another preprocessing that summarises the technical bins even better with regard to the similarities
of the input data. This can be switched off with the flag `--skip-similarity-preprocessing` if it costs too much
runtime.

With `--max-rearrangement-ratio` you can further influence a part of the preprocessing (value between `0` and `1`). If
you set this value to `1`, it is switched off. If you set it to a very small value, you will also need more runtime and
memory. If it is close to `1`, however, just little re-arranging is done, which could result in a less memory-efficient
layout. This parameter should only be changed if the layouting takes to much memory or time, because there it can have a
large influence.

One last observation about these advanced options: If you expect hardly any similarity in the data set, then the
similarity preprocessing makes very little difference.

### Advanced option: tmax

The parameter `--tmax` limits the number of technical bins on each level of the HIBF. Choosing a good \f$t_{max}\f$ is
not trivial. The smaller \f$t_{max}\f$, the more levels the layout needs to represent the data. This results in a higher
space consumption of the index. While querying each individual level is cheap, querying many levels might also lead to
an increased runtime. A good (and current default) \f$t_{max}\f$ is usually the square root of the number of user bins
rounded to the next multiple of `64`. Note that your \f$t_{max}\f$ will be rounded to the next multiple of 64.

At the expense of a longer runtime, you can enable the statistic mode that determines the best \f$t_{max}\f$ using the
option `--determine-best-tmax`.
When this flag is set, the program will compute multiple layouts for \f$t_{max}\f$ in `[64 , 128, 256, ... , tmax]` as
well as `tmax = sqrt(number of user bins)`. The layout algorithm itself only optimizes the space consumption. When
determining the best layout, we additionally keep track of the average number of queries needed to traverse each layout.
This query cost is taken into account when determining the best \f$t_{max}\f$ for your data.
Note that the option `--tmax` serves as upper bound. Once the layout quality starts dropping, the computation is
stopped. To run all layout computations, pass the flag `--force-all-binnings`.
The flag `--force-all-binnings` forces all layouts up to `--tmax` to be computed, regardless of the layout quality. If
the flag `--determine-best-tmax` is not set, this flag is ignored and has no effect.

### Advanced Option: alpha

You should only touch the parameter `--alpha` if you have understood very well how the layout works and you are
dissatisfied with the resulting index, e.g. there is still a lot of space in RAM but the index is very slow.

The layout algorithm optimizes the space consumption of the resulting HIBF but currently has no means of optimizing the
runtime for querying such an HIBF. In general, the ratio of merged bins and split bins influences the query time because
a merged bin always triggers another search on a lower level. To influence this ratio, alpha can be used.

Alpha is a parameter for weighting the storage space calculation of the lower-level IBFs. It functions as a lower-level
penalty, i.e. if alpha is large, the DP algorithm tries to avoid lower levels, which at the same time leads to the
top-level IBF becoming somewhat larger. This improves query times but leads to a bigger index.
