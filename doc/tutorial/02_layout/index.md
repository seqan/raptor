# Create a layout with Raptor {#tutorial_layout}

You will learn how to construct a Raptor layout of large collections of nucleotide sequences.

\attention
This is a tutorial only if you want to use the advanced option `--hibf` for the `raptor index`.
You can skip this chapter if you want to use raptor with the default IBF.

\tutorial_head{Easy, 30 min, \ref tutorial_first_steps,
[<b>Interleaved Bloom Filter (IBF)</b>](https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1interleaved__bloom__filter.html#details)\,
\ref raptor::hierarchical_interleaved_bloom_filter "Hierarchical Interleaved Bloom Filter (HIBF)"}

[TOC]

# IBF vs HIBF

Raptor works with the Interleaved Bloom Filter by default. A new feature is the Hierarchical Interleaved Bloom Filter
(HIBF) (raptor::hierarchical_interleaved_bloom_filter). This uses a more space-saving method of storing the bins. It
distinguishes between the user bins, which reflect the individual samples as before, and the so-called technical bins,
which throw some bins together. This is especially useful when there are samples of very different sizes.

To use the HIBF, a layout must be created

# Create a Layout of the HIBF

To realise this distinction between user bins and technical bins, a layout must be calculated before creating an index.
For this purpose we have developed our own tool [Chopper](https://www.seqan.de/apps/chopper.html) and integrated it into
raptor. So you can simply call it up with `raptor layout` without having to install Chopper separately.

## General Idea & main parameters

\image html hibf.svg width=40%

The figure above shows the storage of the user bins in the technical bins. The resulting tree represents the layout.
The first step is to estimate the number of (representative) k-mers per user bin by computing
[HyperLogLog (HLL) sketches](http://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf) of the input data. These HLL
sketches are stored in a directory and will be used in computing an HIBF layout. We will go into more detail later
\ref HLL. The HIBF layout tries to minimize the disk space consumption of the resulting index. The space is estimated
using a k-mer count per user bin which represents the potential denisity in a technical bin in an Interleaved Bloom
filter.

Using all default values a first call will look like:

```bash
raptor layout --input-file all_bin_path.txt --tmax 64
```

The `input-file` looks exactly as in our previous calls of `raptor index`; it contains all the paths of our database
files.

\todo Chopper braucht ein `input_data.tsv` input, wobei es momentan nur eine Spalte (mit den Pfaden) gibt, also geht auch `.txt`.

The parameter `--tmax` limits the number of technical bins on each level of the HIBF. Choosing a good \f$t_{max}\f$ is
not trivial. The smaller \f$t_{max}\f$, the more levels the layout needs to represent the data. This results in a higher
space consumption of the index. While querying each individual level is cheap, querying many levels might also lead to
an increased runtime. A good \f$t_{max}\f$ is usually the square root of the number of user bins rounded to the next
multiple of `64`. Note that your \f$t_{max}\f$ will be rounded to the next multiple of 64 anyway.

\note
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

We then get the resulting layout (default: `binning.out`) as an output file, which we then pass to Raptor to create the
index. You can change this default with `--output-filename`.

\note
Raptor also has a help page, which can be accessed as usual by typing `raptor layout -h` or `raptor layout --help`.

## Additional parameters

To create an index and thus a layout, the individual samples of the data set are chopped up into k-mers and determine in
their so-called bin the specific bit setting of the Bloom Filter by passing them through hash functions. This means that
a k-mer from sample `i` marks in bin `i` with `j` hash functions `j` bits with a `1`.
If a query is then searched, its k-mers are thrown into the hash functions and looked at in which bins it only points
to ones. This can also result in false positives. Thus, the result only indicates that the query is probably part of a
sample.

This principle also applies to the Hierarchical Interleaved Bloom Filter, except that the bins are then stored even more
efficiently as described above and this is described by the layout. This means that you already have to know some
parameters for the layout, which you would otherwise specify in the index:

With `--kmer-size` you can specify the length of the k-mers, which should be long enough to avoid random hits.
By using multiple hash functions, you can sometimes further reduce the possibility of false positives
(`--num-hash-functions`). We found a useful [Bloom Filter Calculator](https://hur.st/bloomfilter/) to get a calculation
if it could help. As it is not ours, we do not guarantee its accuracy.

Each Bloom Filter has a bit vector length that, across all Bloom Filters, gives the size of the Interleaved Bloom
Filter, which we can specify in the IBF case. Since the HIBF calculates the size of the index itself, it is no longer
possible to specify a size here. But we can offer the option to name the desired false positive rate with
`--false-positive-rate`.

\note These parameters must be set identically for `raptor index`.

\todo This is not checked at the moment?

### Parallelization

Raptor supports parallelization. By specifying `--threads`, for example, the k-mer hashes are processed simultaneously.

### HyperLogLog sketch {#HLL}

The first step is to estimate the number of (representative) k-mers per user bin by computing
[HyperLogLog (HLL) sketches](http://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf) of the input data. These HLL
sketches are stored in a directory and will be used in computing an HIBF layout.

We will also give a short explanation of the HLL sketches here to explain the possible parameters.

\note
Most parameters are advanced and only need to be changed if the calculation takes significantly too long or the memory
usage is too high.

So the question is how many elements of our multiset are identical?
With exact calculation, we need more storage space and runtime than with the HLL estimate.
So, to find this out, we form (binary) 64 bit hash values of the data. These are equally distributed over all possible
hash values. If you go through this hashed data, you can then estimate how many different elements you have seen so far
only by reading leading zeros. (For the i'th element with `p` leading zeros, it is estimated that you have seen
\f$2^p\f$ different elements). You then simply save the maximum of these (\f$2^{p_{max}}\f$ different elements).

<!-- bei p (at least) leading zeros ist die wahrscheinlichkeit dieses vorkommens genau 1/(2^p) -->

However, if we are unlucky and come across a hash value that consists of only `0`'s, then \f$p_{max}\f$ is of course
maximum of all possible hash values, no matter how many different elements are actually present.
To avoid this, we cut each hash value into `m` parts and calculate the \f$p_{max}\f$ over each of these parts. From
these we then calculate the *harmonic mean* as the total \f$p_{max}\f$.

We can influence this m with `--sketch-bits`. `m` must be a power of two so that we can divide the `64` bit evenly, so
we use `--sketch-bits` to set a `b` with \f$m = 2^b\f$.

If we choose our `b` (`m`) to be very large, then we need more memory but get higher accuracy. (Storage consumption is
growing exponentially.) In addition, calculating the layout can take longer with a high `b` (`m`). If we have many user
bins and observe a long runtime, then it is worth choosing a somewhat smaller `b` (`m`).

\todo
Wird bisher ein sketch über alles berechnet oder einzelne sketches die gemerged werden? Laut Felix ist das mergen ebenfalls sehr schnell. (zb 10 genome sketchen und dann mergen)

A call could then look like this:
```bash
raptor layout --input-file all_bin_path.txt \
              --tmax 64 \
              --kmer-size 16 \
              --sketch-bits 5 \
              --num-hash-functions 3 \
              --false-positive-rate 0.25 \
              --threads 4 \
              --output-filename binning.layout
```

An assignment follows in the index tutorial on the HIBF.

#### Advanced options for HLL sketches

The following options should only be touched if the calculation takes a long time.

We have implemented another preprocessing that summarises the technical bins even better with regard to the similarities
of the input data. This can be switched off with the flag `--skip-similarity-preprocessing` if it costs too much
runtime.

\todo
Add parameter `--skip-similarity-preprocessing` instead of `--estimate-union` and `--rearrange-user-bins`. Set it as
advanced.

\todo
`--disable-sketch-output` wahrscheinlich unsinnig, da nur der zwischenstand zwischen count und layout

With `--max-rearrangement-ratio` you can further influence a part of the preprocessing (value between `0` and `1`). If
you set this value to `1`, it is switched off. If you set it to a very small value, you will also need more runtime and
memory. If it is close to `1`, however, just little re-arranging is done, which could be bad. In our benchmarks, however,
we were not able to determine a too great influence, so we recommend that this value only be used for fine tuning.

\todo
Wenn r=1, dann `--rearrange-user-bins` aus, daher die flag nicht nötig.

One last observation about these advanced options: If you expect hardly any similarity in the data set, then the
similarity preprocessing makes very little difference.

### Another advanced Option: alpha

You should only touch the parameter `--alpha` if you have understood very well how the layout works and you are
dissatisfied with the resulting index, e.g. there is still a lot of space in RAM but the index is very slow.

The layout algorithm optimizes the space consumption of the resulting HIBF but currently has no means of optimizing the
runtime for querying such an HIBF. In general, the ratio of merged bins and split bins influences the query time because
a merged bin always triggers another search on a lower level. To influence this ratio, alpha can be used.

Alpha is a parameter for weighting the storage space calculation of the lower-level IBFs. It functions as a lower-level
penalty, i.e. if alpha is large, the DP algorithm tries to avoid lower levels, which at the same time leads to the
top-level IBF becoming somewhat larger. This improves query times but leads to a bigger index.
