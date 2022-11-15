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
sketches are stored in a directory and will be used in computing an HIBF layout. The HIBF layout tries to minimize the
disk space consumption of the resulting index. The space is estimated using a k-mer count per user bin which represents
the potential denisity in a technical bin in an interleaved Bloom filter.

Using all default values a first call will look like:

```bash
raptor layout --input-file all_bin_path.txt --tmax 64
```

The `input-file` looks exactly as in our previous calls of `raptor index`; it contains all the paths of our database
files.

\todo Chopper braucht ein `input_data.tsv` input, wobei es momentan nur eine Spalte (mit den Pfaden) gibt, also geht auch `.txt`.

The parameter `--tmax` limits the number of technical bins on each level of the HIBF. Choosing a good `tmax` is not
trivial. The smaller `tmax`, the more levels the layout needs to represent the data. This results in a higher space
consumption of the index. While querying each individual level is cheap, querying many levels might also lead to an
increased runtime. A good `tmax` is usually the square root of the number of user bins rounded to the next multiple of
`64`. Note that your `tmax` will be rounded to the next multiple of 64 anyway.

\note
At the expense of a longer runtime, you can enable the statistic mode that determines the best `tmax` using the option
`--determine-best-tmax`.
When this flag is set, the program will compute multiple layouts for `tmax` in `[64 , 128, 256, ... , tmax]` as well as
`tmax = sqrt(number of user bins)`. The layout algorithm itself only optimizes the space consumption. When determining
the best layout, we additionally keep track of the average number of queries needed to traverse each layout. This query
cost is taken into account when determining the best `tmax` for your data.
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

Each Bloom filter has a bit vector length that, across all Bloom filters, gives the size of the Interleaved Bloom
filter, which we can specify in the IBF case. Since the HIBF calculates the size of the index itself, it is no longer
possible to specify a size here. But we can offer the option to name the desired false positive rate with
`--false-positive-rate`.

\note These parameters must be set identically for `raptor index`.

\todo This is not checked at the moment?



Thus, for example, a call looks like this:


Parameter Tweaking:
- `--sketch-bits`
      The number of bits the HyperLogLog sketch should use to distribute the values into bins. Default: 12. Value
      must be in range [5,32].
- `--disable-sketch-output`
      Although the sketches will improve the layout, you might want to disable writing the sketch files to disk.
      Doing so will save disk space. However, you cannot use either --estimate-unions or --rearrange-user-bins in
      chopper layout without the sketches. Note that this option does not decrease run time as sketches have to be
      computed either way.

HyperLogLog Sketches:
To improve the layout, you can estimate the sequence similarities using HyperLogLog sketches.
- `--estimate-union`
      Use sketches to estimate the sequence similarity among a set of user bins. This will improve the layout
      computation as merging user bins that do not increase technical bin sizes will be preferred. Attention: Only
      possible if the directory [INPUT-PREFIX]_sketches is present.
- `--rearrange-user-bins`
      As a preprocessing step, rearranging the order of the given user bins based on their sequence similarity may
      lead to favourable small unions and thus a smaller index. Attention: Also enables --estimate-union and is
      only possible if the directory [INPUT-PREFIX]_sketches is present.

<- die machen das layout besser dauern aber laenger

  Parameter Tweaking:
- `--alpha`
      The layout algorithm optimizes the space consumption of the resulting HIBF but currently has no means of
      optimizing the runtime for querying such an HIBF. In general, the ratio of merged bins and split bins
      influences the query time because a merged bin always triggers another search on a lower level. To influence
      this ratio, alpha can be used. The higher alpha, the less merged bins are chosen in the layout. This
      improves query times but leads to a bigger index. Default: 1.2.
- `--max-rearrangement-ratio`
      When the option --rearrange-user-bins is set, this option can influence the rearrangement algorithm. The
      algorithm only rearranges the order of user bins in fixed intervals. The higher --max-rearrangement-ratio,
      the larger the intervals. This potentially improves the layout, but increases the runtime of the layout
      algorithm. Default: 0.5. Value must be in range [0,1].

A call could then look like this:
```bash
raptor layout --input-file all_bin_path.txt \
              --tmax 64 \
              --kmer-size 16 \
              --sketch-bits 5 \
              --num-hash-functions 3 \
              --false-positive-rate 0.25 \
              --estimate-union \
              --rearrange-user-bins \
              --alpha 1.5 \
              --max-rearrangement-ratio 0.25 \
              --threads 4 \
              --output-filename binning.layout
```

An assignment follows in the index tutorial on the HIBF.


### Parallelization

Raptor supports parallelization. By specifying `--threads`, for example, the k-mer hashes are processed simultaneously.
