# raptor layout {#usage_layout}

[TOC]

\attention
A layout is only needed if you want to build an HIBF (\ref hibf_vs_ibf).

# Main Parameters

## -​-input-file {#usage_layout_input}

### List of sequence files

\include{doc} fragments/input_files_sequence.md

\note
Currently, `raptor build` accepts multiple paths per line in its input file. However, `raptor layout` does not.
A possible workaround is to concatenate your sequence files. To avoid changing your files, you can also use
`raptor prepare` to preprocess the files. `raptor prepare` accepts multiple paths per line, but only produces
one output file per bin, which can then be used with `raptor layout`.

### List of preprocessed files
\include{doc} fragments/input_files_preprocessed.md

## -​-kmer-size
See \ref usage_w_vs_k.
\note
This parameter will be used by `raptor build` and hence should be chosen carefully. However, overwriting it in the
`raptor build` call is possible but not recommended, and will result in a warning being emitted.

## -​-num-hash-functions
\include{doc} fragments/ibf_h.md
\note
This parameter will be used by `raptor build` and hence should be chosen carefully. However, overwriting it in the
`raptor build` call is possible but not recommended, and will result in a warning being emitted.

## -​-false-positive-rate
\include{doc} fragments/ibf_fpr.md
\note
This parameter will be used by `raptor build` and hence should be chosen carefully. However, overwriting it in the
`raptor build` call is possible but not recommended, and will result in a warning being emitted.

## -​-output-filename
The output filename may be freely chosen. `raptor build` does not rely on the file name or extension.
The default is `binning.out`.

## -​-threads
The number of threads to use. Currently, the layout algorithm is only partly parallelized.

\note
This option has no effect if `--disable-rearrangement` is set.

## -​-disable-estimate-union
The layout algorithm estimates the sequence similarity between input data. This improves the quality of the resulting
layout at the expense of a higher RAM requirement.

This behaviour can be disabled with this flag.

\note
Setting this flag also sets `--disable-rearrangement`.

## -​-disable-rearrangement
Based on the estimated sequence similarity, the order of sequences is changed to improve layout quality.
The more input sequences there are, the longer this process takes.

This behaviour can be disabled with this flag.

## Advanced Parameters

## -​-tmax
Limits the number of technical bins on each level of the HIBF. Choosing a good tmax is not trivial. The
smaller tmax, the more levels the layout needs to represent the data. This results in a higher space
consumption of the index. While querying each individual level is cheap, querying many levels might also
lead to an increased runtime. A good tmax is usually the square root of the number of user bins/samples
rounded to the next multiple of 64. Note that your tmax will always be rounded to the next multiple of 64.
At the expense of a longer runtime, you can enable the statistic mode that determines the best tmax for your
data set. See the advanced option --determine-best-tmax Default: ≈sqrt(samples).

## -​-alpha
The layout algorithm optimizes the space consumption of the resulting HIBF but currently has no means of
optimizing the runtime for querying such an HIBF. In general, the ratio of merged bins and split bins
influences the query time because a merged bin always triggers another search on a lower level. To influence
this ratio, alpha can be used. The higher alpha, the less merged bins are chosen in the layout. This
improves query times but leads to a bigger index. Default: 1.2.

## -​-max-rearrangement-ratio
When the flag --disable-rearrangement is not set, this option can influence the rearrangement algorithm. The
algorithm only rearranges the order of user bins in fixed intervals. The higher --max-rearrangement-ratio,
the larger the intervals. This potentially improves the layout, but increases the runtime of the layout
algorithm. Default: 0.5. Value must be in range [0.000000,1.000000].

## -​-sketch-bits
The number of bits the HyperLogLog sketch should use to distribute the values into bins. Default: 12. Value
must be in range [5,32].

## -​-determine-best-tmax
When this flag is set, the program will compute multiple layouts for tmax in [64 , 128, 256, ... , tmax] as
well as tmax=sqrt(samples). The layout algorithm itself only optimizes the space consumption. When
determining the best layout, we additionally keep track of the average number of queries needed to traverse
each layout. This query cost is taken into account when determining the best tmax for your data. Note that
the option --tmax serves as upper bound. Once the layout quality starts dropping, the computation is
stopped. To run all layout computations, pass the flag --force-all-binnings.

## -​-force-all-binnings
Forces all layouts up to --tmax to be computed, regardless of the layout quality. If the flag
--determine-best-tmax is not set, this flag is ignored and has no effect.

## -​-output-sketches-to
If you supply a directory path with this option, the hyperloglog sketches of your input will be stored in
the respective path; one .hll file per input file. Default: None.

