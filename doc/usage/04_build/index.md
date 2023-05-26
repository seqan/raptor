# raptor build {#usage_build}

[TOC]

## Main Parameters

### -​-input

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  A layout file produced with `raptor layout`.

- <b class="tab-title">IBF</b>
  The input file contains paths to the sequence data. Each line may contain multiple paths (separated by a whitespace).

  ```
  /absolute/path/to/file1.fasta /absolute/path/to/file2.fasta
  /absolute/path/to/file3.fa.gz
  ```

  Preprocessed files (See `raptor prepare`) are also supported.

  Supported file extensions are (possibly followed by bz2, gz, or bgzf):<br>
    &nbsp;&nbsp;• embl<br>
    &nbsp;&nbsp;• fasta<br>
    &nbsp;&nbsp;• fa<br>
    &nbsp;&nbsp;• fna<br>
    &nbsp;&nbsp;• ffn<br>
    &nbsp;&nbsp;• faa<br>
    &nbsp;&nbsp;• frn<br>
    &nbsp;&nbsp;• fas<br>
    &nbsp;&nbsp;• fastq<br>
    &nbsp;&nbsp;• fq<br>
    &nbsp;&nbsp;• genbank<br>
    &nbsp;&nbsp;• gb<br>
    &nbsp;&nbsp;• gbk<br>
    &nbsp;&nbsp;• sam

</div>

### -​-output
The output file name.

### -​-threads
The number of threads to use. Both IBF and HIBF construction can heavily benefit from parallisation.

### -​-quiet
By default, runtime and memory statistics are printed at the end.

This flag disables this behaviour.

### -​-kmer

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `k`.

- <b class="tab-title">IBF</b>
  Depending on the number of errors that should be accounted for when searching, the `kmer-size` (`k`) has to be chosen
  such that the k-mer lemma still has a positive threshold.

  **K-mer counting lemma**: For a given `k` and number of errors `e`, there are \f$k_p = |p| - k + 1\f$ many k-mers in the
  pattern `p` and an approximate occurrence of `p` in text `T` has to share at least \f$t = (k_p - k \cdot e)\f$ k-mers.

  For example, when searching reads of length 100 and allowing 4 errors, k has to be at most 20
  (100 − 20 + 1 − 4 · 20 = 1).

  Furthermore, k shall be such that a random k-mer match in the database is unlikely.
  For example, we chose k = 32 for the RefSeq data set. In general, there is no drawback in
  choosing the (currently supported) maximum k of 32, as long as the aforementioned
  requirements are fulfilled.

</div>

### -​-window

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `w`.

- <b class="tab-title">IBF</b>
  In the case that minimizers should be used, w has to be bigger than k.
  Canonical k-mers are achieved with w = k. Depending on the choice of k, the choice
  of w has to be made such that we obtain positive thresholds with our probabilistic
  threshold. In general, we aim at having a minimum threshold of 3 for the (w, k) minimizers. Hence, we choose w as large as possible, such that the minimum threshold is 3.
  For example, this is obtained for 2 errors and read length 150 for (29, 20) minimizers, for read length 100 for (24, 20)-minimizers and for read length 250 for (40, 20)
  minimizers, which in turn reduces the amount of k-mers by factors of 5.5, 3, and 11,
  respectively

</div>

### -​-shape
Sets a shape to use. This allows using gaps (do-not-care-positions) for k-mers.

\note
Mutually exclusive with `--kmer`.

### -​-fpr

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `fpr`.

- <b class="tab-title">IBF</b>
  Sets an uper bound for Bloom Filter false positives.
  A low `fpr` limits the number of false-positive results that have to be handled downstream.
  A higher `fpr` can help to reduce memory consumption in cases where false-positive k-mers have little effect.

</div>

[Bloom Filter Calculator](https://hur.st/bloomfilter/)

### -​-hash
The number of hash functions to use for Bloom Filters. Influences the index size.
[Bloom Filter Calculator](https://hur.st/bloomfilter/)

### -​-parts

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Not yet available.

- <b class="tab-title">IBF</b>
  Splits the index. Lower memory consumption, but higher runtime cost when searching.
  Output files will have a suffix `_x`, where `x` is a number.

</div>

