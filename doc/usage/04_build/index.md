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
By default, runtime and memory statistics are printed to stderr at the end.

This flag disables this behaviour.

### -​-kmer

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `k`.

- <b class="tab-title">IBF</b>
  blah

</div>

### -​-window

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `w`.

- <b class="tab-title">IBF</b>
  Blah

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

