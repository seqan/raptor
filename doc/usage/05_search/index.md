# raptor search {#usage_search}

<!--
SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

## Main Parameters

### -​-index
The path to the index. For partitioned indices, the suffix `_x`, where `x` is a number, must be omitted.

### -​-query
File containing query sequences.

<details><summary>Many file types and compressions are supported. Click to show a list.</summary>
Supported file extensions are (possibly followed by bz2, gz, or bgzf):
  * embl
  * fasta
  * fa
  * fna
  * ffn
  * faa
  * frn
  * fas
  * fastq
  * fq
  * genbank
  * gb
  * gbk
  * sam
</details>

### -​-output
The output file name.

<div class="tabbed">

- <b class="tab-title">Format</b>
  ```
  ###<text>                     | Meta-information
  ##<text>                      | Meta-information
  #<number><tab><filepaths>     | Assigns each input file a number. Multiple filepaths are separated by a whitespace
  #QUERY_NAME<tab>USER_BINS     | Header for the results
  <query_id><tab>[<number>...]  | A line for each query, listing matches in input files, if any. Multiple hits are separated by a comma.
  ```

- <b class="tab-title">Example</b>
  \include search_output.txt

</div>

### -​-threads
The number of threads to use. Sequences in the query file will be processed in parallel.
Negligible effect on RAM usage for unpartitioned indices. Moderate effect for partitioned indices.

### -​-quiet
By default, runtime and memory statistics are printed to stderr at the end.

This flag disables this behaviour.

### -​-error
The number of allowed errors.

\note
Mutually exclusive with --threshold.

### -​-threshold
Ratio of k-mers that need to be found for a hit to occur.

\note
Mutually exclusive with --error.

### -​-query_length
The sequence length of a query. Used to determine thresholds. The sequence lengths should have little to no variance.

If not provided:
  * the median of sequence lengths in the query file is used.
  * a warning is emitted if there is a high variance in sequence lengths.
  * an error occurs if any sequence is shorter than the window size.

### -​-tau
The higher tau, the lower the threshold.

\note
Has no effect when using `--threshold` or `w` == `k`.

### -​-p_max
The higher p_max, the higher the threshold.

\note
Has no effect when using `--threshold` or `w` == `k`.

### -​-cache-thresholds
Stores the computed thresholds with a unique name next to the index. In the next search call using this
option, the stored thresholds are re-used.
Two files are stored:
  * threshold_*.bin: Depends on query_length, window, kmer/shape, errors, and tau.
  * correction_*.bin: Depends on query_length, window, kmer/shape, p_max, and fpr.

\note
Has no effect when using `--threshold` or `w` == `k`.
