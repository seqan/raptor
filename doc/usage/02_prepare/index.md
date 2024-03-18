<!--
SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# raptor prepare {#usage_prepare}

[TOC]

Optionally preprocesses files for the use with `raptor layout` and `raptor build`.

Can continue where it left off after a crash or in multiple runs.

When to use:
  * Applying k-mer filtering based on abundance.

# Main Parameters

## -​-input

\include{doc} fragments/input_files_sequence.md

## -​-output {#usage_prepare_output}
A path to the output directory. The directory will be created if it does not exist.

Will create a `minimiser.list` inside the output directory. This file contains a list of generated minimiser
files, in the same order as the input.
This file can be used as input for `raptor layout` or `raptor build`.

Created output files for each input file:
  * `*.header`: Contains the shape, window size, cutoff and minimiser count.
  * `*.minimiser`: Contains binary minimiser values, one minimiser per line.
  * `*.in_progress`: Temporary file to track process. Deleted after finishing computation.

\attention
The window and k-mer sized used for preprocessing are propagated to `raptor layout` and `raptor build` and cannot be
overwritten there.

\note
If `raptor prepare` aborts unexpectedly, you can rerun the same command. Files that have already preprocessed will
be skipped.

\attention
When you manually delete a `.in_progress` file, also delete the corresponding `.header` and `.minimiser` file.

## -​-threads
The number of threads to use. Multiple files will be handled in parallel. While more threads speed up the
preprocessing, the RAM usage also increases.

\note
Use less threads if `raptor prepare` fails due to RAM restrictions.

## -​-quiet
By default, runtime and memory statistics are printed to stderr at the end.

This flag disables this behaviour.

## -​-kmer
See \ref usage_w_vs_k.
\attention
This parameter will be used by `raptor build` and hence should be chosen carefully. The k-mer size cannot be changed
afterwards.

## -​-window
See \ref usage_w_vs_k.
\attention
This parameter will be used by `raptor build` and hence should be chosen carefully. The window size cannot be changed
afterwards.

## -​-kmer-count-cutoff
Only store k-mers with at least (>=) x occurrences.

\note
Mutually exclusive with --use-filesize-dependent-cutoff.

## -​-use-filesize-dependent-cutoff
Apply cutoffs from Mantis(Pandey et al., 2018).

| File size | Cutoff |
|-----------|--------|
| ≤ 300 MiB | 1      |
| ≤ 500 MiB | 3      |
| ≤ 1 GiB   | 10     |
| ≤ 3 GiB   | 20     |
| > 3 GiB   | 50     |

File sizes are based of gzipped FASTQ files. Compression reduces the file size by around factor `3`. FASTA files are
approximately `2` times smaller than FASTQ.

\note
Mutually exclusive with --kmer-count-cutoff.
