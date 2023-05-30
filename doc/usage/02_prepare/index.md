# raptor prepare {#usage_prepare}

[TOC]

Optionally preprocesses files for the use with `raptor layout` and `raptor build`.

Can continue where it left off after a crash or in multiple runs.

When to use:
  * Need to apply filtering
  * ...

## Main Parameters

### -​-input
The input file contains paths to the sequence data. Each line may contain multiple paths (separated by a whitespace).

```
/absolute/path/to/file1.fasta /absolute/path/to/file2.fasta
/absolute/path/to/file3.fa.gz
```

Preprocessed files (See `raptor prepare`) are also supported.

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

### -​-output {#usage_prepare_output}
A path to the output directory. The directory will be created if it does not exist.

Will create a `minimiser.list` inside the output directory. This file contains a list of generated minimiser
files, in the same order as the input.
This file can be used as input for `raptor layout` or `raptor build`.

Created output files for each input file:
  *.header: Contains the shape, window size, cutoff and minimiser count.
  *.minimiser: Contains binary minimiser values, one minimiser per line.
  *.in_progress: Temporary file to track process. Deleted after finishing computation.

\note
If `raptor prepare` aborts unexpectedly, you can rerun the same command. Files that have already preprocessed will
be skipped.

\attention
When you manually delete a .in_progress file, also delete the corresponding .header and .minimiser file.

### -​-threads
The number of threads to use. Multiple files will be handled in parallel. While more threads speed up the
preprocessing, the RAM usage also increases.

\note
Use less threads if `raptor prepare` fails due to RAM restrictions.

### -​-quiet
By default, runtime and memory statistics are printed to stderr at the end.

This flag disables this behaviour.

### -​-kmer
### -​-window
### -​-shape

### -​-kmer-count-cutoff
Only store k-mers with at least (>=) x occurrences.

\note
Mutually exclusive with --use-filesize-dependent-cutoff.

### -​-use-filesize-dependent-cutoff
Apply cutoffs from Mantis(Pandey et al., 2018).

| File size | Cutoff |
|-----------|--------|
| ≤ 300 MiB | 1      |
| ≤ 500 MiB | 3      |
| ≤ 1 GiB   | 10     |
| ≤ 3 GiB   | 20     |
| > 3 GiB   | 50     |

File sizes are based of gzipped FASTQ files. Compression reduces the filesize by around factor `3`. FASTA files are
approximately `2` times smaller than FASTQ.

\note
Mutually exclusive with --kmer-count-cutoff.
