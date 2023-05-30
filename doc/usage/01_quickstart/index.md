# Quickstart {#usage_quickstart}

[TOC]

<!-- * reference usage sections
* copy/paste examples with many/all options
* present subcommands
* present workflow ibf/hibf
* present (w,k) minimiser, (k,k) canonical k-mers -->

## What problems can be solved

* Approximate Membership Query (AMQ)
* Input sequences = bins
* Given a query sequence: Show all bins that possible contain that sequence.
  * Either given a number of errors
  * Or given a percentage (e.g., 0.7 -> 70% of minimiser/k-mers of query must be present in bin)
* Does not model a colored de Bruijn graph:
  * ✔️ Given a k-mer: Tell to which color (bin) it belongs to
  * ❌ Given a color (bin): Show all k-mers that belong to this color (bin)

## HIBF vs IBF
Recommendation: HIBF

Cases in which to consider the IBF:
  * Evenly sized bins
  * Small number of bins (≤ 128)

<details><summary>Click to see an HIBF / IBF comparison.</summary>
\htmlonly
<!--
    Including an SVG like this allows using CSS variables inside the SVG, e.g.,
    fill="var(--page-foreground-color)"
    fig_4_many_user_bins.svg#svg refers to id="svg" in the svg file
    viewBox of svg element is copied from svg file
    svg file needs to be added to HTML_EXTRA_FILES
-->
<svg style="height: auto; width: 100%;" viewBox="0 0 342.22876 198.804504">
  <use href="/fig_4_many_user_bins.svg#svg" width="100%" height="100%"></use>
</svg>
\endhtmlonly
\note
The used data set is the **worst case** for the HIBF. In reality, the index size is usually smaller than the
corresponding IBF, and build times of the HIBF are much closer to IBF build times.

## Workflow

<div class="tabbed">

- <b class="tab-title">HIBF</b>
[`raptor prepare`] -> `raptor layout` -> `raptor build` -> `raptor search`

- <b class="tab-title">IBF</b>
[`raptor prepare`] -> `raptor build` -> `raptor search`

</div>

## Input files

### List of sequence files
Applies to:
  * `raptor prepare`
  * `raptor layout` (without running `raptor prepare` beforehand)
  * `raptor build` (without running `raptor prepare` beforehand)

The input file contains paths to the sequence data. Each line may contain multiple paths (separated by a whitespace).

\attention
`raptor layout` currently does not accept multiple paths, see \ref usage_layout_input for a workaround.

```
/absolute/path/to/file1.fasta /absolute/path/to/file2.fasta
/absolute/path/to/file3.fa.gz
```

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

### List of preprocessed files
Applies to:
  * `raptor layout`
  * `raptor build`

See \ref usage_prepare_output.

\note
As a convenience, `raptor prepare` creates a `minimiser.list` file in the output directory. This file can be used as
input for `raptor layout` and `raptor build`.

\attention
The window and kmer sized used for preprocessing are propagated to `raptor layout` and `raptor build` and cannot be
overwritten there.

## Choosing window and kmer size

### (w,k) minimiser vs (k,k) canonical k-mers

|               | (k,k) | (w,k) |
|---------------|-------|-------|
| Index size    | ❌   | ✔️   |
| Runtime       | ❌   | ✔️   |
| RAM usage     | ❌   | ✔️   |
| Thresholding¹ | Exact | Heuristic |
<small>¹ When searching with a given number of errors.</small><br>

* (w,k) minimiser reduces the number of values to process by roughly \f$\frac{w - k + 2}{2}\f$.
* (w,k) minimiser have a slightly lower accuracy than (k,k). However, the loss of accuracy mainly stems from false
  positves, not false negatives.

#### Recommendation {#usage_w_vs_k_figure}
* (w,k) with gentle compression (w-k=4)

<details><summary>Click to see the differences of (w,k) and (k,k) on different aspects.</summary>
\htmlonly
<svg style="height: auto; width: 100%;" viewBox="0 0 307.29294 196.59375">
  <use href="/fig_6_acc.svg#svg" width="100%" height="100%"></use>
</svg>
\endhtmlonly
</details>

### (w,k)
Requirements:
  * `w > k`
  * `w ≤ query length`

Recommendation:
  * `w - k = 4`
  * `w << query length`

Examples:
  * query length 100: `w = 24`, `k = 20`
  * query length 250: `w = 28`, `k = 24`

Also see the figure in \ref usage_w_vs_k_figure.

### (k,k)
Requirements:
  * `k ≤ query length`
  * k-mer counting lemma satisfied, when searching with a given number of errors.

Recommendation:
  * `k << query length`

Examples (for two errors):
  * query length 100: `k = 20`
  * query length 250: `k = 32`

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

## Choosing IBF parameters

### h
Recommendation: default value (`2`)

### fpr
Recommendation: default value (`0.05`)

* A lower `fpr` limits the number of false-positive results, but increases index size.
* A higher `fpr` can help to reduce memory consumption in cases where false-positive k-mers have little effect.

