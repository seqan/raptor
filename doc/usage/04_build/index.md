<!--
SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# raptor build {#usage_build}

[TOC]

# Main Parameters

## -​-input

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  A layout file produced with `raptor layout`.

- <b class="tab-title">IBF</b>
  <span style="font-weight: 600; font-size: 1.17em;">List of sequence files</span><br>
  \include{doc} fragments/input_files_sequence.md
  <span style="font-weight: 600; font-size: 1.17em;">List of preprocessed files</span><br>
  \include{doc} fragments/input_files_preprocessed.md

</div>

## -​-output
The output file name.

## -​-threads
The number of threads to use. Both IBF and HIBF construction can heavily benefit from parallelisation.

## -​-quiet
By default, runtime and memory statistics are printed to stderr at the end.

This flag disables this behaviour.

## -​-kmer

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `k`.

- <b class="tab-title">IBF</b>
  See \ref usage_w_vs_k.

</div>

## -​-window

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `w`.

- <b class="tab-title">IBF</b>
  See \ref usage_w_vs_k.

</div>

## -​-fpr

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `fpr`.

- <b class="tab-title">IBF</b>
  \include{doc} fragments/ibf_fpr.md

</div>

[Bloom Filter Calculator](https://hur.st/bloomfilter/)

## -​-hash
<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Read from layout file. The read value can be overwritten with this option. However, a warning will be emitted to prevent
  accidentally overwriting `h`.

- <b class="tab-title">IBF</b>
  \include{doc} fragments/ibf_h.md

</div>

## -​-parts

<div class="tabbed">

- <b class="tab-title">HIBF</b>
  Not yet available.

- <b class="tab-title">IBF</b>
  Splits the index. Lower memory consumption, but higher runtime cost when searching.
  Output files will have a suffix `_x`, where `x` is a number.

</div>

