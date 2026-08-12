# raptor update {#usage_update}

<!--
SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

`raptor update` adds user bins to an existing HIBF, or removes them from it, without rebuilding the index from all
input data. This is useful when a collection grows or shrinks over time, e.g. when new samples are added to a
database periodically.

\attention
An index can only be updated if it was prepared for it. See \ref usage_update_requirements.

# Preparing an updatable index {#usage_update_requirements}

Updating needs spare technical bins to place new data into. `raptor layout` only reserves them if you ask for it:

```bash
raptor layout --input all_bins.txt --output layout.txt --empty-bin-fraction 0.1
raptor build --input layout.txt --output raptor.index
```

`--empty-bin-fraction` sets the fraction of `tmax` that is kept empty in every IBF of the layout. The reserved bins
cost space in the index, but they are what an insertion writes into, so choosing `0` gives an index that cannot be
updated at all:

```
[Error] The index does not track the occupancy of its technical bins and hence cannot be updated. Recompute the
layout using 'raptor layout' with --empty-bin-fraction greater than 0, and rebuild the index.
```

There is no way to retrofit this onto an existing index; the layout has to be recomputed and the index rebuilt.

\note
Only an HIBF can be updated, so a layout is always required (\ref hibf_vs_ibf).

The larger the fraction, the more insertions the index absorbs before it has to be rebuilt, and the larger it is to
begin with. `0.1` is a reasonable starting point.

# Workflow

`raptor update` never modifies its input. It reads the index given by `--index` and writes a new one to `--output`,
which must not exist yet. Chaining updates therefore means passing the previous output as the next input:

```bash
raptor update insert --index raptor.index   --output raptor.1.index --insert batch1.txt
raptor update insert --index raptor.1.index --output raptor.2.index --insert batch2.txt
```

An updated index is searched exactly like any other index:

```bash
raptor search --index raptor.2.index --query query.fq --output search.out --error 2
```

# raptor update insert

```bash
raptor update insert --index raptor.index --output updated.index --insert new_bins.txt
```

## -​-insert

Either a single sequence file, which is added as one user bin, or a file listing the files to add, one user bin per
line:

\include{doc} fragments/input_files_sequence.md

Preprocessed input is supported as well:

\include{doc} fragments/input_files_preprocessed.md

The index and the files being inserted do not have to be of the same type: minimiser files can be inserted into an
index built from sequence files and vice versa, as long as they were preprocessed with the same `--kmer` and
`--window` values.

\note
Unlike `raptor build`, insertion currently supports only **one file per user bin**. A line holding several paths is
rejected with `Currently not supporting multiple files per UB for insert.`

## -​-no-partial-rebuild

Inserting degrades an index over time: bins fill up beyond their target false positive rate, and IBFs grow past the
`tmax` their layout assumed. Raptor repairs this automatically, and there are two ways to do so:

* A **partial rebuild** recomputes the layout of a single subtree of the HIBF. It only has to process the user bins
  of that subtree, so it is cheap, but it can only improve that part of the index.
* A **full rebuild** recomputes the layout of the whole index from all user bins. It yields the better layout, but
  has to process everything.

By default, the cheapest repair that fixes the problem is chosen. This flag forces a full rebuild whenever a repair
is needed at all. It trades update time for index quality, which can pay off if you insert rarely and search often.

\note
A full rebuild processes *all* user bins, including the ones of the current call that have not been inserted yet.
Once it happens, the remaining insertions of that call are already covered by it.

# raptor update delete

```bash
raptor update delete --index raptor.index --output updated.index --delete 0 --delete 3
```

## -​-delete

The IDs of the user bins to remove. A user bin ID is the position of its entry in the input file that the index was
built from, starting at `0`. The IDs are also listed in the header of a `raptor search` result:

```
#0	/path/to/bin1.fa
#1	/path/to/bin2.fa
#2	/path/to/bin3.fa
#QUERY_NAME	USER_BINS
query1	0,2
```

Pass the option once per user bin.

## Effect on user bin IDs

Deleting does **not** renumber the remaining user bins: their IDs stay valid across updates, which is what makes it
safe to keep a mapping of your own. The path of a deleted user bin therefore still appears in the header of a search
result, but the bin is never reported as a hit:

```
#0	/path/to/bin1.fa
#1	/path/to/bin2.fa    <-- deleted, never reported
#2	/path/to/bin3.fa
#QUERY_NAME	USER_BINS
query1	0,2
```

The technical bins that a deleted user bin occupied are reclaimed and reused by later insertions, so deleting and
re-inserting does not permanently grow the index.

# Limitations

* The index must have been built from a layout with `--empty-bin-fraction` greater than `0`
  (\ref usage_update_requirements).
* Only the HIBF is supported.
* Insertion supports one file per user bin.
* `--kmer`, `--window`, `--fpr` and `--hash` are taken from the index and cannot be changed by an update.
* Inserting many user bins into an index is not as good as computing a layout for all of them at once. If the index
  has grown far beyond its original size, rebuilding it with `raptor layout` and `raptor build` may produce a
  smaller and faster index.
