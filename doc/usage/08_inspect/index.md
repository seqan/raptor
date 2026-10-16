# raptor inspect {#usage_inspect}

<!--
SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

`raptor inspect` prints the metadata of an index: the parameters it was built with, the layout parameters of an HIBF,
and how far updates have changed it.

```bash
raptor inspect --index raptor.index
```

By default, only the parameters and the file names stored at the beginning of the index file are read. Inspecting is
therefore fast even for very large indices. `--detailed` loads the whole index to also report the structure of its
IBFs.

# Main parameters

## -​-index

The path to the index. For partitioned indices, the suffix `_x`, where `x` is a number, must be omitted.

## -​-output

Write the metadata to this file instead of the standard output. An existing file is overwritten.

## -​-format

Either `text` (default) or `json`. Both contain the same values, see \ref usage_inspect_text and
\ref usage_inspect_json.

## -​-user-bins

Also list the files of every user bin, one line per user bin. For a large collection, this list is long.

## -​-detailed

Also report the structure of the index, see \ref usage_inspect_structure. This loads the whole index and therefore takes
as long and as much memory as loading the index for `raptor search`. A partitioned index is loaded one part at a time.

# Text output {#usage_inspect_text}

The text output uses the same `### Section` and `## Key = value` lines as the header of a `raptor search` result.
For an HIBF built with

```bash
raptor layout --input all_bins.txt --output layout.txt --kmer 19 --empty-bin-fraction 0.1
raptor build --input layout.txt --output raptor.index --window 23
raptor inspect --index raptor.index
```

the output is:

```
### Index
## Index = "raptor.index"
## Index version = 3
## Index type = HIBF
## Index parts = 1
## Index size [bytes] = 4527
### Minimiser parameters
## Window size = 23
## Shape = 1111111111111111111
## Shape size (length) = 19
## Shape count (number of 1s) = 19
### Bloom filter parameters
## False positive rate = 0.05
## Relaxed false positive rate = 0.3
## Hash functions = 2
### Layout parameters
## Number of user bins in layout = 4
## tmax = 64
## Empty bin fraction = 0.109375
## Sketch bits = 12
## alpha = 1.2
## Max rearrangement ratio = 0.5
## Disable estimate union = false
## Disable rearrangement = false
### Update state
## Track occupancy = true
## Number of IBFs = 1
## Original number of IBFs = 1
## Resized IBFs = 0
### User bins
## Number of user bins = 4
```

| Section                 | Content                                                                                       |
|-------------------------|-----------------------------------------------------------------------------------------------|
| Index                   | The path as passed to `--index`, the version of the index format, whether it is an IBF or an HIBF, the number of parts, and the combined size of all its files. |
| Minimiser parameters    | The window size and the shape, as set by `--window` and `--kmer`/`--shape` (\ref usage_w_vs_k). |
| Bloom filter parameters | The false positive rate (`--fpr`). For an HIBF, also the relaxed false positive rate of merged bins (`--relaxed-fpr`) and the number of hash functions (`--hash`). |
| Layout parameters       | HIBF only. The parameters of the layout the index was built from, see \ref usage_layout.      |
| Update state            | HIBF only. How far updates have changed the index, see below.                                 |
| User bins               | The number of user bins, including deleted ones.                                              |

\note
For an IBF, the number of hash functions is stored with the IBF data. It is only reported with `--detailed`.

## Layout parameters

`Number of user bins in layout` is the number of user bins the layout was computed for. Inserting user bins
(\ref usage_update) does not change it, unless the insertion triggers a full rebuild, which computes a new layout from
all user bins.

`Empty bin fraction` is the fraction actually reserved, which is rounded to whole technical bins of `tmax`. In the
example, `--empty-bin-fraction 0.1` of `tmax = 64` reserves 7 of 64 technical bins, i.e. `0.109375`.

## Update state

| Key                     | Meaning                                                                                       |
|-------------------------|-----------------------------------------------------------------------------------------------|
| Track occupancy         | Whether the index tracks the occupancy of its technical bins. Only then can it be updated, see \ref usage_update_requirements. |
| Number of IBFs          | The number of IBFs, including IBFs that updates have retired.                                 |
| Original number of IBFs | The number of IBFs when the index was last built from a layout, by `raptor build` or by a full rebuild. |
| Resized IBFs            | The number of IBFs that insertions have grown beyond the size their layout assigned to them. A partial rebuild counts the IBFs it rebuilt and added as resized. |

The further `Number of IBFs` and `Resized IBFs` have grown beyond `Original number of IBFs`, the more the index has
drifted from a layout computed for its current content.

## User bins

With `--user-bins`, the output additionally lists the files of every user bin, directly after
`## Number of user bins`. The format is the same as in the header of a `raptor search` result:

```
#0	bins/bin1.fa
#1	bins/bin2.fa
#2	bins/bin3.fa
#3	bins/bin4.fa
```

A user bin ID is the number after `#`. It is the ID that `raptor search` reports and that `raptor update delete`
expects.

## Structure {#usage_inspect_structure}

With `--detailed`, the output additionally contains a `### Structure` section and one line per IBF. After deleting
user bin 1 from the example index with `raptor update delete --index raptor.index --output updated.index --delete 1`,
`raptor inspect --index updated.index --detailed --user-bins` reports, from `### User bins` on:

```
### User bins
## Number of user bins = 4
#0	bins/bin1.fa
#1	bins/bin2.fa
#2	bins/bin3.fa
#3	bins/bin4.fa
### Structure
## Levels = 1
## Technical bins = 57
## Size [bits] = 21568
## Number of retired IBFs = 0
## Number of deleted user bins = 1
## Deleted user bins = 1
#IBF	LEVEL	PARENT_IBF	PARENT_BIN	TECHNICAL_BINS	CAPACITY	MERGED_BINS	DELETED_BINS	BIN_SIZE	HASH_FUNCTIONS	BITS
0	0	-	-	57	64	0	15	337	2	21568
```

| Key                         | Meaning                                                                                   |
|-----------------------------|-------------------------------------------------------------------------------------------|
| Levels                      | The number of levels of the HIBF. Always `1` for an IBF.                                  |
| Technical bins              | The number of technical bins in use, summed over all IBFs.                                |
| Size [bits]                 | The size of all IBFs. The index file additionally contains the metadata.                  |
| Number of retired IBFs      | HIBF only. Deleting every user bin of a lower-level IBF retires it: its slot is kept, but it holds no data. `Retired IBFs` lists them. |
| Number of deleted user bins | HIBF only. Deleted user bins keep their ID until a full rebuild drops them (\ref usage_update_ids), but are no longer part of the index and are never reported by `raptor search`. `Deleted user bins` lists them. |

`Retired IBFs` and `Deleted user bins` are omitted if there are none.

The tab-separated lines after `#IBF` describe one IBF each. Retired IBFs are not listed. For a partitioned IBF, there is
one line per part.

| Column         | Meaning                                                                                        |
|----------------|------------------------------------------------------------------------------------------------|
| IBF            | The position of the IBF in the HIBF, or the part of a partitioned IBF.                         |
| LEVEL          | The depth in the HIBF. `0` for the top-level IBF.                                              |
| PARENT_IBF     | The IBF holding the merged bin this IBF belongs to. `-` for the top-level IBF.                 |
| PARENT_BIN     | That merged bin. `-` for the top-level IBF.                                                    |
| TECHNICAL_BINS | The number of technical bins in use.                                                           |
| CAPACITY       | The number of technical bins the IBF can hold without growing: the empty bins reserved by `--empty-bin-fraction` are included, and the number is rounded up to a multiple of 64. |
| MERGED_BINS    | The number of technical bins that point to a lower-level IBF.                                  |
| DELETED_BINS   | The number of technical bins freed by deleting user bins. Later insertions reuse them.        |
| BIN_SIZE       | The size of a single technical bin in bits.                                                    |
| HASH_FUNCTIONS | The number of hash functions.                                                                  |
| BITS           | The size of the IBF in bits, i.e. `CAPACITY` × `BIN_SIZE`.                                     |

In the example, user bin 1 was split across 15 technical bins, which are now free.

# JSON output {#usage_inspect_json}

`--format json` writes the same values as a single JSON object, with one member per section of the text output.
`raptor inspect --index updated.index --detailed --format json` reports:

```json
{
  "index": {
    "path": "updated.index",
    "version": 3,
    "type": "HIBF",
    "parts": 1,
    "size_in_bytes": 4527
  },
  "minimiser": {
    "window_size": 23,
    "shape": "1111111111111111111",
    "shape_size": 19,
    "shape_count": 19
  },
  "bloom_filter": {
    "false_positive_rate": 0.05,
    "relaxed_false_positive_rate": 0.3,
    "hash_functions": 2
  },
  "layout": {
    "number_of_user_bins": 4,
    "tmax": 64,
    "empty_bin_fraction": 0.109375,
    "sketch_bits": 12,
    "alpha": 1.2,
    "max_rearrangement_ratio": 0.5,
    "disable_estimate_union": false,
    "disable_rearrangement": false
  },
  "update": {
    "track_occupancy": true,
    "number_of_ibfs": 1,
    "original_number_of_ibfs": 1,
    "resized_ibfs": 0
  },
  "user_bins": {
    "count": 4
  },
  "structure": {
    "levels": 1,
    "technical_bins": 57,
    "bits": 21568,
    "retired_ibfs": [],
    "deleted_user_bins": [1],
    "ibfs": [
      {"id": 0, "level": 0, "parent_ibf": null, "parent_bin": null, "technical_bins": 57, "capacity": 64, "merged_bins": 0, "deleted_bins": 15, "bin_size": 337, "hash_functions": 2, "bits": 21568}
    ]
  }
}
```

* With `--user-bins`, `user_bins` additionally contains `paths`, an array with one array of file names per user bin.
  The position in `paths` is the user bin ID.
* For an IBF, `layout`, `update`, `bloom_filter.relaxed_false_positive_rate`, `bloom_filter.hash_functions`,
  `structure.retired_ibfs`, and `structure.deleted_user_bins` are omitted.
* `parent_ibf` and `parent_bin` are `null` for the top-level IBF.
* File names are written as stored in the index. A file name that is not valid UTF-8 results in a JSON string that is
  not valid UTF-8 either.
