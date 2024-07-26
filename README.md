<!--
SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Raptor [![build status][1]][2] [![codecov][3]][4] [![install with bioconda][5]][6]

[1]: https://img.shields.io/github/actions/workflow/status/seqan/raptor/ci_linux.yml?branch=main&style=flat&logo=github&label=Raptor%20CI
[2]: https://github.com/seqan/raptor/actions?query=branch%3Amain
[3]: https://codecov.io/gh/seqan/raptor/branch/main/graph/badge.svg?token=SJVMYRUKW2
[4]: https://codecov.io/gh/seqan/raptor
[5]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[6]: #install-with-bioconda-linux

### A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences

## Download and Installation
See the [documentation](https://docs.seqan.de/raptor/main/setup.html) for detailed information on how to download,
install, and compile Raptor.

There may be performance benefits when compiling from source, as Raptor can be optimized for the host system.

### Install with [bioconda](https://bioconda.github.io/recipes/raptor/README.html) (Linux)

```bash
conda install -c bioconda -c conda-forge raptor
```

## Example Data and Usage
See the [documentation](https://docs.seqan.de/raptor/main/usage_quickstart.html) for usage.

`raptor --help` will show available commands. Each command will have a respective help page that can be shown via, e.g.,
`raptor build --help`.

## Authorship and Copyright
Raptor is being developed by [Enrico Seiler](mailto:enrico.seiler@fu-berlin.de), but also incorporates much work from
other members of [SeqAn](https://www.seqan.de).

### Citation
In your academic works (also comparisons and pipelines) please cite:
  * *Raptor: A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences*;
    Enrico Seiler, Svenja Mehringer, Mitra Darvish, Etienne Turc, and Knut Reinert; iScience 2021 24 (7): 102782.
    doi: https://doi.org/10.1016/j.isci.2021.102782
  * *Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries*;
    Svenja Mehringer, Enrico Seiler, Felix Droop, Mitra Darvish, René Rahn, Martin Vingron, and Knut Reinert;
    Genome Biol 24, 131 (2023). doi: https://doi.org/10.1186/s13059-023-02971-4

### RECOMB 2021
Raptor was presented at the [25th International Conference on Research in Computational Molecular Biology][recomb_url]:
  * [Slides][recomb_slides]
  * [Pre-recorded version of the talk][recomb_talk]

Please see the [License](#license) section for information on allowed use.

[recomb_url]: https://www.recomb2021.org/
[recomb_slides]: https://box.fu-berlin.de/s/TtM3Raxixm35Syy
[recomb_talk]: https://box.fu-berlin.de/s/YJFQnwqdE5q2Tym

### Supplementary
The subdirectory `util` contains applications and scripts related to papers.

### License
Raptor is open-source software. However, certain conditions apply when you (re-)distribute and/or modify Raptor,
please see the [license](https://github.com/seqan/raptor/blob/main/LICENSE.md).
