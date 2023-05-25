# Raptor [![build status][1]][2] [![codecov][3]][4] [![install with bioconda][5]][6] [![install with brew][7]][8]

[1]: https://img.shields.io/github/actions/workflow/status/seqan/raptor/ci_linux.yml?branch=main&style=flat&logo=github&label=Raptor%20CI
[2]: https://github.com/seqan/raptor/actions?query=branch%3Amain
[3]: https://codecov.io/gh/seqan/raptor/branch/main/graph/badge.svg?token=SJVMYRUKW2
[4]: https://codecov.io/gh/seqan/raptor
[5]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[6]: #install-with-bioconda-linux
[7]: https://img.shields.io/badge/install%20with-brew-brightgreen.svg?style=flat
[8]: #install-with-brew-linux-macos

### A fast and space-efficient pre-filter for querying very large collections of nucleotide sequences

## Download and Installation
See the [documentation](https://seqan-raptor.vercel.app/setup.html) for detailed information on how to download,
install, and compile Raptor.

There may be performance benefits when compiling from source, as Raptor can be optimized for the host system.

### Install with [bioconda](https://bioconda.github.io/recipes/raptor/README.html) (Linux)

```bash
conda install -c bioconda -c conda-forge raptor
```

### Install with [brew](https://brew.sh/) (Linux, macOS)

```bash
brew install brewsci/bio/raptor
```

## Example Data and Usage
See the [documentation](https://seqan-raptor.vercel.app/tutorial_first_steps.html) for examples and usage.

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

## Sponsorships

[![Vercel][vercel_badge]][vercel_website]

[vercel_badge]: https://raw.githubusercontent.com/seqan/raptor/main/test/documentation/.vercel/powered-by-vercel.svg "Powered by Vercel"
[vercel_website]: https://vercel.com/?utm_source=seqan&utm_campaign=oss

We are grateful to Vercel for building and hosting our documentation, and they even provide preview builds within our
pull requests. Check them out!
