<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Dynamic HIBF benchmarks

## Reports

The interactive reports are available at https://docs.seqan.de/dynamic.

The `reports` directory contains the same content. See its [README](reports/README.md) for how to view it locally.

## Zenodo

You can find raw results, input data information, and the simulated datasets at https://zenodo.org/records/21689968.

## Benchmark workflow

The `dhibf_benchmark` directory is a git subtree and includes the content of https://gitlab.com/pirovc/dhibf_benchmark.

For a general overview regarding `git subtree`, see, for example, https://www.atlassian.com/git/tutorials/git-subtree.

Pulling upstream changes can be achieved by:

```bash
git subtree pull --prefix util/dynamic/dhibf_benchmark https://gitlab.com/pirovc/dhibf_benchmark.git main --squash
```

Similarly, you can push subtree changes to your fork:

```bash
git subtree push --prefix util/dynamic/dhibf_benchmark https://(github.com|gitlab.com)/<user>/dhibf_benchmark.git <branch>
```
