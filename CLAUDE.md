<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC0-1.0
-->

# CLAUDE.md

Guidance for Claude Code when working in this repository. Assume no other instructions exist.

## Working style

Treat the user as an expert-level developer: expert Git and Linux administration, university-level statistics,
proficient in C++, Python, LaTeX and web development. Skip beginner explanations, boilerplate and foundational
concepts. Responses are dense, direct and technical. Focus on architecture, template metaprogramming and
implementation detail rather than basic logic. Do not simplify integrations or abstract away complexity. Use metric
units for any hardcoded physical measurement.

## Verification discipline

This section is not generic advice; every rule below comes from a concrete failure in this repository.

* **Verify, don't infer.** Prefer a compiled probe, an instrumented run or a numeric table over reasoning from
  source. Confirm standard-library and third-party semantics against the actual headers or a test program. Label
  every claim as verified or unverified; never blur the two.
  *Example:* `v = {}` on a `std::vector` selects `operator=(initializer_list)` and keeps the capacity. Only
  `v = std::vector<T>{}` releases it. This was asserted incorrectly until a three-line program was compiled.
* **Exit code 0 is not a pass.** Validate the artifact — search output, resulting index, file contents — not just
  that the command returned. Two `raptor update` tests passed their exit-code checks while asserting the wrong
  thing.
* **Test where the assumption breaks.** A change that depends on an invariant or an optional feature must be
  exercised in the configuration where that thing is absent. A fix relying on `occupancy` was validated only with
  `--empty-bin-fraction 0.1`, the one setting where occupancy is tracked; with the default `0` it silently
  destroyed 48 of 60 user bins.
* **Prove new code paths execute.** Instrument to confirm a branch is reached. Several fixes here were never hit by
  the test workload (deleted-bin reuse fired zero times; subtree retirement required deleting *all* user bins).
  Report unexercised code as unexercised rather than implying coverage.
* **Equivalence checks need a noise floor.** Before diffing outputs to show a change is behaviour-preserving,
  establish whether the tool is deterministic. `raptor update insert` is **not** deterministic across a full
  rebuild — the same input run twice yields different bytes, even with `--threads 1`. Byte comparison is valid only
  up to the first full rebuild; beyond it use a functional oracle (search hits plus the user-bin mapping).
* **Never silently alter tuned heuristics, thresholds or algorithmic policy.** Transcribe faithfully, state the
  concern, let the user decide. If a benchmark contradicts the rationale for a requested change, report the
  measurement before proceeding.
* **Respect the stated scope.** When the user points at a specific item, address that item only. Raise adjacent
  findings as separate suggestions instead of bundling them.

## Project

Raptor is a pre-filter for approximate membership queries over large collections of nucleotide sequences, built on
[seqan/hibf](https://github.com/seqan/hibf), [seqan3](https://github.com/seqan/seqan3),
[sharg](https://github.com/seqan/sharg-parser) and [chopper](https://github.com/seqan/chopper), all pulled via CPM
and pinned in `cmake/package-lock.cmake`.

Subcommands (`src/raptor.cpp`): `build`, `layout`, `prepare`, `search`, `update`, `upgrade`.

Two index types: the **IBF** (flat, `raptor build` straight from sequence files) and the **HIBF** (hierarchical,
requires a layout from `raptor layout`). Workflow: `[prepare] → [layout] → build → search`.

```
include/raptor/        public headers (raptor_index, argument parsing, file_reader, ...)
src/                   one directory per subcommand; private headers live next to the .cpp
test/unit/api/         API tests
test/unit/cli/         CLI tests, one directory per subcommand
test/documentation/    doxygen build and doc tests
test/data/             test fixtures, copied to <build>/data by test/cmake/add_local_data.cmake
doc/usage/NN_<name>/   user documentation, one index.md per page
```

## Build and test

C++23, CMake ≥ 3.25, GCC or Clang, Linux and macOS. `CMAKE_BUILD_TYPE` defaults to `Release`; prefer `Debug`
when working on correctness so the asserts are live.

```bash
# Binary only.
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j        # -> build/bin/raptor

# Unit and CLI tests (this is what CI runs).
cmake -S test/unit -B build_test -DCMAKE_BUILD_TYPE=Debug
cmake --build build_test -j
ctest --test-dir build_test -j --output-on-failure

# Documentation (needs doxygen).
cmake -S test/documentation -B build_doc
cmake --build build_doc -j
ctest --test-dir build_doc -j --output-on-failure    # doc_no_warnings_test fails on any doxygen warning
```

Dependencies are fetched by CPM. Set `CPM_SOURCE_CACHE` to a directory to share them across build directories; that
cache is also where to read the hibf, seqan3 and sharg headers when checking their semantics instead of guessing.

If an existing build directory fails to reconfigure, configure a fresh one rather than fighting it.

Test builds add `-Werror`; pass `-DRAPTOR_WITH_WERROR=OFF` to turn that off while iterating
(`test/raptor-test.cmake`). Verify changed translation units under `-Wall -Wextra -Werror -pedantic` before
reporting completion.

## Conventions

* **Formatting.** `clang-format` with the repository `.clang-format`; the config names `clang-format-18`, and
  versions 18–20 agree in practice. Run it on every changed `.cpp`/`.hpp` and verify with
  `clang-format-18 --dry-run -Werror <file>`. Do not run it on CMake files.
* **Licensing.** Every file starts with the two `SPDX-FileCopyrightText` lines and an `SPDX-License-Identifier`:
  `BSD-3-Clause` for code and `CMakeLists.txt`, `CC-BY-4.0` for documentation, `CC0-1.0` for tooling config
  (`.clang-format`, CI workflows, `test/data/datasources.cmake`). Checked by `reuse` in CI.
* **Doxygen.** File header `/*!\file \brief ... \author ... */`. Members use `//!\brief` or
  `/*!\brief ... \details ... */`. Several `\brief` lines in `src/update/` were copy-paste wrong; check that a
  file's brief names what the file actually provides.
* **Asserts** are bare — `assert(cond);` — with the rationale in an adjacent comment. Exactly one assert in the
  repository carries a `&& "message"`; do not add more of that form.
* **Style.** `auto const` / `size_t const` for locals; designated initialisers for parameter structs; trailing
  `//!\brief`-documented members; `uz`/`u` literal suffixes. Match the surrounding file.
* **Includes** are grouped and separated by a blank line: standard library, then external (`cereal/`, `sharg/`,
  `seqan3/`, `hibf/`), then `raptor/`, then local `"..."` headers; alphabetical within a group.
* Reference code as `path/to/file.cpp:42` — it is clickable in the user's IDE.

## Index format and test fixtures

`raptor_index::version` is bumped **only at release**. Breaking the serialisation format between releases is
acceptable and expected — do not add version handling for it.

Adding a field to `CEREAL_SERIALIZE_FUNCTION_NAME` in `include/raptor/index.hpp` invalidates the checked-in
fixtures, and the failure is not always loud: an `.index` file errors with `Failed to read 8 bytes from input
stream`, but a `.hibf` file can silently misparse and still pass its test. Regenerate all thirteen:

* `test/data/{1,64,128}bins{19,23}window.index` — `raptor build --kmer 19 --window <W> --threads 1 --quiet` from a
  file listing `bin1.fa`..`bin4.fa` repeated 0/16/32 times (0 repetitions means `bin1.fa` alone).
* `test/data/hibf/{1,64,128}bins{19,23}window.hibf` — the same, with `--input <N>bins.layout`.
* `test/data/hibf/three_levels.hibf` — from `three_levels.layout`, `--kmer 19 --window 19`.

The fixtures bake in the absolute paths of the build directory that produced them, because `.layout` files are
`configure_file`d per build directory. Check which paths the current files contain (`strings <fixture> | head`) and
regenerate from a build directory that reproduces them, otherwise the diff carries a path change on top of the
format change. Every file should then grow or shrink by exactly the size of the changed field — verify that before
committing. `test/unit/api/index_size.cpp` computes from IBF bit sizes, not file size, so it usually needs no
update; check rather than assume.

## The dynamic index (`raptor update`)

`raptor update {insert,delete}` adds or removes user bins without rebuilding from all input data. Implementation in
`src/update/`, documented for users in `doc/usage/07_update/index.md`.

* **Only an HIBF built from a layout with `--empty-bin-fraction > 0` can be updated.** hibf only maintains
  `ibf.occupancy` when the layout reserved empty bins (`track_occupancy` is otherwise `false` and the vector stays
  all-zero), and insertion locates free technical bins through it. `src/argument_parsing/update_parsing.cpp`
  rejects unsuitable indices up front; `try_claim_bins` asserts the precondition.
* An update never modifies its input; it writes a new index to `--output`, which must not already exist.
* Deleting does not renumber user bins. Deleted entries keep their IDs and stay in the bin path list; their
  technical bins are reclaimed by later insertions.
* Emptying an IBF completely retires it and its subtree (`src/update/subtree.hpp`). IBFs are never erased from
  `ibf_vector` — that would invalidate every stored IBF index — but their data is released.
* `ibf_bin_to_user_bin_id` is the single source of truth for whether a bin is merged. `next_ibf_id` is only read
  for bins marked `merged`, but hibf documents `next_ibf_id[i][b] == i` for non-merged bins, so keep it consistent.
* Placement is heuristic and tuned: candidate weighting in `find_insert_location`, the `tmax_slack` allowance, the
  per-level FPR factor in `is_fpr_exceeded.hpp`, and the resize budget
  (`number_of_resized_ibfs() < 2 * original_number_of_ibfs()`). Do not adjust these without being asked, and
  benchmark when you do — index size is sensitive to them and the effects are counter-intuitive.

### Testing update changes

`test/unit/cli/update/` builds its own indices per test, because the shipped fixtures have no empty bins. The
oracle: searching `query.fq` with `--error 2` reports every user bin, with `--error 0` every bin except those
holding `bin4.fa`.

For anything algorithmic, the CLI tests are not enough — their four bins are too small to produce multiple IBFs.
Generate a synthetic set (~120 bins of varying size, each carrying a unique marker sequence used as its query),
build with a small `tmax`, insert in batches, and assert that every user bin is reported with no false positives.
That workload exercises partial rebuilds, full rebuilds and subtree retirement; the four-bin CLI tests do not.

## Known issues

* `raptor build` aborts on `hibf/misc/timer.hpp:277`, `assert(count.load() > 0u)` in
  `concurrent_timer::avg_in_seconds()`, when printing statistics for a timer that never ran (Debug builds,
  reproducible with 8 user bins at `--tmax 256`). `--quiet` avoids it. The index is written before the abort.
* `find_insert_location`'s final `return std::nullopt` assumes a layout sizes the top-level IBF such that growing
  it exceeds `tmax`. That holds for every layout `raptor layout` produces but is not guaranteed in general; the
  assert documents it.
