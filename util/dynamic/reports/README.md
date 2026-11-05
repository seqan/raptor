<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Interactive Reports

This directory contains the static HTML reports, which are also published at https://docs.seqan.de/dynamic.

The reports are self-contained, so you can either open `index.html` directly in your browser or serve the directory
locally.

## Quick start

- From the repository root (recommended):

```bash
python3 -m http.server --directory util/dynamic/reports --bind 127.0.0.1 8000
```

- Or change into the reports directory and serve on a chosen port:

```bash
cd util/dynamic/reports
python3 -m http.server --bind 127.0.0.1 8000
```

Then open the reports at `http://127.0.0.1:8000` in your browser.

## Updating reports

After adding or regenerating a report, run `util/dynamic/reports/patch_reports.py`. It adds a licence header and a link
back to `index.html` to each report and sets the page title to the file name. A new report also needs an entry in
`index.html`.

## Troubleshooting

- `Address already in use`: choose a different port, e.g., `8080`.
- Broken links: ensure you are serving the `util/dynamic/reports` directory (not a parent/child directory that lacks
  files).
- No `python3` command: try `python` instead of `python3` if your system's Python 3 is available as `python`.
