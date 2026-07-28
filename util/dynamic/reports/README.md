<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Interactive Reports

This directory contains the static HTML reports. You can either open the HTML files directly in your browser or serve
the directory locally (recommended) so that relative links and assets work correctly.

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

## Troubleshooting

- `Address already in use`: choose a different port, e.g., `8080`.
- Missing assets or broken links: ensure you are serving the `util/dynamic/reports` directory (not a parent/child directory that lacks files).
- No `python3` command: try `python` instead of `python3` if your system's Python 3 is available as `python`.
