#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

"""Adds a licence header, a link to index.html and a link that saves the plot as SVG to every report next to this script
and sets its page title to the file name.

Run this after adding or regenerating a report. Reports that are already patched stay byte-identical.
"""

import re
import sys
from pathlib import Path

# The ignore block keeps reuse from reading this header as the licence of this script.
# REUSE-IgnoreStart
LICENSE_HEADER = """<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-FileCopyrightText: Anaconda, Inc., and Bokeh Contributors
SPDX-License-Identifier: BSD-3-Clause
-->

"""
# REUSE-IgnoreEnd

LINK_STYLE = (
    "font: 14px -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; color: #0066cc; text-decoration: none"
)

SAVE_SVG = """
    document.getElementById('save-svg').addEventListener('click', event => {
      event.preventDefault();
      const plot = Bokeh.documents[0].roots()[0].children[0];
      const svg = Bokeh.index.find_one(plot).export('svg', false).ctx.get_svg();
      // Bokeh also strokes each text with an invisible outline, which turns into Type 3 fonts in a PDF.
      svg.querySelectorAll('text[fill="none"][stroke-opacity="0"]').forEach(text => text.remove());
      // Bokeh places the rasterised legend markers at their device pixel size.
      for (const image of svg.querySelectorAll('image')) {
        for (const side of ['width', 'height']) image.setAttribute(side, image.getAttribute(side) / devicePixelRatio);
      }
      const link = document.createElement('a');
      link.href = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)], {type: 'image/svg+xml'}));
      link.download = `${document.title}.svg`;
      link.click();
    });
  """

NAV = f"""<div id="report-nav" style="padding: 6px 8px">
  <a href="index.html" style="{LINK_STYLE}">&larr; Index</a>
  <a href="#" id="save-svg" style="{LINK_STYLE}; margin-left: 16px">Save SVG</a>
  <script>{SAVE_SVG}</script>
</div>
"""

OLD_LICENSE_HEADER = re.compile(r'\A<!--\nSPDX-.*?-->\n\n', re.DOTALL)
# Line-anchored, because the inlined BokehJS bundle contains "<body" as well.
BODY_TAG = re.compile(r'^[ \t]*<body>[ \t]*\n', re.MULTILINE)
OLD_NAV = re.compile(r'^<div id="report-nav".*?</div>\n', re.MULTILINE | re.DOTALL)
# The first match is the one in <head>, which precedes the bundle.
TITLE = re.compile(r'<title>.*?</title>')


def patch(report: Path) -> bool:
    """Patches a report in place and returns whether its content changed."""
    original = report.read_bytes().decode()

    html, titles = TITLE.subn(lambda _: f'<title>{report.stem}</title>', original, count=1)
    if titles != 1:
        sys.exit(f'{report.name}: no <title> tag found')

    html, bodies = BODY_TAG.subn(lambda match: match.group() + NAV, OLD_NAV.sub('', html))
    if bodies != 1:
        sys.exit(f'{report.name}: expected one <body> line, found {bodies}')

    html = LICENSE_HEADER + OLD_LICENSE_HEADER.sub('', html)

    if html == original:
        return False

    report.write_bytes(html.encode())
    return True


def main():
    for report in sorted(Path(__file__).resolve().parent.glob('*.html')):
        if report.name != 'index.html' and patch(report):
            print(f'Patched {report.name}')


if __name__ == '__main__':
    main()
