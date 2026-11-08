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
      const plot_view = Bokeh.index.find_one(plot);
      const layer = plot_view.export('svg', false);
      const svg = layer.ctx.get_svg();
      // Bokeh also strokes each text with an invisible outline, which turns into Type 3 fonts in a PDF.
      svg.querySelectorAll('text[fill="none"][stroke-opacity="0"]').forEach(text => text.remove());
      // The legend markers are bitmaps in legend order. Draw them again as vectors, the same way Bokeh draws them.
      const legend = plot.above.find(model => model.constructor.__name__ === 'Legend');
      const legend_view = plot_view.views.find_one(legend);
      const {glyph_width: width, glyph_height: height} = legend;
      svg.querySelectorAll('image').forEach((image, index) => {
        const {item, label} = legend_view.entries[index];
        const field = item.get_field_from_label_prop();
        const marker = new layer.constructor('svg', false);
        marker.resize(width, height);
        const ctx = marker.prepare();
        for (const renderer of item.renderers) {
          plot_view.views.find_one(renderer)?.draw_legend(ctx, 0, width, 0, height, field, label, item.index);
        }
        marker.finish();
        const group = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        group.setAttribute('transform', image.getAttribute('transform'));
        group.append(...ctx.get_svg().childNodes);
        image.replaceWith(group);
      });
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
