"""Render a Markdown document in this directory to a paginated PDF.

Written because none of the usual converters (pandoc, weasyprint, wkhtmltopdf,
python-markdown) are installed in ``ocean14``; the only PDF engine available is
the headless Chrome already on the machine. The Markdown subset handled here is
exactly the subset used by the documents in this directory -- headings,
paragraphs, bold/italic/inline code, fenced code blocks, pipe tables with
alignment, images, links, block quotes, horizontal rules and simple lists.

Usage::

    python make_pdf.py                       # WATERHYPERNET.md -> WATERHYPERNET.pdf
    python make_pdf.py respond_to_kevin.md
"""

from __future__ import annotations

import html
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CHROME = ('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')

CSS = """
@page { size: A4; margin: 18mm 16mm 18mm 16mm; }
body { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
       font-size: 10.5pt; line-height: 1.45; color: #1a1a1a; }
h1 { font-size: 19pt; margin: 0 0 2pt 0; line-height: 1.25; }
h2 { font-size: 14pt; margin: 20pt 0 6pt 0; padding-bottom: 3pt;
     border-bottom: 1px solid #ccc; page-break-after: avoid; }
h3 { font-size: 11.5pt; margin: 14pt 0 4pt 0; page-break-after: avoid; }
p { margin: 0 0 7pt 0; }
em.subtitle { color: #444; }
table { border-collapse: collapse; margin: 8pt 0 12pt 0; font-size: 8.6pt;
        width: 100%; page-break-inside: avoid; }
th, td { border: 1px solid #d0d0d0; padding: 3pt 5pt; }
th { background: #f2f2f2; text-align: left; font-weight: 600; }
code { font-family: "SF Mono", Menlo, Consolas, monospace; font-size: 9pt;
       background: #f4f4f4; padding: 1pt 3pt; border-radius: 2px; }
pre { background: #f6f6f6; border: 1px solid #e2e2e2; border-radius: 3px;
      padding: 7pt 9pt; overflow-x: auto; page-break-inside: avoid; }
pre code { background: none; padding: 0; font-size: 8.6pt; line-height: 1.35; }
blockquote { margin: 8pt 0; padding: 4pt 12pt; border-left: 3px solid #999;
             background: #fafafa; }
blockquote p { margin: 0; }
img { max-width: 100%; display: block; margin: 10pt auto; }
figure { margin: 12pt 0; page-break-inside: avoid; }
hr { border: none; border-top: 1px solid #ddd; margin: 16pt 0; }
ul, ol { margin: 0 0 8pt 0; padding-left: 20pt; }
li { margin-bottom: 3pt; }
a { color: #0b5394; text-decoration: none; }
.meta { color: #555; font-size: 10pt; margin-bottom: 10pt; }
"""

_INLINE = (
    (re.compile(r'`([^`]+)`'), lambda m: f'<code>{html.escape(m.group(1))}</code>'),
    (re.compile(r'!\[([^\]]*)\]\(([^)]+)\)'),
     lambda m: f'<img src="{m.group(2)}" alt="{html.escape(m.group(1))}">'),
    (re.compile(r'\[([^\]]+)\]\(([^)]+)\)'),
     lambda m: f'<a href="{m.group(2)}">{m.group(1)}</a>'),
    (re.compile(r'\*\*([^*]+)\*\*'), lambda m: f'<strong>{m.group(1)}</strong>'),
    (re.compile(r'(?<![\w*])\*([^*\n]+)\*(?![\w*])'),
     lambda m: f'<em>{m.group(1)}</em>'),
)


def inline(text):
    """Apply inline Markdown (code, images, links, bold, italic) to one line.

    Code spans are converted first and their contents escaped, so that a
    variable name containing markup characters survives intact.
    """
    out, pos, parts = '', 0, []
    # Protect code spans before escaping the rest of the line.
    for m in re.finditer(r'`[^`]+`', text):
        parts.append(('text', text[pos:m.start()]))
        parts.append(('code', m.group(0)[1:-1]))
        pos = m.end()
    parts.append(('text', text[pos:]))

    for kind, chunk in parts:
        if kind == 'code':
            out += f'<code>{html.escape(chunk)}</code>'
            continue
        chunk = html.escape(chunk)
        for pattern, repl in _INLINE[1:]:
            chunk = pattern.sub(repl, chunk)
        out += chunk
    return out


def _table(rows):
    """Render a collected pipe-table block, honouring the alignment row."""
    cells = [[c.strip() for c in r.strip().strip('|').split('|')] for r in rows]
    header, align_row, body = cells[0], cells[1], cells[2:]

    align = []
    for a in align_row:
        if a.startswith(':') and a.endswith(':'):
            align.append('center')
        elif a.endswith(':'):
            align.append('right')
        else:
            align.append('left')

    out = ['<table>', '<thead><tr>']
    for i, h in enumerate(header):
        out.append(f'<th style="text-align:{align[i] if i < len(align) else "left"}">'
                   f'{inline(h)}</th>')
    out.append('</tr></thead><tbody>')
    for row in body:
        out.append('<tr>')
        for i, c in enumerate(row):
            out.append(f'<td style="text-align:{align[i] if i < len(align) else "left"}">'
                       f'{inline(c)}</td>')
        out.append('</tr>')
    out.append('</tbody></table>')
    return '\n'.join(out)


def to_html(md, title):
    """Convert the supported Markdown subset to a standalone HTML document."""
    lines = md.split('\n')
    body, i = [], 0
    list_stack = []                       # 'ul' / 'ol' currently open

    def close_lists():
        while list_stack:
            body.append(f'</{list_stack.pop()}>')

    while i < len(lines):
        line = lines[i]

        if line.startswith('```'):                       # fenced code
            close_lists()
            i += 1
            buf = []
            while i < len(lines) and not lines[i].startswith('```'):
                buf.append(html.escape(lines[i]))
                i += 1
            body.append('<pre><code>' + '\n'.join(buf) + '</code></pre>')
            i += 1
            continue

        if re.match(r'^\|.*\|\s*$', line) and i + 1 < len(lines) \
                and re.match(r'^\|[\s:|-]+\|\s*$', lines[i + 1]):
            close_lists()
            rows = []
            while i < len(lines) and re.match(r'^\|.*\|\s*$', lines[i]):
                rows.append(lines[i])
                i += 1
            body.append(_table(rows))
            continue

        if re.match(r'^\s*$', line):
            close_lists()
            i += 1
            continue

        if re.match(r'^---+\s*$', line):
            close_lists()
            body.append('<hr>')
            i += 1
            continue

        m = re.match(r'^(#{1,4})\s+(.*)$', line)
        if m:
            close_lists()
            lvl = len(m.group(1))
            body.append(f'<h{lvl}>{inline(m.group(2))}</h{lvl}>')
            i += 1
            continue

        if line.startswith('>'):
            close_lists()
            buf = []
            while i < len(lines) and lines[i].startswith('>'):
                buf.append(lines[i].lstrip('>').strip())
                i += 1
            body.append(f'<blockquote><p>{inline(" ".join(buf))}</p></blockquote>')
            continue

        m = re.match(r'^(\s*)([-*]|\d+\.)\s+(.*)$', line)
        if m:
            kind = 'ol' if m.group(2)[0].isdigit() else 'ul'
            if not list_stack or list_stack[-1] != kind:
                close_lists()
                body.append(f'<{kind}>')
                list_stack.append(kind)
            text = m.group(3)
            i += 1
            # absorb continuation lines of the same item
            while i < len(lines) and re.match(r'^\s{2,}\S', lines[i]) \
                    and not re.match(r'^\s*([-*]|\d+\.)\s', lines[i]):
                text += ' ' + lines[i].strip()
                i += 1
            body.append(f'<li>{inline(text)}</li>')
            continue

        # paragraph: gather until a blank line or a block construct
        buf = [line]
        i += 1
        while i < len(lines) and lines[i].strip() and \
                not re.match(r'^(#{1,4}\s|```|\||>|---+\s*$|\s*([-*]|\d+\.)\s)',
                             lines[i]):
            buf.append(lines[i])
            i += 1
        close_lists()
        text = ' '.join(b.strip() for b in buf)
        cls = ' class="subtitle"' if text.startswith('*') and text.endswith('*') else ''
        body.append(f'<p{cls}>{inline(text)}</p>')

    close_lists()
    return (f'<!doctype html><html><head><meta charset="utf-8">'
            f'<title>{html.escape(title)}</title><style>{CSS}</style></head>'
            f'<body>{"".join(body)}</body></html>')


def render(md_path, pdf_path=None):
    """Convert ``md_path`` to PDF via headless Chrome. Returns the PDF path."""
    md_path = os.path.abspath(md_path)
    stem = os.path.splitext(md_path)[0]
    pdf_path = pdf_path or stem + '.pdf'
    html_path = stem + '.html'

    with open(md_path) as fh:
        md = fh.read()
    title = os.path.basename(stem)
    for line in md.split('\n'):
        if line.startswith('# '):
            title = line[2:].strip()
            break

    with open(html_path, 'w') as fh:
        fh.write(to_html(md, title))

    if not os.path.exists(CHROME):
        raise FileNotFoundError(f'Chrome not found at {CHROME}; HTML written to '
                                f'{html_path} -- print it manually.')

    subprocess.run(
        [CHROME, '--headless', '--disable-gpu', '--no-pdf-header-footer',
         f'--print-to-pdf={pdf_path}', '--virtual-time-budget=20000',
         f'file://{html_path}'],
        check=True, capture_output=True, timeout=180)

    os.remove(html_path)                 # the HTML is a build artifact
    return pdf_path


if __name__ == '__main__':
    target = sys.argv[1] if len(sys.argv) > 1 else 'WATERHYPERNET.md'
    if not os.path.isabs(target):
        target = os.path.join(HERE, target)
    out = render(target)
    print(f'wrote {out} ({os.path.getsize(out) / 1e6:.2f} MB)')
