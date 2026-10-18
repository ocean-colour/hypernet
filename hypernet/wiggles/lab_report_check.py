"""Check the interim report (lab data prompts, task 7).

``reports/wiggles_interim_report.md``:

1. **Links.**  Every Markdown image or link with a relative target must
   resolve on disk, relative to the report's directory.  Figure numbers must
   run 1, 2, ... without gaps or repeats.
2. **Numbers.**  Table 1 (field ΔFWHM per instrument, E shape p) and Table 3
   (lab ΔFWHM) against the committed tables: ``hypernet/data/srf_table.csv``
   (``delta_fwhm`` of the Release 2 rows, ``shape_p``) and
   ``hypernet/wiggles/lab_fwhm_compare_summary.csv`` (``lab_dfwhm``), to the
   report's two decimals.

Exits non-zero if anything fails.

Run from the repository root:  ``python -m hypernet.wiggles.lab_report_check``
"""
import os
import re
import sys

import numpy as np
import pandas as pd

from hypernet.wiggles import DATA_DIR, REPO, WIGGLES_DIR

REPORT = os.path.join(REPO, 'reports', 'wiggles_interim_report.md')


def check_links(text, base):
    bad, n = [], 0
    for m in re.finditer(r'!?\[[^\]]*\]\(([^)\s]+)\)', text):
        target = m.group(1)
        if re.match(r'[a-z]+://', target) or target.startswith('#'):
            continue
        n += 1
        if not os.path.exists(os.path.normpath(os.path.join(base, target.split('#')[0]))):
            bad.append(target)
    figs = [int(x) for x in re.findall(r'!\[Figure (\d+)\]', text)]
    seq_ok = figs == list(range(1, len(figs) + 1))
    return n, bad, figs, seq_ok


def _table(text, title):
    """Rows of the Markdown table that follows ``title``."""
    i = text.index(title)
    rows = []
    for line in text[i:].splitlines()[1:]:
        if line.startswith('|'):
            cells = [c.strip() for c in line.strip('|').split('|')]
            if not set(cells[0]) <= set('-: '):
                rows.append(cells)
        elif rows:
            break
    return rows[1:]                                  # drop the header


def _nums(cell):
    """Numbers in a cell, outside brackets (the brackets hold old values)."""
    cell = re.sub(r'\([^)]*\)', '', cell).replace('−', '-')
    cell = re.sub(r'(\d)-(\d)', r'\1 \2', cell)       # a range a-b, not a minus sign
    return [float(x) for x in re.findall(r'[-+]?\d+\.\d+', cell)]


def check_numbers(text):
    t = pd.read_csv(os.path.join(DATA_DIR, 'srf_table.csv'))
    e = t[(t['channel'] == 'E')]
    lab = pd.read_csv(os.path.join(WIGGLES_DIR, 'lab_fwhm_compare_summary.csv'))
    out = []
    for cells in _table(text, '**Table 1.'):
        inst = 'HYPSTAR_' + cells[0]
        r2 = e[(e['instrument'] == inst) & (e['source'] == 'release2_l2b')]
        want = sorted(set(np.round(r2['delta_fwhm'], 2)))
        got = sorted(_nums(cells[2]))
        ok = (len(got) == 1 and min(abs(got[0] - w) for w in want) <= 0.005) or \
            (len(got) == 2 and abs(got[0] - min(want)) <= 0.005 and abs(got[1] - max(want)) <= 0.005)
        out.append(('Table 1 ΔFWHM field', cells[0], got, want, ok))
        p = sorted(set(np.round(e[e['instrument'] == inst]['shape_p'].dropna(), 2)))
        pg = _nums(cells[4])
        okp = (not p and not pg) or (bool(p) and all(min(abs(x - w) for w in p) <= 0.005 for x in pg))
        out.append(('Table 1 shape p', cells[0], pg, p, okp))
    for cells in _table(text, '**Table 3.'):
        inst = 'HYPSTAR_' + cells[0]
        want = sorted(set(np.round(lab[(lab['instrument'] == inst) &
                                       (lab['channel'] == 'E')]['lab_dfwhm'], 2)))
        got = _nums(cells[1])
        ok = all(min(abs(x - w) for w in want) <= 0.005 for x in got) and bool(got)
        out.append(('Table 3 ΔFWHM lab', cells[0], got, want, ok))
    return out


def main():
    text = open(REPORT).read()
    n, bad, figs, seq_ok = check_links(text, os.path.dirname(REPORT))
    print('links: %d relative targets, %d missing %s' % (n, len(bad), bad or ''))
    print('figures: %s, numbered in sequence: %s' % (figs, seq_ok))
    rows = check_numbers(text)
    nbad = 0
    for what, inst, got, want, ok in rows:
        if not ok:
            nbad += 1
        print('%-20s %-7s report %-16s table %-20s %s' % (what, inst, got, want, 'ok' if ok else 'MISMATCH'))
    fail = bool(bad) or not seq_ok or nbad
    print('FAIL' if fail else 'all checks passed')
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
