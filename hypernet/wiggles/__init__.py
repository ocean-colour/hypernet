"""The wiggles effort: diagnosing and removing spectral wiggles in
WATERHYPERNET reflectance (plan: ``docs/wiggles_planning.md``; prompts:
``claude_prompts/wiggles/``).

Phase 0 scripts, run as modules from the repository root::

    python -m hypernet.wiggles.phase0a_veit          # task 5: VEIT line fits
    python -m hypernet.wiggles.phase0a_template_veit # task 3b: HSRS template SRF
    python -m hypernet.wiggles.phase0b_index         # task 9 ... 13: the delivery

Committed products (small CSV tables, figures) are written next to the code,
in this directory and ``figs/phase0/``; data intermediates (parquet) go to
:data:`OUT`, outside the repository.
"""

import os

#: This directory: committed CSV tables are written here.
WIGGLES_DIR = os.path.dirname(os.path.abspath(__file__))
#: Committed Phase 0 figures.
FIGDIR = os.path.join(WIGGLES_DIR, 'figs', 'phase0')
#: The repository root (for ``docs/`` and ``hypernet/data/``).
REPO = os.path.abspath(os.path.join(WIGGLES_DIR, '..', '..'))
#: The package data directory (SRF model JSON files).
DATA_DIR = os.path.join(REPO, 'hypernet', 'data')
#: Phase 0 data intermediates, outside the repository.
OUT = os.path.join(os.getenv('OS_COLOR', '.'), 'hypernet', 'wiggles', 'phase0')
