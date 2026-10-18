"""Phase 2 task 3 check: Phase 1 still reproduces through the hardened edinterp.

Recomputes the Phase 1 task 9 metrics for one OSOAA scene with
``phase1_twin.run_scene`` (which now goes through ``hypernet.srf``'s cached
convolution matrix and the rewritten ``interpolate_ed_to_l``) and compares
them with the stored ``$OS_COLOR/hypernet/wiggles/phase1/twin_metrics.parquet``.
Writes nothing.

Run from the repository root:
``python -m hypernet.wiggles.phase2_check_phase1_repro [scene_index]``.
"""
import os
import sys

import numpy as np
import pandas as pd

from hypernet.wiggles import phase1_twin as p


def main(i_scene=0):
    rows, _ = p.run_scene((i_scene, i_scene))
    new = pd.DataFrame(rows)
    old = pd.read_parquet(os.path.join(p.ROOT, 'twin_metrics.parquet'))
    old = old[old['scene'] == new['scene'].iloc[0]]
    k = ['instrument', 'case', 'method', 'h']
    m = new.merge(old, on=k, suffixes=('_new', '_old'))
    rel = np.abs(m['rms_line_new'] / m['rms_line_old'] - 1)
    print('scene %s: %d rows; max relative difference in rms_line %.2e'
          % (new['scene'].iloc[0], len(m), rel.max()))
    print(m.loc[rel.idxmax(), k + ['rms_line_new', 'rms_line_old']].to_dict())
    return rel.max()


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
