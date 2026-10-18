"""Print the median normalised ringing curves ring(delta)/min of the lab
inversion check (``lab_inversion_check.py``) per instrument x period x grid,
to judge how sharp each minimum is.

Run from the repository root:  ``python -m hypernet.wiggles.lab_inversion_curves``
"""
import os

import pandas as pd

from hypernet import labcal

OUT = os.path.join(labcal.lab_root(), 'inversion')


def main():
    c = pd.read_csv(os.path.join(OUT, 'ring_scans.csv'))
    c['n'] = c['ring'] / c.groupby(['instrument', 'cal_period', 'site', 'sequence_time',
                                    'grid'])['ring'].transform('min')
    m = c.groupby(['instrument', 'cal_period', 'grid', 'delta'])['n'].median().unstack('delta')
    pd.set_option('display.width', 300)
    print(m[[d for d in m.columns if abs(round(d * 100) % 2) == 0]].round(3).to_string())


if __name__ == '__main__':
    main()
