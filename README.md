# hypernet
WATER HYPERNET related code and more

Analysis code for the [WATER HYPERNET](https://www.hypernets.eu/) network of
automated hyperspectral radiometers (HYPSTAR instruments), which provide
water-leaving radiance and reflectance for satellite ocean-colour validation.

## What it does

- Ingests and analyzes HYPSTAR water-leaving radiance and reflectance products.
- Develops improved wavelength calibration for the instruments.
- Produces new estimates of the uncertainty in the radiance and reflectance
  measurements.
- Generates metrics and diagnostic plots to share with the community.

## Installation

hypernet targets Python >= 3.12 and is developed against the `ocean14` conda
environment.

```bash
git clone https://github.com/ocean-colour/hypernet.git
cd hypernet
pip install -r requirements.txt
pip install -e .
```

To build the documentation locally:

```bash
pip install -r docs/requirements.txt
python -m sphinx -b html docs docs/_build/html
```

## The WATERHYPERNET Release-2 exploration

An indexing and sampling pass over the 56,669-file Release-2 archive, written up
as a user guide to the archive under [`docs/`](docs/) and published with Read the
Docs. Run from the repository root, with the archive reachable at
`$OS_COLOR/WATERHYPERNET/RELEASE_2`:

```bash
python -m hypernet.whn_explore 1       # index all files from their names
python -m hypernet.whn_explore 2       # read a pool, cluster, sample ~100/site
python docs/whn_figures.py             # summary table + figures into docs/
python -m hypernet.whn_simspec_check   # over-subtraction scan at the dark sites
```

Parquet/npz intermediates land outside the repository, under
`$OS_COLOR/hypernet/whn_explore`.

## Tests

```bash
pytest
```

Tests that need the archive skip automatically when it is not mounted.

## Related

- [IOPtics](https://github.com/ocean-colour/IOPtics) — inherent optical
  property algorithm testing; a local checkout lives alongside this repository.
