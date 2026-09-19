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

## Tests

```bash
pytest
```

## Related

- [IOPtics](https://github.com/ocean-colour/IOPtics) — inherent optical
  property algorithm testing; a local checkout lives alongside this repository.
