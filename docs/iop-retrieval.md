# Using these data for IOP retrieval

A reflectance-only archive can exercise inherent-optical-property (IOP) retrieval
algorithms but cannot, on its own, validate what they retrieve. An algorithm can
be run on these spectra and the Rrs residual scored; algorithms can be compared
against one another for consistency and for the spread of their retrieved IOPs.
But nothing in the release constrains retrieved `a_ph`, `a_dg`, `bb_p`,
chlorophyll or the spectral slopes, so no truth-referenced error metric is
available. That is a statement about the intended scope of the network, not a
defect.

Three things nonetheless make the archive valuable for this work:

- **Optical diversity with a stable instrument.** Eleven sites spanning a factor
  of 12 in brightness and 0.52–1.03 in blue–green ratio, measured with two
  well-characterised systems on a fixed geometry, is a good stress test for
  algorithm robustness — particularly the very turbid sites, where many
  ocean-colour IOP algorithms degrade.
- **Hyperspectral resolution.** The HYPSTAR grid (~0.49 nm) is far finer than
  the satellite bands most IOP algorithms were built for, so the same spectrum
  can be convolved to PACE, MODIS, SeaWiFS or SBG bands and the retrieval
  compared across spectral configurations.
- **Repeat measurements at fixed stations.** Thousands of spectra from one
  location allow the *stability* of a retrieval to be assessed — how much of the
  retrieved variance is water and how much is algorithm noise — which
  single-visit field campaigns cannot support.

One such effort is [IOPtics](https://github.com/ocean-colour/IOPtics)
([documentation](https://ioptics.readthedocs.io/en/stable/)), a framework for
testing and comparing IOP retrieval algorithms, which is what motivated this
exploration.

Looking ahead, the measurands Kevin Ruddick would most like to see added to the
sites — automated flow cytometry (phytoplankton composition), turbidimeters and
chlorophyll fluorimeters — would be **more** directly useful for validating
in-water retrievals than IOPs would, and he notes that IOPs are usually an
intermediate, explanatory parameter rather than an end-user product. Chlorophyll
and turbidity/TSS also map onto truth quantities that retrieval-evaluation
frameworks already handle, so such data could be folded in without new machinery.
