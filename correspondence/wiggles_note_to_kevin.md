# Note to Kevin Ruddick: spectral wiggles, first look and plan

**From:** J. Xavier Prochaska (UC Santa Cruz) · **Date:** 28 September 2026 ·
**Status:** sent by JXP (confirmed 2026-09-29)

**Attach:** `docs/wiggles_planning.md` (or its rendered page) and
`docs/wiggles_data_request.csv`

---

Kevin,

Thank you for the wiggles idea and the VEIT sample. We have started on it.
I've put together a plan, and before we go further I'd like to check one
early result with you.

**The symptom is real, but the cause may not be the one we expected.**  On
the one sequence you sent (VEIT, 2026-06-04, HYPSTAR 122304):

- Your processor's L1C irradiance is exactly the linear interpolation of L1A
  E onto the L grid.  The ρw second derivative does track the Fraunhofer
  structure in Ed: corr(ρw'', Ed''/Ed) = −0.56 over 400–700 nm.
- The linear-interpolation error itself is small, though.  The spectra are
  well sampled (~6 pixels per FWHM).  The E–L centre offsets we measure at
  clean solar lines are only 0.03–0.08 nm, and the interpolation accounts for
  only ~3 % of the line structure that survives in Ld/Ed.
- The larger effect appears to be that **E and L have different spectral
  response functions**.  The same solar lines are 0.5–0.9 nm narrower in E
  than in L in the blue, and the two converge by ~690 nm.  Lu and Ld agree
  with each other.  No interpolation on the E grid can remove a resolution
  mismatch.

**The good news is that your eq. (14) needs only a one-symbol change to handle
both cases.**  If the model irradiance in the numerator is convolved with the
*L* SRF instead of the E SRF, the correction also matches resolution.  With
equal SRFs it reduces exactly to your eq. (14), and with a flat model to
linear interpolation.  So we can test your hypothesis and ours side by side on
one code path, starting with the twin experiment you proposed.  The attached
plan sets this out in four phases, each with an explicit go/no-go criterion.
It ends with code for `hypernet/`, a drop-in snippet for the processor, and a
technical note with FIDUCEO-style uncertainties.  All of this is from one
sequence, so treat it as a lead.

**What we would need from you**, when you have time:

1. **L1A_IRR, L1A_RAD, L1C_ALL and L2A_REF** for the sequences in the
   attached CSV.
   - That is 150 primary plus 74 spares, over VEIT, BEFR/THFR, MAFR and
     (optionally) GAIT.
   - The picks were chosen to separate instrument from site effects, and to
     compare before and after recalibration.
   - Could you fill in the `sky` column (clear / overcast / broken), and swap
     in spares where a primary is unsuitable?
2. **The RAD and IRR calibration files** for HYPSTAR 122302, 122303, 122304,
   122305, 121231, 121222 and 120242.  For 122304 and 121222 we would need
   both versions, before and after recalibration.
3. **Any laboratory line-spread or slit-function data** behind those files.
   The `bandwidth` variable is a flat 3.0 nm, so I assume no per-pixel SRF is
   stored.
4. **The processor version** that produced Release 2 L1C.
5. **Air or vacuum?**  Is the HYPSTAR wavelength calibration on an air or a
   vacuum scale?  The two differ by 0.1–0.2 nm in the visible, which matters
   when we compare line centres with laboratory wavelengths and with TSIS-1
   HSRS (vacuum).
6. **Your view on the hypothesis.**  Is there a known reason the irradiance
   and radiance paths would give different line widths, e.g. diffuser vs
   fore-optics filling the slit differently?

I can confirm the calibration-attribute bug you mentioned: `_rad` points to
the IRR file in L1A_RAD, L1C and L2A.  We will pair calibration files by
instrument and date instead.

No rush on any of this.  The twin experiment and the SRF fitting on your
sample need nothing further from you, so we will get on with those in the
meantime.

Best,
Xavier
