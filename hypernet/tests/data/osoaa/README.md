OSOAA output fixtures for `hypernet/tests/test_osoaa.py`.

One real run of `hypernet.rt.osoaa` (2026-10-03, JXP's OSOAA fork built with
gfortran 16.2.0): 550 nm, SZA 40°, relative azimuth 90°, AOT(550) = 0.1
(mono-modal log-normal, `AER.Waref` = 0.55 µm), Chl 1 mg m⁻³, wind 5 m/s,
100 m sea.  `ListParam.txt` lists every parameter.  The two Advanced files
are trimmed to levels 0 (TOA), 25, 26 (0+), 27 (0−) and 28 to keep the repo
small; their headers are unchanged.
