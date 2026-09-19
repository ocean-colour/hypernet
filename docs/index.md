# WATERHYPERNET Data Release 2 — an independent data user's guide

**J. Xavier Prochaska** (University of California, Santa Cruz) · 18 September 2026

*An independent guide to the WATERHYPERNET Release 2 archive, written from the
perspective of a data user preparing the dataset for ocean-colour algorithm
work. It is intended to complement — not replace — the official release notes.
A draft was reviewed by Kevin Ruddick, whose comments are incorporated
throughout; any remaining errors are mine.*

Every number here was measured from the archive itself or read from the release
notes. Scope of each measurement is stated in
[How these numbers were produced](provenance.md), and all of
the code is available (see the same section) so any figure can be reproduced or
extended.

**Where to start.** If you just want to read a spectrum, go to
[Quickstart](quickstart.md) and then
[Getting Rrs right](rrs.md) — between them they cover the two things
most likely to go silently wrong. If you are choosing which sites to work with,
start at [Choosing your data](choosing.md).

---

```{toctree}
:maxdepth: 2
:caption: Using the archive

about
quickstart
rrs
reading-many-files
traps
quality
choosing
over-subtraction
iop-retrieval
licence
provenance
```

<!-- Slot for the future package-API section. When `hypernet` grows its own
     documentation, add a second toctree here, e.g.

```{toctree}
:maxdepth: 2
:caption: The hypernet package

api/index
```

     Sphinx reads .rst and .md in one project, so those pages may be written in
     reStructuredText with autodoc + napoleon, as the sibling IOPtics docs are. -->
