// Build docs/slides/wiggles_phase0_report.pptx (Phase 0 task 8b).
// Numbers come from wiggles_phase0_report_data.json, written by
// wiggles_phase0_report_data.py from the committed Phase 0a products.
// Rule for this deck: no text smaller than 20 pt.
//
//   python docs/slides/wiggles_phase0_report_data.py
//   node docs/slides/build_wiggles_phase0_report.js
const path = require('path');
const pptxgen = require('pptxgenjs');
const D = require(path.join(__dirname, 'wiggles_phase0_report_data.json'));

const C = {
  navy: '0B2545', teal: '13505B', tint: 'E6F0F1', coral: 'D9542C', ink: '1F2A30',
  muted: '4A5A60', white: 'FFFFFF', line: 'C9D8DA',
};
const HEAD = 'Cambria', BODY = 'Calibri';
const W = 13.333, H = 7.5, M = 0.5;

const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';
pres.title = 'WATERHYPERNET wiggles: Phase 0a report';
pres.author = 'J. Xavier Prochaska';

function title(s, text) {
  s.addText(text, { x: M, y: 0.35, w: W - 2 * M, h: 0.85, fontFace: HEAD, fontSize: 36,
    bold: true, color: C.navy, margin: 0, valign: 'middle', isTextBox: true });
}

function badge(s, x, y, d, label, fill) {
  s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: fill || C.teal },
    line: { color: fill || C.teal } });
  s.addText(label, { x, y, w: d, h: d, fontFace: HEAD, fontSize: 24, bold: true,
    color: C.white, align: 'center', valign: 'middle', margin: 0, isTextBox: true });
}

function bullets(s, items, box, size) {
  s.addText(items.map((t, i) => ({ text: t, options: { bullet: true,
    breakLine: i < items.length - 1 } })),
  Object.assign({ fontFace: BODY, fontSize: size || 22, color: C.ink, valign: 'top',
    paraSpaceAfter: 14, margin: 0.05, isTextBox: true }, box));
}

function image(s, key, box) {
  const im = D.images[key];
  let w = box.w, h = w / im.aspect;
  if (h > box.h) { h = box.h; w = h * im.aspect; }
  s.addImage({ path: im.path, x: box.x, y: box.y, w, h });
  return { w, h };
}

const ROW_H = 0.5;
function table(s, header, rows, box, colW) {
  const hdr = header.map(t => ({ text: t, options: { bold: true, color: C.white,
    fill: { color: C.teal }, align: 'center', valign: 'middle' } }));
  const body = rows.map((r, i) => r.map((t, j) => ({ text: t, options: {
    fill: { color: i % 2 ? C.white : C.tint }, align: j === 0 ? 'left' : 'center',
    valign: 'middle', bold: j === 0 } })));
  s.addTable([hdr].concat(body), Object.assign({ fontFace: BODY, fontSize: 20,
    color: C.ink, border: { type: 'solid', pt: 0.75, color: C.line }, colW,
    rowH: ROW_H, margin: [0.04, 0.1, 0.04, 0.1] }, box));
  return box.y + (rows.length + 1) * ROW_H;   // bottom of the table
}

function note(s, text, y) {
  s.addText(text, { x: M, y, w: W - 2 * M, h: 0.5, fontFace: BODY, fontSize: 20,
    italic: true, color: C.muted, margin: 0, isTextBox: true });
}

// 1. Title -------------------------------------------------------------------
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addText('Spectral wiggles in WATERHYPERNET', { x: M + 0.3, y: 1.7, w: 12, h: 1.0,
    fontFace: HEAD, fontSize: 44, bold: true, color: C.white, margin: 0, isTextBox: true });
  s.addText('Phase 0a report: what one VEIT sequence tells us about the E and L spectral '
    + 'response functions', { x: M + 0.3, y: 2.8, w: 11.5, h: 1.2, fontFace: BODY,
    fontSize: 26, color: 'CFE3E6', margin: 0, isTextBox: true });
  s.addText('HYPSTAR 122304 · VEIT · 2026-06-04 08:45 (SEQ20260604T084543)',
    { x: M + 0.3, y: 4.4, w: 12, h: 0.5, fontFace: BODY, fontSize: 22, color: C.white,
      margin: 0, isTextBox: true });
  s.addText('J. Xavier Prochaska · UC Santa Cruz · 29 September 2026',
    { x: M + 0.3, y: 5.9, w: 12, h: 0.5, fontFace: BODY, fontSize: 20, color: 'CFE3E6',
      margin: 0, isTextBox: true });
  s.addNotes('Phase 0a covers the single VEIT sequence Kevin sent. Everything here is one '
    + 'sequence from one instrument: a lead, not a result. Phase 0b repeats it on the ~150 '
    + 'requested sequences.');
}

// 2. Two hypotheses ------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Two explanations for the wiggles in ρw');
  const cards = [
    ['H1', 'Interpolation (Kevin)', ['E and L sit on different grids', 'Linear interpolation of '
      + 'Ed onto the L grid mis-represents the Fraunhofer lines'], C.teal],
    ['H2', 'SRF mismatch (ours)', ['E and L take different optical paths', 'Different line '
      + 'widths leave the lines in L/E, whatever the interpolation'], C.coral],
  ];
  cards.forEach(([b, head, body, col], i) => {
    const x = M + i * 6.3, y = 1.6, w = 5.9, h = 4.2;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12,
      fill: { color: C.tint }, line: { color: C.tint } });
    badge(s, x + 0.35, y + 0.35, 0.9, b, col);
    s.addText(head, { x: x + 1.45, y: y + 0.35, w: w - 1.7, h: 0.9, fontFace: HEAD,
      fontSize: 26, bold: true, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    bullets(s, body, { x: x + 0.35, y: y + 1.6, w: w - 0.7, h: h - 1.8 });
  });
  s.addText('Both are tested on one code path: a generalised Ruddick et al. (2023) eq. 14',
    { x: M, y: 6.2, w: W - 2 * M, h: 0.6, fontFace: BODY, fontSize: 22, italic: true,
      color: C.muted, margin: 0, isTextBox: true });
  s.addNotes('Kevin\'s eq. 14 uses the E SRF on both sides, so it fixes grid offsets but not '
    + 'a width mismatch. Putting the L-SRF-convolved model in the numerator handles both; '
    + 'equal SRFs recover eq. 14, and a flat model recovers linear interpolation.');
}

// 3. What was built --------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'What Phase 0a built');
  const steps = [
    ['1', 'Read', 'L1A/L1C/L2A readers, checked against the processor (whn_l1a)'],
    ['2', 'Fit lines', 'Gaussian fits of 13 solar lines in E, Ld, Lu; FWHM(λ) models (srf)'],
    ['3', 'Template SRF', 'TSIS-1 HSRS × Gaussian SRF in 10 nm windows (refspec, srf)'],
    ['4', 'Budget', 'Residual line structure vs the H1 and H2 predictions'],
  ];
  steps.forEach(([n, head, body], i) => {
    const y = 1.5 + i * 1.18;
    badge(s, M, y, 0.85, n);
    s.addText(head, { x: M + 1.15, y, w: 2.6, h: 0.85, fontFace: HEAD, fontSize: 24,
      bold: true, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(body, { x: M + 3.8, y, w: 8.5, h: 0.85, fontFace: BODY, fontSize: 22,
      color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
  });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: M, y: 6.3, w: W - 2 * M, h: 0.8,
    rectRadius: 0.1, fill: { color: C.tint }, line: { color: C.tint } });
  s.addText('Phase 0b pipeline (tasks 9–13) written and tested on VEIT; waiting on the '
    + 'L1A delivery', { x: M + 0.3, y: 6.3, w: W - 2 * M - 0.6, h: 0.8, fontFace: BODY,
    fontSize: 22, color: C.teal, bold: true, valign: 'middle', margin: 0, isTextBox: true });
  s.addNotes('All code is in the hypernet package: srf.py, whn_l1a.py, whn_srf.py, '
    + 'refspec.py, and the Phase 0 scripts in hypernet/wiggles/. 55 tests pass.');
}

// 4. What the files taught us -------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Four things the VEIT files taught us');
  const cards = [
    ['vza < 90', 'is the water view (Lu): vza is measured from nadir'],
    ['interpolated', 'L1C Ld: the two sky series, interpolated to the water-view time'],
    ['0', 'every L1C u_rel_random value: a placeholder, not an uncertainty'],
    ['5 % vs 0.4 %', 'broadband scan-to-scan change vs per-pixel noise'],
  ];
  cards.forEach(([big, small], i) => {
    const x = M + (i % 2) * 6.3, y = 1.5 + Math.floor(i / 2) * 2.75, w = 5.9, h = 2.45;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12,
      fill: { color: C.tint }, line: { color: C.tint } });
    s.addText(big, { x: x + 0.35, y: y + 0.25, w: w - 0.7, h: 1.0, fontFace: HEAD,
      fontSize: 40, bold: true, color: C.coral, margin: 0, valign: 'middle', isTextBox: true });
    s.addText(small, { x: x + 0.35, y: y + 1.3, w: w - 0.7, h: 1.0, fontFace: BODY,
      fontSize: 22, color: C.ink, margin: 0, valign: 'top', isTextBox: true });
  });
  s.addNotes('The first-look scripts had Ld and Lu swapped; the plan table headers are '
    + 'corrected. Because the scan scatter is broadband, raw std/sqrt(N) overstates the '
    + 'per-pixel error 7-14x; we use scatter after removing a smooth scale per scan.');
}

// 5. Task 5 figure ------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'L lines are wider than E lines (task 5)');
  const im = image(s, 'veit_fwhm_vs_lambda', { x: M, y: 1.35, w: 7.4, h: 5.8 });
  bullets(s, [
    'Ld − E: +0.4 to +0.6 nm from 430 to 590 nm, each line at 5–15σ',
    'The gap closes to the red: ~0.1 nm at Ca II 854–866',
    'Ld and Lu agree: one L SRF',
  ], { x: M + im.w + 0.4, y: 1.5, w: W - 2 * M - im.w - 0.4, h: 5.4 });
  s.addNotes('Gaussian fits of each solar line in the mean E, Ld and Lu spectra, weighted by '
    + 'the flattened scan scatter. Filled symbols are clean lines, open ones blends. The '
    + 'quadratic fits have chi2_nu 5-26 because each line has its own intrinsic width.');
}

// 6. Task 5 table ---------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Gaussian line widths, FWHM in nm (task 5)');
  const y5 = table(s, ['Line', 'E', 'Ld', 'Lu', 'Ld − E', 'Signif.'], D.task5,
    { x: M, y: 1.4, w: W - 2 * M }, [3.0, 1.6, 1.6, 1.6, 2.6, 1.93]);
  note(s, 'These widths include each line\'s own profile; the template fit removes it',
    y5 + 0.3);
  s.addNotes('Mean spectra of 6 scans per channel. Errors from curve_fit with absolute sigma '
    + '= flattened scan-to-scan error of the mean. Ld is the sky view, Lu the water view.');
}

// 7. Task 6 figure ---------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Per-scan fits check the errors (task 6)');
  const im = image(s, 'veit_scan_scatter', { x: M, y: 1.3, w: 7.0, h: 5.95 });
  bullets(s, [
    'Scan spread ≈ fit error: the error model is right',
    'Widths steady between the two sky series, 50 s apart',
    'Ld = Lu to about 0.1 nm at the clean lines',
    'L lines 10–17 % shallower, but same equivalent width: resolution, not Ring',
  ], { x: M + im.w + 0.4, y: 1.5, w: W - 2 * M - im.w - 0.4, h: 5.6 });
  s.addNotes('Every line fitted in each of the 6 scans per channel. Panels: FWHM and centroid '
    + 'scan std over median fit error; Ld minus Lu FWHM; L/E depth (open) and equivalent '
    + 'width (filled). Lu is 0.16-0.33 nm wider at Ca H and the G band, plausibly water-leaving '
    + 'content; Ld is the reference L channel.');
}

// 8. Task 6 table ------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Scan scatter, Ld vs Lu, and the Ring check (task 6)');
  const y6 = table(s, ['Line', 'Ld scan std / fit err', 'Ld − Lu (nm)', 'Depth Ld/E',
    'EW Ld/E'], D.task6, { x: M, y: 1.4, w: W - 2 * M }, [2.8, 2.9, 2.6, 2.02, 2.01]);
  note(s, 'EW Ld/E ≈ 1: a wider SRF lowers depth but keeps equivalent width', y6 + 0.3);
  s.addNotes('Ring filling-in would lower both depth and equivalent width. EW ratios of '
    + '0.94-1.09 leave room for at most a few per cent of filling. In Lu the red EW ratios fall '
    + 'to 0.64-0.85, suggesting an additive line-free component there (stray light or dark).');
}

// 9. Task 7 figure -------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'H2 accounts for the line structure (task 7)');
  const im = image(s, 'veit_budget', { x: M, y: 1.3, w: 6.5, h: 5.95 });
  bullets(s, [
    'Measured Ld/Ed residual matches the H2 prediction in size and shape',
    'H1 is ~' + D.task7_stats.h2_over_h1 + '× smaller and has the opposite sign',
    'ρw follows H2 from Ca H/K to Na D',
  ], { x: M + im.w + 0.4, y: 1.5, w: W - 2 * M - im.w - 0.4, h: 5.6 });
  s.addNotes('Top: residual profiles in Ld/Ed with the H1 (linear minus cubic interpolation) '
    + 'and H2 (E line re-observed at the L width, equivalent width conserved) predictions. '
    + 'Then rms per line for the ratios and for rho_w, and alpha2, the fraction of the H2 '
    + 'profile present.');
}

// 10. Big numbers ------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Fraction of the H2 profile present, α₂');
  const st = D.task7_stats;
  const vals = [['Ld / Ed', st.LdEd], ['Lu / Ed', st.LuEd], ['ρw', st.rho]];
  vals.forEach(([lab, [v, e]], i) => {
    const x = M + i * 4.15, w = 3.85;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.6, w, h: 2.9, rectRadius: 0.12,
      fill: { color: C.tint }, line: { color: C.tint } });
    s.addText(v, { x, y: 1.8, w, h: 1.4, fontFace: HEAD, fontSize: 66, bold: true,
      color: C.coral, align: 'center', valign: 'middle', margin: 0, isTextBox: true });
    s.addText('± ' + e + '   ' + lab, { x, y: 3.3, w, h: 0.8, fontFace: BODY, fontSize: 24,
      color: C.ink, align: 'center', valign: 'middle', margin: 0, isTextBox: true });
  });
  bullets(s, [
    '1 = the SRF width mismatch alone explains the feature',
    'Caveat: partly circular, since H2 is built from fits to the same spectra',
  ], { x: M, y: 5.0, w: W - 2 * M, h: 1.9 });
  s.addNotes('Weighted means over the SRF lines. What it establishes: the residual is what '
    + 'the width difference gives with equivalent width conserved (no Ring needed), the rho_w '
    + 'wiggle has the H2 shape, and interpolation has the wrong sign. The independent tests '
    + 'are the Phase 1 twin experiment and the Phase 3 correction.');
}

// 11. Task 7 table ------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'rms residual within ±5 nm of each line (task 7)');
  const y7 = table(s, ['Line', 'Measured Ld/Ed', 'H2 predicted', 'H1 predicted', 'α₂ Ld/Ed',
    'α₂ ρw'], D.task7, { x: M, y: 1.4, w: W - 2 * M }, [2.5, 2.1, 1.9, 1.9, 1.965, 1.965]);
  note(s, 'H2 matches the measured residual; H1 falls 10–100× short', y7 + 0.3);
  s.addNotes('Relative residual after a linear continuum through the pixels 5-7 nm from the '
    + 'line. Measured quantities are exactly those the processor divides: Ld interpolated in '
    + 'time over Ed linearly interpolated onto the L grid.');
}

// 12. Template SRF (3b) ------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'The SRF itself, from TSIS-1 HSRS (task 3b)');
  const im = image(s, 'veit_template_fwhm', { x: M, y: 1.3, w: 6.6, h: 5.95 });
  const x0 = M + im.w + 0.4, w0 = W - M - x0;
  table(s, ['λ (nm)', 'E', 'Ld', 'Ld − E'], D.task3b, { x: x0, y: 1.5, w: w0 },
    [w0 * 0.2, w0 * 0.2, w0 * 0.2, w0 * 0.4]);
  bullets(s, [
    'E ≈ 2.3 nm; Ld 0.5 nm wider in the blue',
    'Ld − E offset +0.045 ± 0.011 nm',
  ], { x: x0, y: 5.0, w: w0, h: 2.1 });
  s.addNotes('HSRS v2 in air convolved with a Gaussian SRF, fitted in 27 windows of 10 nm '
    + 'over 390-880 nm with telluric bands excluded. Constrained over 390-680 and 850-870 nm '
    + 'only. The empirical widths were 0.3-0.5 nm too large. These models are shipped as '
    + 'hypernet/data/veit_srf_model.json for Phase 1.');
}

// 13. Gate status --------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Gate G0 on this one sequence');
  const cards = [
    ['a', 'E = L?', 'No: H2 stays', 'Ld − E = 0.5 nm at 400–550 nm'],
    ['b', 'SRF stable?', 'Needs the data', 'Stable within one sequence; Phase 0b tests '
      + 'SZA, sky, season, recalibration'],
    ['c', 'Offset?', 'Relative: no', 'L − E = +0.045 nm; the absolute scale '
      + 'waits on air vs vacuum'],
  ];
  cards.forEach(([b, q, verdict, why], i) => {
    const x = M + i * 4.15, y = 1.5, w = 3.85, h = 4.2;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12,
      fill: { color: C.tint }, line: { color: C.tint } });
    badge(s, x + 0.3, y + 0.3, 0.8, b);
    s.addText(q, { x: x + 1.25, y: y + 0.3, w: w - 1.45, h: 0.8, fontFace: HEAD, fontSize: 24,
      bold: true, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(verdict, { x: x + 0.3, y: y + 1.35, w: w - 0.6, h: 0.8, fontFace: HEAD,
      fontSize: 28, bold: true, color: C.coral, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(why, { x: x + 0.3, y: y + 2.3, w: w - 0.6, h: 1.7, fontFace: BODY, fontSize: 22,
      color: C.ink, valign: 'top', margin: 0, isTextBox: true });
  });
  note(s, 'Criteria: (a) FWHM_E = FWHM_L; (b) spread < 0.2 nm; (c) offset > 0.1 nm', 6.1);
  s.addNotes('G0(a): if FWHM_E = FWHM_L on most instruments, drop H2. G0(b): spread below '
    + '0.2 nm across SZA, sky, season and no unexplained jump at recalibration means a '
    + 'per-instrument SRF table. G0(c): offsets above 0.1 nm mean a wavelength recalibration '
    + 'first. On the air scale all three channels are within 0.1 nm.');
}

// 14. Next steps ---------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addText('Next', { x: M + 0.3, y: 0.6, w: 12, h: 0.9, fontFace: HEAD, fontSize: 40,
    bold: true, color: C.white, margin: 0, isTextBox: true });
  s.addText([
    'Run Phase 0b on the ~150 requested sequences when they arrive',
    'Phase 1 twin experiment with the template SRFs',
    'From Kevin: cal files, lab line-spread data, air or vacuum?',
    'Look into a possible additive signal in the red of Lu',
  ].map((t, i, a) => ({ text: t, options: { bullet: true, breakLine: i < a.length - 1 } })),
  { x: M + 0.3, y: 1.8, w: 12, h: 4.8, fontFace: BODY, fontSize: 26, color: C.white,
    paraSpaceAfter: 18, valign: 'top', margin: 0, isTextBox: true });
  s.addNotes('The data request has been sent. The delivery runs through phase0b_index, '
    + 'phase0b_fit_all, phase0b_stability, phase0b_budget and phase0b_calfiles, then gate G0.');
}

pres.writeFile({ fileName: path.join(__dirname, 'wiggles_phase0_report.pptx') })
  .then(f => console.log('wrote', f));
