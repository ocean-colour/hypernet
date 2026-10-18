// Build slides/wiggles_phase1.pptx (Phase 1 task 13: the twin experiment).
// Numbers come from wiggles_phase1_data.json, written by wiggles_phase1_data.py
// from the committed Phase 1 tables.  Rule for this deck: no text below 20 pt.
//
//   python slides/wiggles_phase1_data.py
//   NODE_PATH=$(conda run -n slides npm root -g) conda run -n slides node slides/build_wiggles_phase1.js
const path = require('path');
const pptxgen = require('pptxgenjs');
const D = require(path.join(__dirname, 'wiggles_phase1_data.json'));

const C = {
  navy: '0B2545', teal: '13505B', tint: 'E6F0F1', coral: 'D9542C', ink: '1F2A30',
  muted: '4A5A60', white: 'FFFFFF', line: 'C9D8DA', green: '1B7F5A', pale: 'CFE3E6',
};
const HEAD = 'Cambria', BODY = 'Calibri';
const W = 13.333, H = 7.5, M = 0.5;

const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';
pres.title = 'WATERHYPERNET wiggles: Phase 1, the twin experiment';
pres.author = 'J. Xavier Prochaska';

function title(s, text) {
  s.addText(text, { x: M, y: 0.3, w: W - 2 * M, h: 0.85, fontFace: HEAD, fontSize: 34,
    bold: true, color: C.navy, margin: 0, valign: 'middle', isTextBox: true });
}

function badge(s, x, y, d, label, fill, size) {
  s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: fill || C.teal },
    line: { color: fill || C.teal } });
  s.addText(label, { x, y, w: d, h: d, fontFace: HEAD, fontSize: size || 24, bold: true,
    color: C.white, align: 'center', valign: 'middle', margin: 0, isTextBox: true });
}

function bullets(s, items, box, size) {
  s.addText(items.map((t, i) => ({ text: t, options: { bullet: true,
    breakLine: i < items.length - 1 } })),
  Object.assign({ fontFace: BODY, fontSize: size || 22, color: C.ink, valign: 'top',
    paraSpaceAfter: 12, margin: 0.05, isTextBox: true }, box));
}

function image(s, key, box) {
  const im = D.images[key];
  let w = box.w, h = w / im.aspect;
  if (h > box.h) { h = box.h; w = h * im.aspect; }
  const x = box.center ? box.x + (box.w - w) / 2 : box.x;
  s.addImage({ path: im.path, x, y: box.y, w, h });
  return { x, w, h };
}

const ROW_H = 0.5;
function table(s, header, rows, box, colW, opts) {
  opts = opts || {};
  const hdr = header.map(t => ({ text: t, options: { bold: true, color: C.white,
    fill: { color: C.teal }, align: 'center', valign: 'middle' } }));
  const body = rows.map((r, i) => r.map((t, j) => ({ text: t, options: {
    fill: { color: i % 2 ? C.white : C.tint }, align: j === 0 ? 'left' : 'center',
    valign: 'middle', bold: j === 0 || (opts.boldCol === j),
    color: opts.boldCol === j ? C.green : C.ink } })));
  s.addTable([hdr].concat(body), Object.assign({ fontFace: BODY, fontSize: 20,
    color: C.ink, border: { type: 'solid', pt: 0.75, color: C.line }, colW,
    rowH: ROW_H, margin: [0.03, 0.1, 0.03, 0.1] }, box));
  return box.y + (rows.length + 1) * ROW_H;
}

function note(s, text, y, h) {
  s.addText(text, { x: M, y, w: W - 2 * M, h: h || 0.5, fontFace: BODY, fontSize: 20,
    italic: true, color: C.muted, margin: 0, isTextBox: true });
}

// A callout card: a big number over a short label.
function stat(s, x, y, w, h, big, label, col) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.1,
    fill: { color: C.tint }, line: { color: C.tint } });
  s.addText(big, { x: x + 0.2, y: y + 0.12, w: w - 0.4, h: 0.85, fontFace: HEAD,
    fontSize: 36, bold: true, color: col || C.coral, margin: 0, valign: 'middle',
    isTextBox: true });
  s.addText(label, { x: x + 0.2, y: y + 0.95, w: w - 0.4, h: h - 1.05, fontFace: BODY,
    fontSize: 20, color: C.ink, margin: 0, valign: 'top', isTextBox: true });
}

// Figure on the left, callout cards stacked on the right.
function figureWithStats(s, key, stats, figW) {
  figW = figW || 8.4;
  const im = image(s, key, { x: M, y: 1.3, w: figW, h: H - 1.3 - 0.4 });
  const x = M + im.w + 0.3, w = W - M - x, n = stats.length;
  const gap = 0.25, h = (H - 1.3 - 0.4 - gap * (n - 1)) / n;
  stats.forEach(([big, label, col], i) => stat(s, x, 1.3 + i * (h + gap), w, h, big,
    label, col));
}

const T9 = D.t9, T10 = D.t10, T10D = D.t10d, T11 = D.t11;

// 1. Title ---------------------------------------------------------------------------
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addText('Spectral wiggles in WATERHYPERNET', { x: M + 0.3, y: 1.6, w: 12, h: 1.0,
    fontFace: HEAD, fontSize: 44, bold: true, color: C.white, margin: 0, isTextBox: true });
  s.addText('Phase 1: a twin experiment. Which interpolation removes which wiggle, '
    + 'and what does it do to real features?', { x: M + 0.3, y: 2.75, w: 11.8, h: 1.3,
    fontFace: BODY, fontSize: 26, color: C.pale, margin: 0, isTextBox: true });
  s.addText('Synthetic HYPSTAR scenes · 25 sky and water states · 2 instruments · '
    + '16 cases · 5 methods', { x: M + 0.3, y: 4.4, w: 12, h: 0.5, fontFace: BODY,
    fontSize: 22, color: C.white, margin: 0, isTextBox: true });
  s.addText('J. Xavier Prochaska · UC Santa Cruz · 5 October 2026', { x: M + 0.3, y: 5.9,
    w: 12, h: 0.5, fontFace: BODY, fontSize: 20, color: C.pale, margin: 0,
    isTextBox: true });
  s.addNotes('Phase 1 tests the interpolation methods on synthetic data where the truth is '
    + 'known, before any correction touches real WATERHYPERNET spectra. It closes with '
    + 'gate G1, which sets the default method for Phase 2.');
}

// 2. Why a twin experiment -------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Why a twin experiment first');
  const cards = [
    ['H1', 'Grid offset (Kevin)', 'E and L sample different wavelengths; linear interpolation '
      + 'of Ed onto the L grid misplaces the Fraunhofer lines', C.teal],
    ['H2', 'SRF mismatch (ours)', 'E sees the lines through a narrower SRF than L, so the '
      + 'lines survive in L/E however Ed is interpolated', C.coral],
  ];
  cards.forEach(([b, head, body, col], i) => {
    const x = M + i * 6.3, y = 1.4, w = 5.9, h = 3.0;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12,
      fill: { color: C.tint }, line: { color: C.tint } });
    badge(s, x + 0.3, y + 0.3, 0.85, b, col);
    s.addText(head, { x: x + 1.35, y: y + 0.3, w: w - 1.6, h: 0.85, fontFace: HEAD,
      fontSize: 26, bold: true, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(body, { x: x + 0.3, y: y + 1.35, w: w - 0.6, h: h - 1.5, fontFace: BODY,
      fontSize: 22, color: C.ink, valign: 'top', margin: 0, isTextBox: true });
  });
  bullets(s, [
    'Real ρw has no truth to compare against, and real features (fluorescence, Raman, '
      + 'gas-band peaks) look like wiggles',
    'A twin has a known truth, so every method can be scored, and its failure modes '
      + 'mapped, before Phase 3 applies it to data',
  ], { x: M, y: 4.75, w: W - 2 * M, h: 2.2 });
  s.addNotes('Phase 0 measured the SRFs and found H2 on the narrow-E instruments. Phase 1 '
    + 'asks which correction fixes it, how robust that correction is, and how much '
    + 'reduction Phase 3 should expect on real data.');
}

// 3. Three methods ---------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Three ways to put Ed on the L grid');
  const rows = [
    ['linear', 'E(λL) = (1 − w) Eᵢ + w Eᵢ₊₁', 'The processor today (paper eqs. 11–12)',
      C.coral],
    ['ruddick2023', 'E(λL) = Emod⊗SRF_E(λL) × [interpolated E / Emod⊗SRF_E]',
      'Paper eq. 14: corrects the grid offset (H1)', C.teal],
    ['srf', 'E(λL) = Emod⊗SRF_L(λL) × [interpolated E / Emod⊗SRF_E]',
      'This work: L SRF in front, so H1 and H2', C.green],
  ];
  rows.forEach(([name, eq, what, col], i) => {
    const y = 1.45 + i * 1.65;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: M, y, w: W - 2 * M, h: 1.4,
      rectRadius: 0.1, fill: { color: C.tint }, line: { color: C.tint } });
    s.addText(name, { x: M + 0.3, y, w: 2.6, h: 1.4, fontFace: HEAD, fontSize: 26,
      bold: true, color: col, valign: 'middle', margin: 0, isTextBox: true });
    s.addText([{ text: eq, options: { breakLine: true, fontFace: 'Cambria Math' } },
      { text: what, options: { color: C.muted, italic: true } }],
    { x: M + 3.0, y, w: W - 2 * M - 3.3, h: 1.4, fontFace: BODY, fontSize: 22,
      color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
  });
  note(s, 'Emod: a 0.01 nm model irradiance.  With SRF_L = SRF_E, srf is eq. 14; with a '
    + 'flat Emod, it is linear.  Cubic and sinc are nulls.', 6.45, 0.8);
  s.addNotes('The bracket is the measured-over-model ratio, which is smooth if the model is '
    + 'good, so it is interpolated linearly. Eq. 14 multiplies by the model seen through the '
    + 'E SRF, so it restores E-width lines on the L grid; dividing L by that still leaves '
    + 'the width mismatch. The srf form uses the L SRF in front.');
}

// 4. How the twin is built ---------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'How the twin is built');
  const steps = [
    ['1', 'Irradiance', 'TSIS-1 HSRS × HITRAN O₂/H₂O × O₃ on a 0.01 nm grid, '
      + 'direct and diffuse air masses'],
    ['2', 'Smooth fields', 'OSOAA at 5 nm: direct/diffuse Ed, sky Ld, Lw; SZA 30/50/70 '
      + '(+ a VEIT-like 40°), two aerosol loads, four waters'],
    ['3', 'Scene', 'Ed, Ld, Lw at 0.01 nm, with chlorophyll fluorescence and water Raman '
      + 'as controls'],
    ['4', 'Instrument', 'Gaussian SRFs from the Release 2 fits, the real 1536/1538-pixel '
      + 'grids, VEIT per-scan noise'],
    ['5', 'Score', 'ρw from each method; rms of ρw″ − ρw″_true within ±5 nm of ten lines'],
  ];
  steps.forEach(([n, head, body], i) => {
    const y = 1.35 + i * 1.12;
    badge(s, M, y + 0.05, 0.8, n);
    s.addText(head, { x: M + 1.05, y, w: 2.6, h: 0.9, fontFace: HEAD, fontSize: 24,
      bold: true, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(body, { x: M + 3.7, y, w: W - 2 * M - 3.7, h: 0.9, fontFace: BODY,
      fontSize: 22, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
  });
  s.addNotes('OSOAA has no gas absorption, so the split into a high-resolution irradiance '
    + 'times smooth radiative-transfer fields is clean. The full twin (25 scenes x 2 '
    + 'instruments x 16 cases x 5 methods) runs in about a minute on 12 cores. ρw″ is the '
    + 'second difference with h = 1 nm, as in Ruddick et al. (2023).');
}

// 5. Scene check -------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'A twin scene, and the features we must not erase');
  figureWithStats(s, 'scene_check', [
    ['683 nm', 'fluorescence: a 25 nm peak with no Fraunhofer structure'],
    ['O₂-A, O₂-B', 'Raman and fluorescence are emitted below the atmosphere, so ρw peaks '
      + 'in the gas bands'],
  ], 7.9);
  s.addNotes('Top: the 0.01 nm Ed and Ld at the G band and O2-A. Bottom: true ρw for the '
    + 'four waters and the two controls. A correction that smooths these peaks away would '
    + 'be removing real signal.');
}

// 6. Instruments and cases ---------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Two instruments, sixteen cases');
  const inst = [
    ['122304', 'narrow-E: FWHM_L − FWHM_E ≈ 0.5 nm (also 122302, 120242)', C.coral],
    ['122305', 'E ≈ L within 0.07 nm (also 121222, 121231, 122303)', C.teal],
  ];
  inst.forEach(([name, body, col], i) => {
    const y = 1.35 + i * 0.85;
    s.addText([{ text: name + '  ', options: { bold: true, color: col,
      fontFace: HEAD } }, { text: body }], { x: M, y, w: W - 2 * M, h: 0.7,
      fontFace: BODY, fontSize: 22, color: C.ink, valign: 'middle', margin: 0,
      isTextBox: true });
  });
  table(s, ['Case', 'What changes'], [
    ['(i)', 'grid offset only, equal SRFs: H1 alone'],
    ['(ii) / (iii)', 'SRF mismatch alone / offset + mismatch: H2'],
    ['(iv)', 'E wavelength error: shift ±0.1, ±0.3 nm; stretch ±0.1 nm'],
    ['(v)', 'correction told FWHM_E or FWHM_L ±0.3 nm'],
    ['(vi)', 'wrong Emod: SZA + 10°, water vapour × 2, low aerosol'],
    ['(vii)', 'VEIT per-scan noise, mean of 6 scans'],
  ], { x: M, y: 3.2, w: W - 2 * M }, [2.4, W - 2 * M - 2.4]);
  s.addNotes('The SRFs are the clear-sky medians of the Release 2 template fits (Phase 0 '
    + 'task 8c). Each case runs on all 25 scenes. Cases iv to vi ask what happens when the '
    + 'correction is told something wrong; vii asks whether the gain survives noise.');
}

// 7. rho'' at the lines ------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'What the methods do at Ca H/K, Hα and O₂-A');
  figureWithStats(s, 'twin_rho2_lines', [
    [T9.ratio_lin_true, 'linear error vs the true line structure, case (iii), 122304'],
    ['srf = truth', 'green lies on black in (iii); ruddick2023 tracks linear'],
    ['(i) small', 'offset alone: ' + T9.lin_i + ', ' + T9.offset_frac
      + ' of the mismatch error'],
  ], 8.2);
  s.addNotes('Top row: case iii (offset plus mismatch). Bottom row: case i (offset only). '
    + 'In case i every model-based method, and even cubic, removes the error; in case iii '
    + 'only srf does.');
}

// 8. Reduction overview ------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Every case, every method');
  const im = image(s, 'twin_reduction', { x: M, y: 1.25, w: W - 2 * M, h: H - 1.25 - 0.95,
    center: true });
  note(s, 'Median rms ρw″ error near the lines (h = 1 nm).  Grey bars: srf away from the '
    + 'lines.', 1.25 + im.h + 0.15);
  s.addNotes('Top: 122304, narrow-E. Bottom: 122305, E about L. Green (srf) drops by '
    + 'orders of magnitude in i to iii and vi; the cases that hurt it are iv (wavelength '
    + 'error) and v (wrong SRF width).');
}

// 9. 122304 table -------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Narrow-E (122304): srf fixes it, eq. 14 does not');
  const bottom = table(s, ['Case', 'linear error', 'ruddick2023', 'srf'], D.t9_rows,
    { x: M, y: 1.3, w: W - 2 * M }, [4.6, 2.6, 2.6, 2.533], { boldCol: 3 });
  note(s, 'Reductions relative to linear; error = median rms ρw″ − ρw″_true near the '
    + 'lines, h = 1 nm.  True line structure: ' + T9.true_line + '.', bottom + 0.15, 0.8);
  s.addNotes('Ruddick2023 cannot fix a width mismatch: it removes 0 % in ii and makes iii '
    + 'slightly worse. srf removes essentially all of it, is robust to a wrong model '
    + 'irradiance, and is sensitive to the E wavelength scale and the SRF widths. With '
    + 'noise the line error falls to the away-from-lines level: line/away ' + T9.lo_lin_vii
    + ' for linear, ' + T9.lo_srf_vii + ' for srf.');
}

// 10. 122305 -------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'E ≈ L (122305): little to fix, and easy to break');
  const st = [
    [T9.lin5_iii, 'linear error in (iii): already below the true line structure ('
      + T9.true_line + ')', C.teal],
    [T9.srf5_iii, 'srf reduction when told the right SRFs', C.green],
    [T9.srf5_v, 'srf error vs linear when an SRF width is 0.3 nm off', C.coral],
    [T9.rud5_iii, 'ruddick2023 in (iii): linear\'s smoothing hid part of the mismatch',
      C.coral],
  ];
  st.forEach(([big, label, col], i) => stat(s, M + (i % 2) * 6.25, 1.35
    + Math.floor(i / 2) * 2.15, 6.0, 1.9, big, label, col));
  note(s, 'H1 alone gives ' + T9.rho5_i + ' rms ρw error near the lines, offset + residual '
    + 'mismatch ' + T9.rho5_iii + '.  The 0.5–2 % floor seen on these instruments needs '
    + 'something else.', 5.75, 1.0);
  s.addNotes('On E about L instruments the correction gains nothing measurable and can do '
    + 'harm if the SRFs are slightly wrong. This is why Phase 2 applies srf only above a '
    + 'width-difference threshold. The twin also says the grid offset (H1) cannot explain '
    + 'the unexplained floor Phase 0 found on these instruments: candidates are the Ring '
    + 'effect in Ld, SRF shape, or the wavelength scale. The metrics differ, so this is '
    + 'indicative.');
}

// 11. Controls -----------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'srf restores the controls; the others distort them');
  image(s, 'twin_controls', { x: M, y: 1.3, w: 7.0, h: H - 1.3 - 0.4 });
  const x = M + 7.3, w = W - M - x;
  table(s, ['122304 (iii)', 'fluor.', 'Raman'], T10.rows, { x, y: 1.3, w },
    [2.25, 1.43, 1.4], { boldCol: 0 });
  s.addText('Change in each control\'s ρw″, as a fraction of the control itself',
    { x, y: 3.45, w, h: 0.9, fontFace: BODY, fontSize: 20, italic: true, color: C.muted,
      margin: 0, isTextBox: true });
  bullets(s, [
    'No method moves a control by more than ' + T10.max_over_noise + ' of the noise floor',
    'But the controls are themselves below one sequence\'s noise (' + T10.amp_fl
      + ', ' + T10.amp_raman + ')',
  ], { x, y: 4.45, w, h: 2.6 }, 20);
  s.addNotes('Each scene is composed with both controls, without fluorescence and without '
    + 'Raman; a method sees a control as the difference. Ed is the same, so only the '
    + 'irradiance error at the lines distorts it, and both controls peak in the O2 bands '
    + 'where that error is largest. The G1 control wording was sharpened for this reason.');
}

// 12. Degradation -----------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'How srf degrades');
  const im = image(s, 'twin_degradation', { x: M, y: 1.2, w: W - 2 * M, h: H - 1.2 - 0.85,
    center: true });
  note(s, 'Case (iii) with one error swept; medians over 25 scenes.  Top 122304, bottom '
    + '122305.', 1.2 + im.h + 0.1);
  s.addNotes('The model irradiance can be quite wrong (SZA off by 10 degrees, water vapour '
    + 'off by a factor 3) because only its line shapes matter. The E wavelength scale and '
    + 'the SRF widths are what count. On 122304 srf is never worse than linear; on 122305 '
    + 'it is, once an SRF is a few tenths of a nm off.');
}

// 13. Tolerances ------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'What srf needs to keep 80 % (122304)');
  const st = [
    ['≤ 0.05 nm', 'E wavelength scale vs the model (shift: ' + T10D.shift05 + '; 0.1 nm → '
      + T10D.shift10 + ', 0.3 nm → ' + T10D.shift30 + ')', C.coral],
    ['≲ 0.07 nm', 'FWHM_E error (0.05 nm → ' + T10D.e05 + ', 0.1 nm → ' + T10D.e10 + ')',
      C.coral],
    ['≲ 0.1 nm', 'FWHM_L error (0.1 nm → ' + T10D.l10 + ')', C.teal],
    ['≥ ' + T10D.emod_min, 'kept with the Emod at SZA ± 10° or water vapour × 0.5–3',
      C.green],
  ];
  st.forEach(([big, label, col], i) => stat(s, M + (i % 2) * 6.25, 1.35
    + Math.floor(i / 2) * 2.3, 6.0, 2.05, big, label, col));
  note(s, 'Phase 0: SRFs stable to 0.01–0.06 nm (meets the width tolerances); a dispersion '
    + 'trend > 0.1 nm at the ends (does not meet the wavelength one).', 6.15, 0.9);
  s.addNotes('The SRF tolerances are met by per-instrument tables. The wavelength tolerance '
    + 'is the binding one, so Phase 2 adds a fit of a shift and stretch of E against the '
    + 'reference solar spectrum (task 4b).');
}

// 14. Prediction figure ----------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Prediction for Phase 3');
  const im = image(s, 'twin_prediction', { x: M, y: 1.3, w: W - 2 * M, h: 4.3,
    center: true });
  bullets(s, [
    'Phase 3 has no truth, so G3 uses the observable excess √(rms²near − rms²away) of ρw″',
    'Per line, the gain follows signal-to-noise: Ca H/K best, H₂O and O₂-A noise-limited',
  ], { x: M, y: 1.3 + im.h + 0.2, w: W - 2 * M, h: H - 1.5 - im.h - 0.3 }, 20);
  s.addNotes('Case iii on the 122304 grids, with the E SRF narrower than L by 0 to 1 nm. '
    + 'Three scenarios: ideal, VEIT noise, and noise plus 0.05 nm errors in the wavelength '
    + 'scale and FWHM_E. The perfect-correction curve is the ceiling: real line structure '
    + 'leaves an excess of its own.');
}

// 15. Prediction table ------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'The yardstick for gate G3');
  const bottom = table(s, ['srf, one sequence, h = 1 nm'].concat(T11.dfwhm.map(
    x => x + ' nm')), T11.rows, { x: M, y: 1.35, w: W - 2 * M },
  [4.333].concat(T11.dfwhm.map(() => 8.0 / 6)));
  note(s, 'Reduction of the observable excess vs linear, by FWHM_L − FWHM_E.  At the '
    + 'narrow-E 0.5 nm: 60–75 %.  Below 0.1 nm: little to gain.', bottom + 0.15, 0.8);
  const y0 = bottom + 1.15;
  s.addText('Per line at 0.5 nm (noise only):', { x: M, y: y0, w: W - 2 * M, h: 0.5,
    fontFace: HEAD, fontSize: 22, bold: true, color: C.ink, margin: 0, isTextBox: true });
  s.addText(T11.lines.map(([n, , r], i) => ({ text: n + ' ' + r + (i < T11.lines.length - 1
    ? '   ·   ' : ''), options: {} })), { x: M, y: y0 + 0.55, w: W - 2 * M, h: 1.0,
    fontFace: BODY, fontSize: 22, color: C.ink, margin: 0, valign: 'top', isTextBox: true });
  s.addNotes('Averaging N sequences lowers the noise and moves the result toward the '
    + 'noise-free curve, which is close to 100 %. h = 5 nm is dominated by real structure, '
    + 'so it is not a useful G3 metric.');
}

// 16. Gate G1 ---------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Gate G1: passed for narrow-E instruments');
  const g = [
    ['1', '> 80 % in (ii)–(iii)', 'srf ~100 % on both instruments', 'Pass', C.green],
    ['2', 'Graceful in (iv)–(vi)', 'narrow-E: yes, within the tolerances; E ≈ L: no',
      'Pass*', C.green],
    ['3', 'Controls kept', 'srf exact in (iii); never worse than linear', 'Pass', C.green],
    ['4', 'Does eq. 14 alone pass (iii)?', '0 % in (ii), −6 % in (iii)', 'No: H2 needed',
      C.coral],
  ];
  g.forEach(([n, crit, res, v, col], i) => {
    const y = 1.35 + i * 1.12;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: M, y, w: W - 2 * M, h: 0.95,
      rectRadius: 0.1, fill: { color: C.tint }, line: { color: C.tint } });
    badge(s, M + 0.15, y + 0.1, 0.75, n, col, 22);
    s.addText(crit, { x: M + 1.1, y, w: 4.1, h: 0.95, fontFace: HEAD, fontSize: 22,
      bold: true, color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(res, { x: M + 5.25, y, w: 4.6, h: 0.95, fontFace: BODY, fontSize: 20,
      color: C.ink, valign: 'middle', margin: 0, isTextBox: true });
    s.addText(v, { x: M + 9.9, y, w: 2.3, h: 0.95, fontFace: HEAD, fontSize: 22,
      bold: true, color: col, align: 'right', valign: 'middle', margin: 0,
      isTextBox: true });
  });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: M, y: 5.95, w: W - 2 * M, h: 1.0,
    rectRadius: 0.1, fill: { color: C.navy }, line: { color: C.navy } });
  s.addText('Phase 2 default: srf where FWHM_L − FWHM_E ≥ 0.15 nm, linear otherwise',
    { x: M + 0.3, y: 5.95, w: W - 2 * M - 0.6, h: 1.0, fontFace: HEAD, fontSize: 24,
      bold: true, color: C.white, valign: 'middle', margin: 0, isTextBox: true });
  s.addNotes('* conditional on the E wavelength scale to 0.05 nm and the SRF widths to '
    + '0.07 to 0.1 nm. The 0.15 nm threshold is where srf breaks even with linear under '
    + 'realistic knowledge errors. Full report: claude_prompts/wiggles/gate_G1.md.');
}

// 17. Next ------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, 'Open concerns and next steps');
  const cols = [
    ['Open concerns', [
      'E wavelength scale: the binding requirement',
      'E SRFs need the L1A delivery (L2 widths read ~0.1 nm wide)',
      'Gaussian SRFs only; Ring and the Lu NIR excess are outside the twin',
      'The E ≈ L floor is still unexplained',
    ], C.coral],
    ['Phase 2', [
      'Harden interpolate_ed_to_l, with the 0.15 nm rule',
      'Fit E\'s wavelength shift and stretch (task 4b)',
      'SRF table, punpy uncertainties, processor snippet',
      'Non-Gaussian SRF twin check (task 9b)',
    ], C.teal],
  ];
  cols.forEach(([head, items, col], i) => {
    const x = M + i * 6.3, w = 5.9;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.35, w, h: 4.6, rectRadius: 0.12,
      fill: { color: C.tint }, line: { color: C.tint } });
    s.addText(head, { x: x + 0.3, y: 1.5, w: w - 0.6, h: 0.7, fontFace: HEAD, fontSize: 26,
      bold: true, color: col, margin: 0, valign: 'middle', isTextBox: true });
    bullets(s, items, { x: x + 0.3, y: 2.35, w: w - 0.6, h: 3.5 }, 22);
  });
  s.addNotes('Phase 3 then applies the corrections to the requested sequences and judges '
    + 'them against the yardstick on slide 15 (gate G3). Phase 0 tasks 9 to 14, and so '
    + 'gate G0, still wait on Kevin\'s L1A delivery.');
}

pres.writeFile({ fileName: path.join(__dirname, 'wiggles_phase1.pptx') })
  .then(f => console.log('wrote ' + f));
