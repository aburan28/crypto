#set page(paper: "a4", margin: 17mm)
#set text(font: "Arial", size: 10pt)
= Fourth BMSS candidate control
2026-10-15 · 96/96 exact outputs; 384/384 counters unchanged

#image("CONTROL.svg", width: 100%)

The small-precision Weierstrass recurrence now takes a dedicated direct-division path. The batch-inverted path removes per-iteration Option dispatch and unused denominator conversion. This timing hypothesis follows the retained 256-bit/degree-7 Atkin regressions in the preceding AMD campaign; its timing benefit remains unverified.

All eight BMSS constructors at field widths 192/224/256/511 and prime degrees 7/31/101 reproduce the original and previous-rung full coefficients, codomain and kernel. Each case retains independent Kohel and public point-transport checks. All four generic Field counters are identical to the previous rung in all 96 cells. Inversion internals remain opaque. Setup, sigma, reference generation and verification stay outside constructor counts. These are synthetic supersingular field-width fixtures. Named P192/P224 whole-search and replay costs remain open.

The standalone release library, binary, integration and example tests pass; strict all-target Clippy and formatting pass. The original tested crate and continuation checkout compare byte for byte. Both comparison run hashes are checked before generating this report. The initial accidental selection of an older header-only baseline was rejected by the strict reader and remains documented in the frozen run record.

The next native wave retains separate original/candidate and previous/candidate twenty-command campaigns, each with five A/A and five interleaved A/B pairs. No timing ratio is established by this control run.

Source revision:
#text(size: 7pt)[9a55ce6d22b9473c1d2afbf398bfdf51ae22b544]
Binary SHA256:
#text(size: 7pt)[1e375ac4a49d317ca8375b6a7cc22c9c4b99b06102e41afea1847c00c454c8ab]
Incremental comparison SHA256:
#text(size: 7pt)[334a2dc8c64c6e34aa7cbfc7e83b651f72ab15995dc3d79d822d741dd8e1c520]
