# Exact F5 rank-table width screen

## Hypothesis and frozen inputs

Build on fixed cut 19 at `676f2941d`. The support certificate's exact
rank groups use row-basis elimination with eight Gray-code tables per
block, each table holding up to 256 entries. Reducing the pattern width
may save table construction work enough to lower total certificate XORs.
Freeze the four n24, m24, degree-4 planted systems with seed XORs `0`,
`badc0de1`, `5eed2026`, `f5c02a28`; all seven F5 cases, one Rayon thread,
the same direct packed-row and unpack settings, and one release binary.
Compare explicit rank pattern widths 5, 6, 7, and 8 in separate
processes. Width 8 is the unchanged reference.

## Exactness and accounting

Add an experiment-only `KIC_GF2_ROW_BASIS_BITS` selector to
`rank_row_basis_counted` while leaving `echelon_counted`, RREF and the
default row-basis width unchanged. Use a native Rust screen driver that
preserves every call, failure and timeout, source/binary hashes, host,
rank, canonical row-space and raw fingerprints, route, counted word
XORs, exclusive phases and complete-call wall time. Each width must
produce all seven cases; the primary must match rank, canonical row
space, F5 criterion, row and column counts, and original-row output.
All smaller cases must match every nontiming output field except the
explicit requested width label. A failed certificate or different row
space rejects that width.

The primary decision unit is counted 64-bit word XORs, including all
table construction and row clears. Promote a width only if it certifies
all four seeds and uses at most 90% of width 8's reduction XORs on
each. Preserve nonpromoting Apple ARM64 phase and complete-call times,
but do not choose a width by their noisy wall ratios. If none passes,
stop and retain the negative cells. A promoted width still needs the
unchanged four-seed physical isolated Linux complete-call median and
bootstrap lower bound both above 2× against selective echelon, plus
smaller-case and two-thread nonregression controls. This is a solver
stage study, not an IC online or rho speedup.
