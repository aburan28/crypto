# 05-qr-window: PROTOCOL (frozen before data)

Pair-46 factor-base campaign, item 05. Written 2026-10-15 before any data was generated.
Not edited after the run starts; changes go in a dated amendment at the bottom.

## Question

The twist-aware quadratic-character filter (FB-QR-WINDOW, `research/factor_base_lab.md`
in crypto-autoresearcher) keeps only window x with chi(f(x)) = 1, f(x) = x^3+ax+b,
chi the Legendre symbol mod p. This item is the campaign's **positive control**: the
predicted effect is the ~2x point-density gain from discarding twist x-values, and
nothing more (identical relations, half the candidate columns). Does the instrument
detect exactly that?

## Claim classes

- **Exact** (proved, then checked by computation):
  - E1. If #E is odd (no 2-torsion, f has no root), #{x in F_p : chi(f(x)) = 1} = (#E-1)/2,
    and #{chi(f(x)) = -1} = (#E'-1)/2 for the quadratic twist E', #E + #E' = 2p+2.
    For pair 46 (crater and floor both have #E = N prime): (N-1)/2 liftable x, global
    density gain 2p/(N-1).
  - E2. Every E-point with x < B has chi(f(x)) = 1 (r = 0), so the plain window and the
    qr_window admit exactly the same E-points: every decomposition found by one is found
    by the other, target by target. The filter changes the candidate count, not the relations.
  - E3. The twist-only filter chi(f(x)) = -1 admits no E-point: zero relations.
- **Empirical**: the measured window densities (M1) at pair-46 size, the density of the
  literal spec reading chi(x) in {0,1} (predicted ~1/2, no gain), random-half control, and the
  toy-scale relation rates (M2) of the reduced arms.

## Arms (one instrument, one membership oracle per arm, same windows, same targets)

- A0 plain: V = {x in W}.
- A1 qr_window (twist-aware): V = {x in W : chi(f(x)) = 1}.
- A2 spec literal: V = {x in W : chi(x) in {0,1}} (the factor_base_lab.md formula read
  as the character of x itself). Predicted no density gain (negative control 1).
- A3 random-half: V = {x in W : splitmix64(x xor seed) even}. Negative control 2,
  same thinning rate as A1 with no arithmetic content.
- A4 twist-only: V = {x in W : chi(f(x)) = -1}. Instrument self-test: M1 must be 0 exactly.

## Metrics

- M1 point density: #{x in V : x is the x-coordinate of an E-point} / |V|.
  Gain G(arm) = M1(arm)/M1(A0). Prediction: G(A1) = |W|/#{x in W: chi(f)=1} ~ 2;
  G(A2), G(A3) ~ 1; M1(A4) = 0.
- M1' coverage: fraction of liftable x in W kept by the arm (A1: 1 exactly; A2, A3 ~ 1/2).
- M2 (toy only): Semaev m=2 relation rate, fraction of random targets R with
  R = ±P1 ± P2, x(P1), x(P2) in V; plus mean decompositions per target (D).
  Prediction: A1 identical to A0 per target (exact E2); A2, A3 lower (~ 1/4 of D);
  A4 = 0.

## Sizes and seeds

- Pair 46, both curves (crater, floor): first verify N·P = O on 3 random points per curve.
  - Canonical small-x window W = [0, B), B = ceil(sqrt(2N)) (~1.09e8, full enumeration).
  - 32 random windows of length 2^16 at uniform offsets in [0, p - 2^16), seed 46005.
- Toy: curves from `input_curves_from_04.csv` (copied verbatim from item 04,
  sha256 88c26329...7c12): all 100 F-arm curves (family p=5t^2+19t+19, j=-96^3, prime N,
  20..30 bits) and the first 100 N2-arm rows (random prime-order curves). For each:
  N·P = O check on 3 points; B = ceil(sqrt(2N)); T = 40 random targets R = k·G per curve,
  shared by all arms; seed 46005 xor curve index.
  Full-field exact check of E1 on every toy curve with p < 2^23.
- Intervals: Wilson 95% for single proportions (M1 on the canonical window); for ratios
  across windows/curves, percentile bootstrap (2000 resamples, seed 46005) over windows
  (pair 46) or curves (toy, paired by curve).

## Decision rule (fixed in advance)

- **Positive control detected** iff the 95% interval of G(A1) on pair 46 (both the canonical
  window and the 32-window bootstrap) lies inside [1.9, 2.1] and excludes 1, AND the A2 and
  A3 gain intervals contain 1, AND M1(A4) = 0 exactly, AND every toy curve has A1 relation
  sets identical to A0 (E2), AND the E1 count matches on every exact-checked curve.
- "Nothing more": the A1/A0 ratio of M2 and D equals 1 exactly on every curve.
- Any failure of an exact check is reported as an instrument bug, not a result.

## Stop rule

Fixed sizes above, one run, no adaptive extension. If the wall-clock budget (~40 min, 1 core)
is exceeded, cut the toy curve set from the top bit sizes down and report the cut.

## Implementation

Rust (AGENTS.md "no Python" for results), `src/main.rs`, single thread. Python only
plots recorded values (`plot.py`) and lays out the PDF (`make_pdf.py`, copied from item 04).
