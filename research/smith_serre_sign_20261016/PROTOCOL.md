# Protocol: the Prop. 4 sign follows #H(F_p)

Frozen before any `p = 11` data of this study. Follows up
[`problems/smith-kernel-isogeny-class-walk.md`](../../problems/smith-kernel-isogeny-class-walk.md):
`full_success` clusters within `L`-classes beyond isomorphism-class
granularity ([2026-10-15](../smith_iso_granularity_20261015/REPORT.pdf)), and
no 2-adic invariant of `J_H` explains it
([2026-10-16](../smith_class_invariant_20261015/REPORT.pdf)).

## What the exploration at p = 7 found (exploratory, not evidence here)

Scratch analysis of `p = 7` (same instrument as below, run before this
freeze) showed:

1. **Twist complementarity.** Every `L`-class has as many Prop. 2 passes as
   its twist class `L(−T)`, and Prop. 4 flips under the quadratic twist
   except on degenerate subgroups where it fails on both sides (Smith's
   construction fails at small `p`). This matches Serre's obstruction: the
   quotient `J_H/S` is either `Jac(X)` (Prop. 4 passes, `X` in class `L`)
   or its quadratic twist (`X` in class `L(−T)`).
2. **Sign law.** Among non-degenerate Prop. 2 passes, the Prop. 4 pass rate
   rises with `N1 = #H(F_p) = p + 1 + a1` from 0.016 at `N1 = 0` to 0.984 at
   `N1 = 16`. It is close to `N1 / (2p + 2)` (0.500 vs 0.500 at `N1 = 8`), the
   rate if `X` were chosen in proportion to its number of rational points
   between the two candidates (`#X(F_p) = N1` or `2p + 2 − N1`).
3. **Clustering.** Conditioning the iso-level null on `(stratum, N1, any
   Prop. 2 pass)` moves its median from 0.4579 to about 0.442 against
   observed `C = 0.4304`.

Since `N1` is an `L`-class invariant, the law would make Prop. 4 outcomes
correlated within a class, and empty classes would be classes with few
rational points or no Prop. 2 pass.

## Confirmatory measurements at p = 11 (and p = 7 rerun, labelled exploratory)

Instrument: `smith_iso_granularity --p P --out DIR --subgroups`, which also
writes per-isomorphism-class Prop. 2 / Prop. 4 subgroup counts
(`iso_sub_p<P>.csv.gz`), with controls C1–C3 of the 2026-10-15 protocol
plus C4 (counts constant on every orbit). C5 (exact twist complementarity)
is reported but expected to fail on degenerate subgroups. Its replacement
C5′ is that the curve-side plus twist-side Prop. 4 passes never exceed the
Prop. 2 passes. Analysis: `smith_iso_granularity summarize-serre`, in Rust.

- **H1 (sign law).** For each `N1`, `r(N1) = I(N1) / (I(N1) + I(2p+2−N1))`.
  Here `I(n)` is the model-weighted count of Prop. 4 passes on curves with
  `#H(F_p) = n`. Pass if both hold: `r` is non-decreasing in `N1` over bins
  holding ≥ 1% of the Prop. 2 passes, and `max |r(N1) − N1/(2p+2)| ≤ 0.10` over
  those bins.
- **H2 (clustering).** Run 999 permutations (seed 20261016) of `full_success`
  over isomorphism classes, within cells `(stratum, N1, any Prop. 2 pass)`.
  Verdicts:
  - observed `C` inside the central 95% → "N1 and Prop. 2 explain the
    clustering";
  - observed `C` below the central 95% → "partial: residual clustering
    remains".

  Report the share of the gap `U − C` closed, `(U − median)/(U − C)`, next
  to the 2026-10-15 stratum-only null.
- **Control (random coarse label).** Give each `L`-class a uniform random
  label with as many values as `N1` takes (seed 20261016), and rerun the H2
  null with that label in place of `N1`. It should close far less of the
  gap than `N1` does. If it closes ≥ 75% of what `N1` closes, H2's gain is
  attributed to coarse conditioning alone and H2 is reported as
  inconclusive.

## Scope and limits

`p = 7, 11`, degree-8 models only. H1 is an empirical law at two primes. The
`N1 / (2p+2)` form is a heuristic suggested by the data, not derived. If it
holds, the bias vanishes like `|a1|/(2p+2) = O(1/√p)` and Smith's fair coin
is recovered asymptotically. That is an extrapolation, not tested here.
Nothing here builds or prices an isogeny walk.
