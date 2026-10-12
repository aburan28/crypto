# ecbench calibration, follow-up: the two cells whose predictions failed

**Preregistered 2026-10-02, after round 1 and before this follow-up ran.**
Round 1 ([`PROTOCOL.md`](PROTOCOL.md), results in [`README.md`](README.md))
failed two of its predictions as written. Those failures stand. This
follow-up does not revise them; it tests a stated explanation on fresh
data.

## What failed in round 1

| cell | prediction | observed (24 runs) |
|---|---|---|
| `rho.plain` on `icv1-f2m41-tm2308219-7f48b14a` | `S/theory ∈ [0.85, 1.35]` | `0.792`; the 95 % interval `[0.688, 1.307]` of the mean contains the theory constant |
| `kangaroo.vow` on `icv1-fp26-tm1775-7e8fb6df` | `S/theory ∈ [0.8, 1.5]` | `0.601`; the 95 % interval `[0.684, 1.761]` of the mean excludes `2.0` |

## Explanation under test

1. **Rho:** the band was set on a point estimate without allowing for its
   sampling spread. A rho run's cost has a coefficient of variation of
   about 0.5, so a mean of 24 runs has about ±20 % at 95 %, wider than the
   band's lower margin. The cell is a band-design error, not an
   accounting error.
2. **Kangaroo:** `kangaroo.vow` is deterministic given its target (the
   seed is used only on a restart), so round 1's three rounds per target
   were replicates and its effective sample was 8 targets. Its cost is
   driven by the distance `|k − r/2|`. The 8 targets' mean distance was
   `0.119 r` against the uniform `0.25 r`, about −2.6σ. The derivation of
   targets is not biased: a native test over 20 000 draws
   (`workload::tests::planted_scalars_are_uniform`) passes a 1 %
   Kolmogorov–Smirnov test with mean distance `0.2499`.

## Design

- [`specs/followup-prime.json`](specs/followup-prime.json): the fp26 curve,
  arms `rho.negation` (reference) and `kangaroo.vow`.
- [`specs/followup-koblitz.json`](specs/followup-koblitz.json): the m = 41
  curve, arms `rho.signed_frobenius` (reference) and `rho.plain`.
- **64 fresh targets** per curve (`target_seed` 20261003, disjoint in
  derivation from round 1's 20261002), 2 rounds each (algorithm seed
  20261003), the same binary and accounting as round 1.

## Predictions

On the fresh data, with the same metric (mean `S`, 95 % two-stage
bootstrap interval, 4 000 resamples, seed 20261002):

- **F1.** `rho.plain` on m = 41: `√(π/2)` lies inside the interval, and
  `S/theory ∈ [0.85, 1.35]`.
- **F2.** `kangaroo.vow` on fp26: `2.0` lies inside the interval, and
  `S/theory ∈ [0.8, 1.5]`.
- **F3.** Every execution verifies, and `verify --replay 12` is
  identical.

A note derived now, before the data: this implementation's mean jump is
`(2^k − 1)/k` for the largest `k` not above `√r/2`, which is about
`0.67·√r/2` on fp26. The van Oorschot–Wiener cost `2(E|k − r/2|/μ + μ)`
then gives about `2.16`, not `2.0`. F2 still tests the preregistered
`2.0`. A mean well below `2.0` on 64 targets would indicate an
implementation or accounting defect.

## Interpretation, fixed in advance

- **F1 and F2 pass:** round 1's two failures are attributed to sample
  size (a band-design error for rho, an 8-target sample for the
  kangaroo), and remain recorded as failures of round 1's predictions.
  Calibration on the methods tested is accepted, with the lesson that a
  tolerance on a mean must be set from that mean's sampling spread.
- **F2 fails, with the interval excluding 2.0:** the kangaroo's cost
  differs from its analysis. Investigate the jump schedule and the
  distinguished-point handling, and keep `kangaroo.vow` out of any
  comparison until it is resolved.
- **F1 fails:** investigate the tuned walk's set-up and cycle accounting
  on large `r`.
