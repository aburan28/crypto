# Amendment 2: an ECM-completed E1 instrument (E1b) and an unbiased E4 height (E4r)

**Status: frozen 2026-10-12, before any v3 run.** It is dated after two
dry runs on scratch output (disclosed below). `PROTOCOL.md` and
`AMENDMENT_1.md` are not edited. This amendment adds two experiments. It
reports them as **v3** alongside v1 and v2 and does not replace them.

## Why

- **E1b.** E1's M1 is defined on what budgeted rho returns. A reader can
  fairly ask whether P-256 `n − 1`'s rank (1992/2000 in v2) is an artefact
  of that incompleteness. Analytically, Brent rho with 2^18 iterations per
  constant misses a given prime below 2^32 with probability of about
  `exp(−2^36/2^33) ≈ 3·10⁻⁴` per constant, so M1 should already be
  essentially complete. E1b tests that with a second factoring stage.
- **E4r.** v1/v2's E4 instrument has a floor. Random pairs pile up at
  `|cos|` = 0.0013, and exactly dependent pairs read 1.0011–1.0025.
  `height_sequence` takes the naive height on the denominator-cleared model
  `(u²x, u³y)`, where `u = lcm` of all denominators including the point's
  own. That inflates `h_x` by about `2 log u`, and 4⁵ doublings do not
  remove it. E4r takes the naive height of the reduced rational `x` on the
  original model. The canonical limit is unchanged, the bias is removed,
  and the cost is the same.

## E1b

- **Targets.** The 9 values of E1 (`n − 1`, `n + 1` and twist order for
  three curves).
- **Null.** The first **N = 800** Hasse-interval primes of v2's E1 stream
  (seed `20261011 + curve index`, drawn sequentially), so the same inputs
  as v2's first 800 controls.
- **Instrument (identical on every arm).** `arith::factor(x, 2^18)` as in
  E1. Then every composite it leaves gets ECM with `B1 = 2000` and 20
  curves (`cryptanalysis::ecm::ecm_factor`, seeds fixed per value),
  recursively on the pieces until each is prime or ECM fails.
- **Metric.** M1 as in E1, on the ECM-completed factor list. Rows also
  record the rho-only M1 on the same input, so completeness is measured
  directly. M2 is recorded but not tested.
- **Positive control.** 50 values per curve of the form `q32 · q100 · r`.
  `q32` is a random prime in `[2^31, 2^32)`, `q100` a random 100-bit prime
  and `r` a random prime bringing the value to the size of `n` (seed
  `20261011 + 400 + curve index`). **Pass:** M1 equals `log2 q32` exactly
  on at least 49/50 per curve.
- **Decision rule.** Two-sided add-one rank as in E1. The family is 9
  targets × M1, so α = 0.05/9 = 0.0056. Power check: `2/(N + 1)` = 0.0025
  < 0.0056, and `N + 1 ≥ 4m/α = 720` holds. Verdict wording as in
  `PROTOCOL.md`.
- **What counts as "robust".** P-256 `n − 1`'s v2 verdict ("within the
  control distribution") is robust to instrument completeness if E1b
  reaches the same verdict, and if fewer than 1% of control values change
  M1 under ECM.

## E4r

- **Targets, null, seeds, positive controls.** As in v2's E4: the same 500
  random pairs per curve (seed `20261011 + 300 + curve index`), the target
  `(G̃, R̃)` with `R = [0x1234567890abcdef]G`, and the dependent pairs
  `(G̃, [2]G̃)` and `(G̃, [3]G̃)` on `E_1`.
- **Instrument.** `height_sequence_on(…, 5, rational = true)`: the same
  Jacobian doubling, with `h_x` taken of `X/(Z²u²)` in lowest terms. The
  metrics are `|cos|` and det/(ĥĥ) at k = 5, as before. The k = 4 values
  are recorded for the accuracy check.
- **Accuracy gates (pre-registered; both must pass, else no orthogonality
  verdict).**
  1. **Dependent pairs.** On every curve, abs(`|cos|` − 1) ≤ 1e-4 and
     abs(det/(ĥĥ)) ≤ 1e-4 at k = 5.
  2. **Parallelogram law.** On the first 50 control pairs per curve,
     `max |ĥ(P+Q) + ĥ(P−Q) − 2ĥ(P) − 2ĥ(Q)|` ≤ 1e-2 nats.
- **Decision rule.** Rank as in E4. The family is 3 curves × 2 metrics, so
  α = 0.05/6. N = 500 passes the power check (`N + 1 ≥ 480`).
- **What this licenses.** If both gates pass, the measured `|cos|` of
  `(G̃, R̃)` and its rank among random pairs replace the withdrawn
  "cos < 0.002" sentence. If either gate fails, the note keeps
  "orthogonality unresolved".

## Pre-flight runtime estimate (from the dry runs, 4 vCPU)

| family | per-sample cost | fixed cost | estimate at protocol N |
|:--|:--|:--|:--|
| E1b | ≈ 2.3 s CPU per value (one ECM pass of 20 curves on the composite rho leaves) | negligible | 3 × (2400 + 3 + 50) × 2.3 s / 4 ≈ 70 min |
| E4r | ≈ 13 / 18 / 23 s CPU per pair (P-192 / P-224 / P-256); parallelogram adds a third for 50 pairs per curve | dependent-pair heights (`[3]G̃`, `[4]G̃`), single-threaded, ≈ 5 / 7 / 9 min | ≈ 2.5 h |
| total | | | ≈ 3.6 h |

The actual run time will be recorded next to this estimate.

## Disclosure: what the dry runs showed before this amendment

Both dry runs used N = 40 (E1b) and N = 10 (E4r) with 5 E1b positives per
curve.

- **E1b.**
  - No control value's M1 changed under ECM (0/120 per curve).
  - The targets' M1 were unchanged: P-256 `n − 1` stays 125.0.
  - P-256 `n − 1` ranked 40/40 among the 40 dry-run controls.
  - 15/15 planted 32-bit primes were recovered.
  - 20–27 composites per curve stayed unsplit after ECM. Those are large
    cofactors, irrelevant to M1.
- **E4r.**
  - Dependent pairs gave `|cos|` = 1 and det = 0 to about 1e-16, and the
    parallelogram law held to about 1e-13 on most pairs.
  - The canonical pairings of random two-point lifts came out as **exact
    small multiples of `log ℓ` or exactly 0**: ½ ln 3, ½ ln 5, ½ ln 7,
    ½ ln 11, and about (1/6) ln 2 with a residual of 1.4e-3 nats.
  - Target `|cos|`: 5.9e-4 (P-192), 0 (P-224), 4.4e-4 (P-256).
  - So `|cos|` is a discrete-valued statistic with many ties. The rank
    rule handles ties conservatively.

None of this changes the design above, which v2 fixed except for the
instrument. It is recorded so readers can discount v3 for having been
specified after seeing it.
