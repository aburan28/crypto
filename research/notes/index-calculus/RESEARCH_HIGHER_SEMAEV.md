# Higher-order Semaev polynomials over prime fields

**Module:** `src/cryptanalysis/semaev_higher.rs`
**Demo:**  `cargo run --release --example semaev_higher_demo`

## What's there

Computes `S_3, S_4, S_5, S_6` for prime-field curves
`y² = x³ + a x + b`, returned as univariate polynomials in the *last*
variable with the first `n − 1` slots evaluated at fixed `F_p`
constants.  Built via Semaev's recurrence and an inline 8×8
Sylvester-matrix determinant for `S_6`.

## Observed degrees

| `n` | `deg(S_n)` in last variable | # coefficients |
|---:|---:|---:|
| 3 | 2 | 3 |
| 4 | 4 | 5 |
| 5 | 8 | 9 |
| 6 | 16 | 17 |

Matches the theoretical doubling progression exactly (Semaev 2004 §3
proves `deg_{X_n}(S_n) = 2^{n−2}`).

## What this unlocks

- **5- and 6-decomposition index calculus** — needed for any prime-
  field IC asymptotically better than Pollard ρ.
- **Symmetrised-Semaev pipeline** — `S_5` over the `e_1..e_5` basis
  reduces monomial count by `~5! ≈ 120`, the Faugère–Gaudry–Huot–
  Renault speed-up.
- **Cross-checking the FFD harness** — applying the FFD measurement
  to `S_5` rather than `S_3` tests the FFD conjecture in a more
  demanding regime.

## Limits

- All `S_n` are returned with the first `n−1` slots already evaluated
  at `F_p` constants.  Full symbolic `S_5` / `S_6` would need a
  multivariate polynomial container; the [`MPoly`] type in
  `symmetrized_semaev.rs` can be used for that route.
- No specialised modular reduction or Schönhage-Strassen multiplication;
  arithmetic is naïve.  Adequate for `p ≤ 2³²` toy curves; for
  cryptographic `p` a single `S_6` resultant takes minutes.

## Tests

```bash
cargo test --release --lib cryptanalysis::semaev_higher
```

All 6 pass:
- `s3_in_last_matches_definition`
- `s3_vanishes_on_known_collinear_triple`
- `s4_in_last_is_degree_four`
- `s5_in_last_is_degree_eight`
- `s5_in_last_is_symmetric_in_first_four_args`
- `s6_in_last_runs`
