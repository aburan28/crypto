# Item 6: exact decomposition census against Trimoska's thesis Table 7.7

**Instrument.** [`src/main.rs`](src/main.rs) (standalone Rust crate, no
dependencies).  **Frozen output.** [`results/census_n24_l4_5.txt`](results/census_n24_l4_5.txt),
[`results/census_n24_l6_9.txt`](results/census_n24_l6_9.txt).
**Status.** Exploratory: run before a prediction was registered, so it
classes as **accounting** and nothing here is a result against a boundary.
A bug in the first build (the lifting solver read pivot transforms before
later columns were cleared; `assert!(on_curve)` caught it at `x = 1`) was
fixed before any number below was written; the buggy output was discarded.

## Question

The survey's §6 item 6: thesis Table 7.7 (`n = 24`, `m = 3`, `ℓ = 4…8`)
reports a decomposition probability `P` about `2×` the heuristic
`P_approx = 2^{mℓ}/(m!·2^n)`.  This repository's yield law
`E[#decompositions] = C(|F|, m)/#E` was confirmed on twelve toy cells.  Which
is right, and what is the `2×`?

## Method

Curve `E₁ : y² + xy = x³ + x² + 1` over `F_{2^24}` (the thesis' curve
equation; the repository's first irreducible of degree 24, `0x100001b`, since
the thesis generator's modulus is not vendored), factor base
`F_ℓ = {P : x(P) ∈ V_ℓ}`, `V_ℓ = ⟨1, z, …, z^{ℓ−1}⟩`.  Every unordered triple
of factor-base points is summed and `x(R)` counted, so for every abscissa of
the field the number of x-triples the x-only `S₄` system admits is exact.
The same census is run on the quadratic twist (`a = 1 + δ`, `Tr δ = 1`),
because an x-only instance whose target abscissa does not lift to `E₁` is
satisfiable exactly when the twist point decomposes over the twist's base.
Checked: `#E₁ + #E₁' = 2(2^24 + 1)`; every lifted point and every sum
re-checked on the curve; the census total equals the triple count.

## Table

`P_exact` is the fraction of abscissae that lift to the curve and have at
least one x-triple.  `P_thesis` is Table 7.7's column converted to a
probability: the column is `#sat / #unsat` (`256/688 = 0.37209` reproduces
the `ℓ = 8` entry exactly), so `P = #sat / (#sat + #unsat)`.  `λ_model =
C(|F|, 3)/#E` is the repository's yield law with the exact `|F|`.

| ℓ | `|F|` curve / twist | `P_exact` curve | `P_exact` twist | `P_approx` | `1 − e^{−λ_model}` curve | `P_exact / P_approx` | `P_thesis` (Table 7.7) | `P_thesis / P_exact` |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| 4 | 15 / 17 | 2.26e-5 | 3.45e-5 | 4.07e-5 | 2.71e-5 | 0.56 | 1.07e-4 | 4.7 |
| 5 | 29 / 35 | 1.94e-4 | 3.55e-4 | 3.26e-4 | 2.18e-4 | 0.60 | 6.46e-4 | 3.3 |
| 6 | 59 / 69 | 1.83e-3 | 2.98e-3 | 2.60e-3 | 1.94e-3 | 0.70 | 5.54e-3 | 3.0 |
| 7 | 123 / 133 | 1.73e-2 | 2.19e-2 | 2.08e-2 | 1.79e-2 | 0.83 | 4.10e-2 | 2.4 |
| 8 | 255 / 257 | 0.1375 | 0.1405 | 0.1667 | 0.1502 | 0.83 | 0.2712 | 1.97 |
| 9 | 507 / 517 | 0.4611 | — | 1.333 | 0.7238 | 0.35 | — | — |

Per-target x-triple counts (curve, `ℓ = 9`, `λ_x = 2.57`): `0: 53.9 %`,
`2: 9.9 %`, `4: 12.7 %`, `6: 10.8 %`, `8: 6.9 %`, `≥ 9: 5.9 %`.  Counts are
always even because `{−P₁, −P₂, −P₃}` sums to `−R`, which has the same
abscissa.

## Reading

1. **The repository's yield law holds for the mean.**  `λ_measured / 2`
   (x-triples per lifting abscissa, halved for the `±R` pairing) equals
   `C(|F|, 3)/#E` to three digits at every `ℓ` (`1.283` vs `1.287` at
   `ℓ = 9`).
2. **The exact satisfiable fraction is below the heuristic, not above it.**
   `P_exact / P_approx = 0.56 → 0.83` on the curve, `0.85 → 1.1` on the
   twist, for two reasons that are both in the yield law and not in the
   heuristic: `|F| ≈ 0.9 · 2^ℓ` (not every `x ∈ V` lifts; the two curves
   share the `2^ℓ` abscissae), and `P ≤ λ`.
3. **At `λ ≳ 1` the hit probability is well below Poisson.**  At `ℓ = 9`,
   `P_exact = 0.46` against `1 − e^{−1.29} = 0.72`: the x-triple counts are
   strongly over-dispersed (54 % of targets have none, 6 % have nine or
   more).  The mean is right; the fraction of decomposable targets at
   `λ ≈ 1` is about `0.64×` the Poisson figure.  This is an accounting
   correction to the "targets per relation" box of the product law in the
   regime `C(|F|, m) ≈ #E`, which is exactly the regime an attack would
   choose; it moves a constant, not an exponent.
4. **The thesis' `2×` is not reproduced.**  Converted to a probability, the
   thesis column sits `2.0×` (`ℓ = 8`) to `4.7×` (`ℓ = 4`) above the exact
   census on the same curve equation, standard subspace and either target
   stratum.  The yield law cannot be the source.  What remains is the
   thesis experiment's instance generation (the target sampling and the
   modulus, neither vendored here) or a model that admits solutions the
   point question does not; deciding that needs the upstream generator,
   which this environment cannot fetch (`github.com` → 403).  The
   repository's own audit of the upstream `-U` corpus (1 satisfiable of 30
   random targets where the heuristic predicts ≈ 6) points the same way as
   this census, not as Table 7.7.

## Reproduce

```sh
cargo build --release
./target/release/yield_census 24 4,5,6,7,8,9     # ≈ 10 min; ℓ = 9 dominates
```
