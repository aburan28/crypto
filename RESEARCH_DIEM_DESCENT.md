# Diem-style index calculus on `E / F_{p^k}` — toy implementation

**Module:** `src/cryptanalysis/diem_descent.rs`
**Demo:**  `cargo run --release --example diem_descent_demo`

## What's there

Demonstrates the Diem 2011 framework on the smallest meaningful
example, `E / F_{25}`:

1. Hand-rolled `F_{25} = F_5[θ]/(θ² − 2)` arithmetic.
2. Curve `y² = x³ + (1+θ) x + 2`.
3. Factor base `FB = { P : x(P) ∈ F_5 }`.
4. Brute-force 2-decomposition search `R = ±P_1 ± P_2` over `FB²`.

## Headline result on the toy

On `E / F_{25}`:
- Factor base = 4 points
- 18 / 25 random targets `R = a G + b Q` (with `a, b ∈ [1, 5]`)
  admitted a 2-decomposition over `FB`.

This is the constructive demonstration that the *subfield-x-coord*
factor base captures a large fraction of group elements — the
combinatorial ingredient Diem 2011 turns into an asymptotic
sub-exponential algorithm.

## Why this matters for the "Holy Grail"

The screenshot's §10 ("Index calculus for prime fields, the Holy
Grail") asks: *can the descent be made to work over a prime field
`F_p` rather than an extension `F_{p^k}`?*

The Diem framework requires a non-trivial subfield, which prime
fields don't have.  But many proposed prime-curve attacks attempt to
**simulate** a subfield via:

- hidden isogenies (→ `ec_trapdoor.rs` already audits for this);
- generalised-Mersenne prime structure (→ `pkm_criterion.rs`,
  `p256_structural.rs`);
- char-of-discriminant trickery (→ `fght_snfs.rs`).

Until one of those routes lands a real attack, the Diem framework is
*only* applicable to `F_{p^k}` curves.  This module is the working
existence-proof.

## What's missing (relative to real Diem)

- **Brute-force decomposition.**  We sweep `FB × FB` exhaustively
  instead of solving the Semaev polynomial system.  Real Diem uses
  the `S_3 = 0` polynomial system with `(x_1, x_2)` symbolic over
  `F_p`, splits by F_p-basis, and solves the resulting `k = 2`
  equations.  Wire `groebner_f4` here for the next step.
- **No `(p, k) = (5, 3)` extension.**  The factor base for `k = 3`
  would be `S_3`-decompositions giving 1 equation in 2 unknowns over
  `F_5`, which has many solutions — and that's exactly Diem's win
  (asymptotic cost `~p^{2/3}` rather than `~p`).  Pending a
  performant Gröbner solver.
- **No linear algebra mod `n`.**  We *find* the decompositions but
  do not currently chain them into a Gaussian-elimination DLP
  recovery.  The infrastructure to do so already exists in
  `ec_index_calculus::gaussian_eliminate_mod_n`.

## Tests

```bash
cargo test --release --lib cryptanalysis::diem_descent
```

All 5 pass:
- `fpk_arithmetic_basic`         (θ² ≡ 2, (1+θ)² check)
- `fpk_inverse_round_trip`
- `factor_base_has_some_points`
- `curve_arithmetic_consistency` (2G + (−G) = G)
- `end_to_end_decomposition`     (R = G + H reverse-engineered)

## References

- C. Diem, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011).
- J.-C. Faugère, A. Joux, V. Vitse, *Improving the complexity of
  index calculus algorithms in elliptic curves over small extension
  fields*, EUROCRYPT 2012.
- A. Amadori, F. Pintore, M. Sala, *On the discrete logarithm
  problem for prime-field elliptic curves*, FFA 51 (2018).
