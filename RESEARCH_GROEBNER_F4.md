# Minimal Buchberger / matrix-F4 Gröbner-basis solver

**Module:** `src/cryptanalysis/groebner_f4.rs`

## What's there

- Multivariate polynomial type ([`MPoly`] from `symmetrized_semaev.rs`).
- Three monomial orderings: `Lex`, `Grlex`, `Grevlex`.
- S-polynomial computation, multivariate polynomial reduction.
- **Buchberger's algorithm** with Gebauer–Möller coprime-LM pruning.
- **`reduce_basis`** to canonicalise to the reduced Gröbner form.

## What this unlocks

- The "Gröbner solve" step that the PKM, FFD, symmetrised-Semaev, and
  Diem-descent modules all eventually need — without an external CAS
  dependency (Magma, Singular, msolve).
- The minimum-viable backend for the symmetrised `S_5` / `S_6`
  pipeline.

## Performance envelope

- ~6 variables, degree ≤ 4, ≤ 5 polynomials → seconds.
- 8+ variables: starts to hit the 5,000-step Buchberger safety brake.
- True F4 (signature reduction, sparse Macaulay matrices) would push
  this to ~12-15 vars; not yet implemented.

## Tests

```bash
cargo test --release --lib cryptanalysis::groebner_f4
```

All 6 pass:
- `buchberger_circle_intersect_line_over_f7` (the canonical 2-var
  textbook example)
- `buchberger_solves_linear_system`
- `buchberger_constant_ideal`
- `s_polynomial_of_coprime_polys_pre_reduction`
- `reduction_of_x_squared_plus_y_by_x`
- `grevlex_ordering_within_a_degree`

## What's missing (relative to a research-grade F4)

- F5 signature-based reduction (avoids redundant zero-reductions).
- FGLM monomial-order conversion (lex → grevlex switch).
- Sparse linear algebra (Lanczos / Wiedemann).
- Hilbert-driven D selection (currently we just iterate until
  Buchberger condition holds).
- Field-specific specialisations (e.g. F_2 with XOR-only arithmetic).

## References

- B. Buchberger, *An algorithm for finding the basis elements of the
  residue class ring of a zero dimensional polynomial ideal*, 1965.
- J.-C. Faugère, *A new efficient algorithm for computing Gröbner
  bases (F4)*, J. Pure Appl. Algebra 1999.
- R. Gebauer, H. M. Möller, *On an installation of Buchberger's
  algorithm*, J. Symb. Comp. 1988.
