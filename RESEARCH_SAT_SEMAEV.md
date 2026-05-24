# SAT-encoded binary Semaev systems

**Module:** `src/cryptanalysis/semaev_sat.rs`
**Provenance:** Bridges the existing CDCL SAT solver
(`cryptanalysis::sat`) to the binary-Semaev `S_3` Weil-descent
pipeline.  Direct response to the "Hybrid SAT / algebraic" entry in
the ECDLP research-direction map.

## Pipeline

1. **Weil-descend** binary Semaev `S_3` into `n` quadratic
   polynomials over `F_2` in `2n` bit-variables (re-uses
   `ffd_harness::weil_descend_s3`).
2. **Tseitin-encode**: one auxiliary boolean `z_{i,j} = xᵢ ∧ xⱼ`
   per quadratic monomial appearing in any equation.  Three 3-literal
   clauses per AND.
3. **XOR-encode** each equation `f = c ⊕ ℓ₁ ⊕ … ⊕ ℓ_w = 0`:
   direct truth-table expansion for `w ≤ 4`, intermediate-aux chain
   for wider parities.
4. **Solve** with the existing CDCL.
5. **Decode** bit-variables → `(X_1, X_2) ∈ F_{2^n}²` and verify
   `S_3(X_1, X_2, x_3) = 0`.

## Status

- `n = 4` round-trip works in ms (`semaev_sat_round_trip_n4`).
- `n = 5` solvable but takes minutes with the textbook CDCL; marked
  `#[ignore]`, run with `--ignored`.
- `n ≥ 6` not currently tractable without XOR-clause-native
  reasoning.

## Why XOR-native matters

Our CDCL solver doesn't have XOR-clause Gaussian elimination
(CryptoMiniSat's killer feature).  For Semaev-after-Weil-descent the
system is *literally* a bunch of XORs over `F_2`, so a XOR-aware
solver would do Gaussian elimination on the linear part first and
hand the small residual to the CDCL search.  That single feature
is what would push tractability from `n = 4` to `n ≈ 12`.

## Tests

```bash
cargo test --release --lib cryptanalysis::semaev_sat
cargo test --release --lib cryptanalysis::semaev_sat -- --ignored   # for the n=5 round-trip
```

Default suite (4 tests):
- `direct_xor_emits_correct_clause_count`
- `xor_only_system_is_satisfiable`
- `semaev_sat_round_trip_n4`
- `inconsistent_equation_yields_unsat`

Marked `#[ignore]`:
- `semaev_sat_round_trip_n5`

## What's missing

- XOR-clause-native propagation (CryptoMiniSat-style Gaussian
  elimination on XOR constraints).  Without it, the SAT search has
  to *learn* the linear relations the hard way.
- Symmetry-breaking clauses (the F_2-coefficient symmetry between
  X_1 and X_2 doubles the search space unnecessarily).
- DIMACS round-trip for external solvers (kissat, CaDiCaL).  The
  `to_dimacs` function in `sat.rs` is already wired; just need a
  CLI shim.

## References

- M. Soos, K. Nohl, C. Castelluccia, *Extending SAT solvers to
  cryptographic problems*, SAT 2009 (CryptoMiniSat).
- G. Tseitin, *On the complexity of derivation in propositional
  calculus*, 1968 (AND/XOR gate-CNF encoding).
- C. Cid, S. Murphy, M. Robshaw, *Algebraic aspects of the Advanced
  Encryption Standard*, 2006 (ANF↔CNF).
