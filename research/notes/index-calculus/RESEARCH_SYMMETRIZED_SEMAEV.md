# Symmetrised summation polynomials (Faugère–Gaudry–Huot–Renault)

**Module:** `src/cryptanalysis/symmetrized_semaev.rs`
**Provenance:** Direct response to the *Symmetrised Summation
Polynomials* paper the user uploaded (§4 of the ECDLP research-
direction map: "Exactly the paper you uploaded").

## What's there

- Sparse multivariate polynomial type `MPoly` — `BTreeMap<Vec<u32>,
  FieldElement>` exponent → coefficient.
- **Elementary symmetric polynomials** `e_k(X_1, …, X_n)` and
  **power sums** `p_k(X_1, …, X_n)`.
- **Decomposition algorithm** (`decompose_symmetric`): given a
  symmetric polynomial in `n` variables, expresses it in
  `F_p[e_1, …, e_n]`.  Uses the standard leading-monomial-cancel
  recursion (Sturmfels §1.1; Macdonald's "Symmetric Functions").
- `symmetrise_s3` — applies the decomposition to the prime-field
  Semaev `S_3`, yielding the Faugère et al. compact form in 3
  elementary-symmetric variables.

## What this enables

The whole point of symmetrisation is **density reduction**:

| Polynomial | # X-monomials | # e-monomials (claim) | Ratio |
|---|---:|---:|---:|
| `S_3` | 17 | ≤ 10 | ~1.7× |
| `S_4` | ~500 | ~100 | ~5× |
| `S_5` | ~10⁵ | ~10⁴ | ~10× |
| `S_8` | ~10¹⁵ | ~10¹³ | ~100× |

For `S_8` over a degree-5 extension (the paper's headline result)
this is the difference between "memory-bound" and "computable in
hours".  This module makes the **first** of those reductions
verifiable in CI (the `s3_symmetrisation_reduces_monomial_count`
test asserts it).

## What's missing (relative to the paper)

- The actual `S_5` / `S_8` symmetrisations.  The infrastructure
  (`decompose_symmetric` + `semaev_higher::s5_in_last`) is in place;
  applying it to fully-symbolic `S_5` requires extending
  `s5_in_last` to return an `MPoly` in 5 variables instead of a
  univariate in the last.
- Newton-identity inverse change of variables (e ↔ p_k) — the building
  blocks are there (`power_sum`) but not the explicit Newton iteration.
- Block-Vieta variables: the paper uses a finer change of variables
  (splitting `e_k` by orbit) for `S_8` specifically.

## Tests

```bash
cargo test --release --lib cryptanalysis::symmetrized_semaev
```

All 5 pass:
- `elementary_symmetric_sanity` (cardinalities of e_0, e_1, e_2, e_3)
- `mpoly_arithmetic_sanity`
- `decompose_power_sum_2`  (`p_2 = e_1² − 2 e_2` round-trip)
- `symmetrise_s3_runs`
- `s3_symmetrisation_reduces_monomial_count`  (the headline claim)

## References

- J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using symmetries
  in the index calculus for elliptic curves discrete logarithm*,
  J. Cryptology 2014.
- I. Macdonald, *Symmetric Functions and Hall Polynomials*, Oxford
  1995, Ch. 1.
- B. Sturmfels, *Algorithms in Invariant Theory*, Springer 1993, §1.1.
