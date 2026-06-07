# initial_minors — testing the Mahalanobis–Abdullah–Mallick conjecture

Pure-Python (no numpy/sympy/Sage) experiments for §4 of
[`RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md`](../../RESEARCH_INDEX_CALCULUS_FACTOR_BASE.md).

The "minors" method solves ECDLP by searching for a vanishing minor of the
matrix `M[i][j] = φ_j(R_i)`, `R_i = a_iP + b_iQ`, where `φ_j` is the
pole-order-ordered Riemann–Roch basis `[1, x, y, x², xy, …]`. A `k×k` minor
vanishes **iff** the `k` rows (points) sum to `O`, which yields a relation
`Σ(a_i + m b_i) ≡ 0 (mod n)` and hence the discrete log `m`.

| file | what it does |
|------|--------------|
| `ecmin.py` | toy EC arithmetic, RR basis, mod-`p` determinant / LU, leading-minor pivot scan |
| `verify_correspondence.py` | confirms `sum==O ⇔ minor singular` (299/300) and recovers a real discrete log |
| `scaling.py` | measures cost-to-first-relation vs `n`: **`Θ(n)`** (0.98–1.07·n over a 245× range) |
| `crossover.py` | shows `√n`, `L[1/2]`, `L[1/3]` are numerically indistinguishable at `n ≤ 2^50` |

```bash
python3 verify_correspondence.py   # instant
python3 scaling.py                 # a few minutes (point-counts toy curves)
python3 crossover.py               # instant
```

**Finding.** The natural minor search has relation density exactly `1/n`, so it
costs `Θ(n)` — worse than Pollard rho's `Θ(√n)`. The conjecture's claimed
subexponentiality must come entirely from the "initial minors" sub-family
finding a vanishing minor far faster than the density allows; the published
`2^50` ceiling is rho-feasible (`2^25` work) and below the regime where
subexponential behaviour is distinguishable from `√n`. Reframed falsifiable
test: does `E[initial-minors-to-relation]` drop below `√n`? Measured `Θ(n)` for
the leading-principal family says no.
