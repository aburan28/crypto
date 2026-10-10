---

## 7. Experimental research program

Common protocol. Every cell records: field (q, n), curve (a, b, j), N and its
factorisation, r, h, D_K, f and v_ℓ(f) for each ℓ | f, the isogeny (ℓ, kernel
rationality |M^π|, cocycle class distribution over a sample of targets),
fiber statistics (histogram of rational fiber sizes, π-orbit lengths over a
sample), system size (variables, equations, monomials, density), first-fall
degree and degree of regularity observed, F4/SAT wall and peak memory,
factor-base build time, verified relation yield per target, duplicate and
failed attempts, linear-algebra cost model (|F|² · weight), and the projected
end-to-end cost (collection + algebra + descent). Primary metric:
**amortised cost per new independent verified relation** (independent = after
cofactor projection and τ-/K-orbit identification, measured by incremental
rank over Z/r). Secondary: projected end-to-end cost. Baseline for each cell:
the fastest existing arm in this repository on the *same* curve (binary:
`icv1-*` F4 cells; prime: enumeration PDP of PR #1560; rho with
automorphisms where applicable). Each arm runs the same target list (seeded)
and the same wall budget; three repetitions; report medians and the paired
bootstrap interval the ecbench harness already produces.

### A. Binary Koblitz curves, m ∈ {31, 51, 53, 83}

A1 (H1 + H4, the main experiment). K₀/F_{2^m}, K = ⟨(0, 1)⟩ (rational, F₂-
defined), E₂ = K₀/K by the char-2 Vélu formula, u = x + 1/x (b = 1).
Arms: (i) plain S₃^{K₀}(x) with τ-orbit factor base (ledger baseline);
(ii) S₃^{E₂}(u) with the ⟨τ⟩ × K orbit table of §4.5; (iii) as (ii) with the
FHJRV elementary-symmetric variables. V of F₂-dimension ⌈m/3⌉ (two-summand
cells as in the ledger), identical target lists. Prediction: identical F4
wall per target (Prop. 3.2), identical yield per target, ×2 fewer relations
needed ⇒ primary metric improves by ≤ 2 for (ii) and by 2·(1+ε) for (iii)
with ε the symmetrisation sparsity gain. Kill criterion: < 1.5× at m = 31
and m = 51.

A2 (H9 negative control). Enumerate Φ₂, Φ₃ roots of j(K₀) over F₂ and over
F_{2^m} for m = 31, 51; list the F₂-rational members of the class (expected:
K₀ and its twist only); for one non-F₂-rational neighbour measure the factor
base size without τ-compression. Prediction: loss factor m/|M^π| ≥ 7.

A3 (toward ECC2K-130, m = 83 and 131). Cost model only: take A1's measured
per-target F4 cost at m = 51, 53, extrapolate with the ledger's degree-4
Macaulay growth, and fold in the ×2 relation saving (×4 for K₀, a = 0, whose cofactor is 4
and whose F₂-rational point of order 4 gives K = ⟨T₄⟩, k = 4; ECC2K-130 is
of this type).
Report the projected end-to-end ratio against rho on ECC2K-130
(the published 2^{60.9}-ish iteration estimate).

### B. Ordinary prime-field curves with controlled CM data

B1 (H3). `csd/volcano.py`: D_K = −3, ℓ ∈ {2, 3, 5}, height 2, bits ∈ {28, 32,
36}; counted rho on level 0 (j = 0, automorphism arm) vs level 2 (plain
arm) on the same N, r. Prediction: ratio √3 ± 5 %; plain arms equal.

B2 (Prop. 3.3). Same families; measure the cost of computing the ascending
ℓ-isogeny from the floor (kernel polynomial via ψ_ℓ factorisation) for
ℓ ∈ {2, 3, 5, 7, 11, 13} and fit the exponent in ℓ; compare evaluation by
Vélu vs √élu. This calibrates the "reachability" constant in H2.

B3 (H5 negative). On a 24-bit curve with a rational 3-torsion point:
small-x factor base of size B inside the index-3 image vs unrestricted,
same B; yield ratio must be 3^{−(m−1)} after target correction.

B4 (Prop. 3.2). Horizontal 2-, 3-, 5-isogenous neighbours of a 32-bit
curve with h(D_K) ≥ 20: enumeration-PDP relation yield and timing per
target with B = 120 (PR #1560 harness); prediction: identical to within
|F| variation (< 5 %).

### C. Prime-field extension curves, n ∈ {2, 3, 5}

C1 (H1 over F_{p²}). p ≈ 2^{20}, E with N even (rational T₂) and a
2-isogenous E₂; Gaudry factor base {x ∈ F_p} on E vs {u ∈ F_p} on E₂ with
the K-coset table; S₃ systems after Weil descent (2 variables each). Primary
metric and degree of regularity; prediction ×2 on relations needed, same
F4 cost per target.

C2 (H8). n = 3, p ≈ 2^{14}: trace-zero index calculus (Gorla–Massierer)
on E and on E′ 2-isogenous; costs equal within |F| variation.

C3 (Prop. 2.3 statistics at scale). n = 2, 3, 5, p ≈ 2^{16}: for 2- and
3-isogenies with M^π ∈ {1, 2, 3}, histogram rational fiber sizes and orbit
lengths over 10⁴ targets; must match {0, |M^π|} and the affine-map orbit
law exactly (a deviation would falsify Prop. 2.3 and is the cheapest
possible sanity check of the whole framework).

### Reproducible pseudocode (SageMath; the pure-Python probe is the fallback)

```python
# cell(E1, K, V, targets): one arm of A1/C1
E2 = E1.isogeny(K.gens()[0]); phi = E1.isogeny(K.gens()[0])
F2 = {Q for Q in E2.points() if Q[0] in V}            # or lifted from x in V
pre = {}                                               # x(phi(P)) -> one preimage
for P in E1.points():
    if phi(P)[0] in V: pre.setdefault(phi(P)[0], P)
rows = []
for R in targets:
    Rq = phi(R)
    sols = solve_S3_weil_descent(E2, V, Rq)            # F4 / SAT, symmetrised vars
    for (u1, u2) in sols:                              # m = 2 shown
        Q1, Q2 = lift_x(E2, u1), lift_x(E2, u2)
        for s1, s2 in signs:
            if s1*Q1 + s2*Q2 == Rq:
                P1, P2 = pre[u1], pre[u2]
                assert h*(s1*P1 + s2*P2) == h*R        # Prop. 4.1
                rows.append(((u1, s1), (u2, s2), R))
rank = incremental_rank_mod_r(rows)                    # independent relations
report(cost_per_independent_relation = wall / rank)
```
