# Linearization reach of the bilinear S₃ core, and the gap to 2^61

**Stage diagnostic.** `S`, end-to-end cost and speedup are **unset**. No work
below `2^61` is claimed. **Class: boundary (negative result).** Nothing here
moves a ratio to rho. The note measures how far one algebraic lever reaches,
and states in bits how far short of rho it falls.

## The question

A decomposition oracle for ECC2K-130 that grows more slowly than enumeration
already exists: pairs-and-solve, and more generally resolving one summand by
linear algebra ([`RESEARCH_SEMAEV_DECOMPOSITION.md`](../notes/index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md),
[`RESEARCH_CHAIN_SPLIT_ORDER.md`](../notes/ecc2k130/RESEARCH_CHAIN_SPLIT_ORDER.md)).
To go below `2^61`, an oracle must stay cheap on systems well past `n/3`, where
linearization stops. SAT and Gröbner are measured to be exponential there
([`../frobenius_quotient_sat_20260930`](../frobenius_quotient_sat_20260930/README.md)).
This note makes both halves exact:

1. how many bits of a decomposition one trial can resolve by linearization,
   and where that stops;
2. how many bits an oracle must resolve, per trial, at every arity at `n = 131`
   to come in under rho.

## Boundaries, stated in `PROTOCOL.md` before measuring

- **Reference.** Matched rho on ECC2K-130, `2^60.8090` (RHO = 60.8090 below),
  including the `⟨−1⟩ × ⟨π⟩` speedup.
- **Floor for this family.** Per relation, the guess-and-linearize oracle costs
  about `2^{n − l − b}`. Once `Pr[decompose] ≈ 1`, the per-attempt budget under
  rho is about `2^{RHO − l}`. So the gap is `n − RHO − b_max = 70.19 − b_max`
  bits, **independent of arity**. `70.19` is the same constant as the product-law
  floor in the autoresearcher's KN-FIND-aa2efc.

## Result 1: the bilinear core has a symmetry defect of exactly `b(b+1)/2`

With `t` fixed, `S₃(X₂, X₃, t) = 0` is `F₂`-bilinear. `X₂` ranges over `V`
(`dim l`) and `X₃` over a coset `w₀ + W` with `W < V` (`dim b`). The coset must lie
in `V`, because `X₃` is a factor-base element. Linearizing gives
`N = lb + l + b` monomials in `n` bit equations. Two families of column identities
hold for every instance:

- `col(c_w) + col(d_w) + Σ_{v_i ∈ supp(w₀)} col(c_i d_w) = 0` for each basis
  vector `w` of `W` (`b` identities);
- `col(c_w d_{w′}) = col(c_{w′} d_w) = w²w′² + ww′t` for `w ≠ w′` (`C(b,2)`
  identities).

Both come from `S₃` being symmetric in `X₂, X₃` when the two share the
directions of `W`. So the rank is at most `N − b(b+1)/2`.

**Confirmed (Amendment 1, A1).** Below the corrected cliff the defect equals
`b(b+1)/2` in at least 93.5% of trials in every one of 17 cells. Above it, the
defect exceeds `b(b+1)/2` in 100% of trials in 10 cells. **The original protocol
predicted a defect of `b`, and that prediction (P1) failed.** It is kept as
failed; the correction is additive.

**Resolving the `2^{b(b+1)/2}`-point family costs polynomial time.** The
identities leave `c_w + d_w`, `c_w d_w` and the pair sums determined. So at most
two product-consistent points survive: one assignment and its `X₂ ↔ X₃` swap on
`W`. Measured: at most 2 in 99.99% of 16,190 fully enumerated consistent solves
(A2). A polynomial structured resolver agrees with full enumeration on every
canonical-defect solve, with 0 disagreements (A3). 1,227 solves with an extra,
accidental defect next to the cliff fell back to enumeration.

## Result 2: per-trial success follows `2^{l+b−n}`, and linearization past `n/3` adds ≤ 1 bit

**Reach.** Linearization is determined when `l(b+1) + b − b(b+1)/2 ≤ n`.
Hence `b_max ≥ 2` iff `l ≤ (n+1)/3`. Past that, the lever adds at most one bit
beyond guessing `X₃`. At `n = 131`, `b_max(44) = 2` and `b_max(45) = 1`. The lever
is largest at `l = b = 14` (`b_max = 14`), where the condition becomes
`(u² + 3u)/2 ≤ n` with `u ≈ √(2n)`.

**Law (confirmatory run, seed `20261001`, P3 and P4 passed).** Every success
below is verified on the curve, with 0 algebra mismatches and 0 curve rejects.
The `b` range run was set by the pre-amendment formula `⌊(n−l)/l⌋`, which agrees
with the corrected one on these cells.

| n | l | b | verified / trials | log₂ p | law `l+b−n` | deviation |
|---:|---:|---:|---:|---:|---:|---:|
| 13 | 3 | 0 | 20 / 40000 | −10.97 | −10 | −0.97 |
| 13 | 3 | 1 | 45 / 40000 | −9.80 | −9 | −0.80 |
| 13 | 3 | 2 | 112 / 40000 | −8.48 | −8 | −0.48 |
| 13 | 3 | 3 | 140 / 40000 | −8.16 | −7 | −1.16 |
| 13 | 4 | 0 | 75 / 40000 | −9.06 | −9 | −0.06 |
| 13 | 4 | 1 | 147 / 40000 | −8.09 | −8 | −0.09 |
| 13 | 4 | 2 | 262 / 40000 | −7.25 | −7 | −0.25 |
| 17 | 4 | 0 | 2 / 60000 | −14.87 | −13 | −1.87 |
| 17 | 4 | 1 | 15 / 60000 | −11.97 | −12 | +0.03 |
| 17 | 4 | 2 | 28 / 60000 | −11.07 | −11 | −0.07 |
| 17 | 4 | 3 | 35 / 60000 | −10.74 | −10 | −0.74 |
| 17 | 5 | 0 | 11 / 60000 | −12.41 | −12 | −0.41 |
| 17 | 5 | 1 | 40 / 60000 | −10.55 | −11 | +0.45 |
| 17 | 5 | 2 | 56 / 60000 | −10.07 | −10 | −0.07 |
| 19 | 5 | 0 | 1 / 60000 | −15.87 | −14 | −1.87 |
| 19 | 5 | 1 | 5 / 60000 | −13.55 | −13 | −0.55 |
| 19 | 5 | 2 | 9 / 60000 | −12.70 | −12 | −0.70 |

Slopes of `log₂ p` on `b` (cells with at least 8 successes): 0.97, 0.90, 0.61,
1.17. The frozen band was `[0.6, 1.4]`, so `(17, 4)` passed only barely. Each
further linearized bit doubles the success rate, up to `b_max`. The intercept
sits about 0.5 bits below the law, a constant from the half-trace root and sign
losses. These are toy fields, `n ≤ 19`. The law is a heuristic count with a
measured slope, not a theorem.

## Result 3: the gap at `n = 131` is at least 56 bits at every arity

`costmap_a1.py` uses the corrected `b_max`. Heuristics are stated in the script:
`Pr = min(1, X^m/(m!·2^n))`, relations = columns, and linear algebra
`= m·columns²`. The table gives the best cell per arity; the base dimension is
`c` (`l` above).

| base | m | c | Pr[dec] | attempts | LA | oracle budget | b_max | linearize / relation | gap (bits) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| subspace | 3 | 29 | 2^-46.58 | 2^75.58 | 2^59.58 | 2^-15.58 | 3 | 2^99 | 114.58 |
| subspace | 4 | 29 | 2^-19.58 | 2^48.58 | 2^60.0 | 2^11.0 | 3 | 2^99 | 88.0 |
| subspace | 5 | 28 | 1 | 2^28.0 | 2^58.32 | 2^32.53 | 3 | 2^100 | 67.47 |
| subspace | 8 | 19 | 1 | 2^19.0 | 2^41.0 | 2^41.81 | 7 | 2^105 | 63.19 |
| subspace | 12 | 14 | 1 | 2^14.0 | 2^31.58 | 2^46.81 | 14 | 2^103 | 56.19 |
| orbit union | 3 | 29 | 2^-25.48 | 2^54.48 | 2^59.58 | 2^5.52 | 3 | 2^99 | 93.48 |
| orbit union | 4 | 27 | 1 | 2^27.0 | 2^56.0 | 2^33.76 | 4 | 2^100 | 66.24 |
| orbit union | 8 | 14 | 1 | 2^14.0 | 2^31.0 | 2^46.81 | 14 | 2^103 | 56.19 |

The full table, with more arities, is in `results/costmap_a1.md`.
`results/costmap.md` is the pre-amendment table. Its floor of `60.19` is
superseded, not deleted.

## Reading

- **The goal is half met.** The slower-growing method is real and measured:
  each trial resolves `l` bits by one linear solve, plus up to `b_max` more by
  linearization at polynomial cost.
- **The rest is not in reach.** In the arity ≥ 5 cells, the oracle budget per
  attempt is `2^{32.5}` to `2^{46.8}`. To come in under `2^61`, one trial within
  that budget must resolve `n − log₂(budget) ≈ 84–99` bits. Linearization
  resolves `l + b_max = 28–31` there. The remaining 56–67 bits must
  come from something that stays cheap well past `n/3`. There, linearization is
  exhausted by the reach formula, and SAT and Gröbner solvers are measured to be
  exponential.
- **What would close the gap.** An oracle whose algebraically resolved bits per
  trial grow linearly in `n` on underdetermined descended systems (`m·l ≳ n`),
  at cost `2^{o(n)}` per trial. Nothing in this repository or in the literature
  cited in these notes provides one.
- **Not ruled out.** Methods outside the guess-and-linearize family: higher-degree
  XL on the bilinear core, structured Gröbner on bilinear systems, and
  representation-style filtering. The gap is a boundary for this family, not a
  lower bound for the ECDLP.

## Files and reproduction

| file | role |
|---|---|
| `PROTOCOL.md` | frozen before the confirmatory run; predictions P1–P4 |
| `AMENDMENT-1.md` | additive correction after P1 failed; predictions A1–A3 |
| `lr.py` | field, curve, the oracle, the rank-defect cells (frozen) |
| `amend1.py` | A1–A3: corrected cliff, family size, structured resolver |
| `costmap.py`, `costmap_a1.py` | `n = 131` cost map (pre- and post-amendment `b_max`) |
| `score.py` | scores P1–P4 from `results/results.json` |
| `test_lr.py` | curve arithmetic, `S₃` on real sums, column formula, the defect identity |
| `results/` | `results.json`, `score.json`, `amend1.json`, `amend1_score.json`, logs, cost maps |

```sh
python3 -m unittest test_lr -v
python3 lr.py results 20261001 && python3 score.py      # about 8 min, one core
python3 amend1.py results 20261002                      # about 3 min
python3 costmap_a1.py
```

Scores: P1 **failed**; P2 passed (uninformative for `b ≥ 2`, see Amendment 1);
P3 and P4 passed; A1–A3 passed. Pure Python 3, no dependencies. Single-trial
timing is not reported: every cost here is a count of trials or of bits, not
of seconds.
