# Research Note: Point-Decomposition Speed-ups and the Isogenous-Curve Hypothesis

**Date.** 2026-10-08.  **Status.** Measured on public synthetic small-field instances.
**Goal addressed.** (1) a drastic speed-up for the point-decomposition problem
(PDP) in Semaev-polynomial index calculus; (2) whether isogenous curves give
better relation yield or lower Gröbner degrees.
**Code.** `examples/isogeny_pdp_probe.rs`, `examples/macaulay_syzygy_probe.rs`;
raw outputs in `research/pdp_isogeny_probe_20261008/`. Prime-field PDP numbers
are from `RESEARCH_PRIME_FAST_ECDLP.md` (PR #1560).

---

## 1. PDP speed-ups: what is in hand, by regime

### 1.1 Prime fields — enumeration PDP (landed in PR #1560)

Over `F_p` the summation-polynomial roots are exactly the x-coordinates of
`R ∓ F_i` (for `S₃`) and of `R ∓ F_i` matched against `x(F_j ± F_k)` (for
`S₄`), so the PDP is a batched subtraction plus a lookup rather than a root
finding. Measured on the a = −3 ladder rungs, same relations as the Semaev
sweeps:

| PDP | before | after | factor |
|---|---:|---:|---:|
| `S₃` sweep per target, 20-bit rung, B = 120 | ≈ 73 ms (352 s / 4788 targets) | ≈ 104 µs | ≈ 700× |
| `S₄` sweep per relation, 16-bit rung, B = 120 | 42 s (brute-force roots) | 5 µs | ≈ 10⁷× |
| `S₄` end to end at 32 bits | not reachable | 20.5 s | — |

This is drastic in the per-operation sense, and it is also the end of the
road for the enumeration PDP: §2.4 of `RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`
shows that any enumeration/table PDP stays on the generic `S·T² ≳ p` line
whatever its constant.

### 1.2 Binary fields — algebraic PDP (F4 on the Weil-descended system)

The ledger's binary cells solve the PDP with matrix-F4 on the
Weil-restricted Semaev system; the cost is dense `F₂` elimination of the
Macaulay matrix at the first useful degree. The one "drastic" candidate in
`docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md` (G1, Joux–Vitse F4 trace
replay, 10–100× claimed in the literature) can only save rows that reduce to
zero, so its ceiling is `rows / rank`. Measured at the landed cell shape
(`icv1-f2m31-tm90707-c95f16f5`, 16-dimensional subspace, two summands, three random
targets; `macaulay_syzygy_probe`):

| degree | rows | cols | rank | rows / rank | syzygies | stable across targets |
|---:|---:|---:|---:|---:|---:|---|
| 2 | 31 | 289 | 31 | 1.000 | 0 | yes |
| 3 | 1 023 | 4 369 | 1 021–1 022 | 1.001 | 1–2 | yes |
| 4 | 16 399 | 37 809 | 15 723 | 1.043 | 676 | yes (676 every time) |

**Only 4.3 % of the degree-4 rows are redundant**, identically on every
target. Trace replay would save at most that — the matrix is almost all
useful rows. The three-summand system at this shape exceeds the matrix
limits before degree 3. So G1 is closed as "not drastic here": the cost is
the elimination of a 16k × 38k matrix, not reductions to zero. What remains
on the binary side is structure that shrinks the matrix itself
(Frobenius-equivariant block-diagonalisation, G2 in the same note), which is
a different and larger project.

## 2. The isogenous-curve hypothesis

### 2.1 Prediction

Binary `S₃(x₁, x₂, x₃) = (x₁x₂ + x₁x₃ + x₂x₃)² + x₁x₂x₃ + b` contains the
curve coefficient `b` only as a constant, and `a` not at all. The
Weil-descended system therefore differs between isogenous curves only in its
constant terms, so the Macaulay rank profile, the first-fall degree and the
F4 behaviour should be identical across an isogeny class. Relation yield can
differ only through `|F|` (the on-curve fraction of the subspace, which
depends on `b`) and through the loss of Frobenius symmetry when the
neighbour is not defined over `F₂`.

### 2.2 Measurement

**Koblitz `K_0` base, n = 17, V = ⟨1, z, …, z^7⟩ (2^8 x-values), 2-decomposition.**

| curve | a | b | Frobenius-symmetric | \|F\| | 2-sum yield | FFD | Macaulay rank (deg 2/3/4) vs generic | F4 ms decomposable | F4 ms random | F4 found | mean reductions |
|---|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| K_0 (base) | 0x0 | 0x1 | true | 271 | 0.1218 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 42.6 | 7.3 | 6/6 | 11.9 |
| 2-isogenous #0 | 0x0 | 0x1 | true | 271 | 0.1218 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 30.9 | 9.2 | 6/6 | 11.9 |
| 2-isogenous #1 | 0x1 | 0x1 | true | 241 | 0.0993 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 29.5 | 11.8 | 6/6 | 13.1 |

**Koblitz `K_0` base, n = 19, V = ⟨1, z, …, z^8⟩ (2^9 x-values), 2-decomposition.**

| curve | a | b | Frobenius-symmetric | \|F\| | 2-sum yield | FFD | Macaulay rank (deg 2/3/4) vs generic | F4 ms decomposable | F4 ms random | F4 found | mean reductions |
|---|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| K_0 (base) | 0x0 | 0x1 | true | 527 | 0.1161 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 23.5 | 24.3 | 6/6 | 20.1 |
| 2-isogenous #0 | 0x0 | 0x1 | true | 527 | 0.1161 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 16.8 | 16.9 | 6/6 | 19.6 |
| 2-isogenous #1 | 0x1 | 0x1 | true | 497 | 0.1052 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 24.3 | 22.9 | 6/6 | 21.8 |

**random-`b` base (seed 7), n = 17, V = ⟨1, z, …, z^7⟩ (2^8 x-values), 2-decomposition.**

| curve | a | b | Frobenius-symmetric | \|F\| | 2-sum yield | FFD | Macaulay rank (deg 2/3/4) vs generic | F4 ms decomposable | F4 ms random | F4 found | mean reductions |
|---|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| random-b base | 0x0 | 0x1d43 | false | 251 | 0.1066 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 15.4 | 19.5 | 6/6 | 11.5 |
| 2-isogenous #0 | 0x0 | 0x42fd | false | 235 | 0.0958 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 26.7 | 9.3 | 6/6 | 13.3 |
| 2-twist of neighbour #1 | 0x1 | 0x42fd | false | 277 | 0.1273 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 35.0 | 11.1 | 5/6 | 10.7 |
| 2-isogenous #2 | 0x0 | 0x115ed | false | 269 | 0.1206 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 22.5 | 21.0 | 6/6 | 13.4 |
| 2-twist of neighbour #3 | 0x1 | 0x115ed | false | 243 | 0.1018 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 19.5 | 22.3 | 6/6 | 14.3 |
| 3-isogenous #0 | 0x0 | 0x3e41 | false | 245 | 0.1030 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 8.1 | 22.9 | 6/6 | 11.8 |
| 3-twist of neighbour #1 | 0x1 | 0x3e41 | false | 267 | 0.1196 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 50.7 | 15.7 | 6/6 | 12.5 |
| 3-isogenous #2 | 0x0 | 0xea39 | false | 273 | 0.1242 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 11.3 | 35.3 | 6/6 | 12.8 |
| 3-twist of neighbour #3 | 0x1 | 0xea39 | false | 239 | 0.0985 | 3 | d2: 17/17, d3: 594/595, d4: 9945/10132 | 39.1 | 16.9 | 5/6 | 9.8 |

**random-`b` base (seed 7), n = 19, V = ⟨1, z, …, z^8⟩ (2^9 x-values), 2-decomposition.**

| curve | a | b | Frobenius-symmetric | \|F\| | 2-sum yield | FFD | Macaulay rank (deg 2/3/4) vs generic | F4 ms decomposable | F4 ms random | F4 found | mean reductions |
|---|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| random-b base | 0x0 | 0x39445 | false | 499 | 0.1055 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 15.5 | 18.3 | 6/6 | 19.0 |
| 2-isogenous #0 | 0x0 | 0x2046b | false | 529 | 0.1177 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 31.2 | 26.1 | 6/6 | 17.3 |
| 2-twist of neighbour #1 | 0x1 | 0x2046b | false | 495 | 0.1046 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 22.9 | 18.4 | 6/6 | 18.0 |
| 2-isogenous #2 | 0x0 | 0x64cf1 | false | 515 | 0.1118 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 28.3 | 24.1 | 6/6 | 21.2 |
| 2-twist of neighbour #3 | 0x1 | 0x64cf1 | false | 509 | 0.1096 | 3 | d2: 19/19, d3: 740/741, d4: 13870/14098 | 16.8 | 15.9 | 6/6 | 19.2 |


### 2.3 Reading

- **Gröbner degrees do not move.** Every curve in every class — the base,
  each distinct-`j` 2- and 3-isogenous neighbour, and each neighbour's
  twist — has the *identical* Macaulay rank at degrees 2, 3 and 4 and the
  same first-fall degree (3). This is the §2.1 prediction confirmed to the
  last rank: `b` only shifts constants. F4 wall and reduction counts sit
  inside run-to-run noise (8–50 ms, 10–21 reductions, 5–6/6 decomposable
  targets recovered).
- **Yield moves only with `|F|`.** Across each class `|F|` ranges 239–277
  (n = 17) and 495–529 (n = 19), i.e. ±8 %, and the exact 2-sum yield
  follows `|F|²`: 0.099–0.127 and 0.105–0.118. That is a factor ≤ 1.3 at
  best, chosen by the on-curve fraction inside `V`, not an exponent.
- **For the Koblitz campaign the class is rigid.** `K_0` has class number
  one for its maximal order, so the only modular-polynomial roots at
  `ℓ = 2` are itself (the Verschiebung) and its twist `K_1` (same `j`,
  not isogenous, smaller `|F|`); no `ℓ = 3` neighbour exists at n = 17 or 19.
  Any neighbour reached by walking down the volcano is not defined over
  `F₂` and loses the `τ`-orbit multiplier the Koblitz factor base relies on.
- The modular-polynomial table in `binary_isogeny` covers `ℓ ∈ {2, 3}`
  only; n = 23 was dropped because the full-descent Macaulay rank at
  degree 4 (46 variables) did not finish inside the run budget.

### 2.4 Prime fields

For the small-x factor base the enumeration PDP never sees the curve model,
and the yield `(2B)^m / (m!·p)` counts signed sums of a set of size `B`, so it
is model-independent up to coincidences. Isogenous prime-field curves cannot
change either; only extra automorphisms (`j = 0, 1728`) would, and those are
isogeny-class invariants the deployed curves lack.

## 3. Decision

PDP constants over prime fields are done (PR #1560); the binary algebraic
PDP has no cheap drastic lever left (G1 closed above), and the isogeny
hypothesis is answered in §2. The remaining large-factor candidate on the
binary side is the equivariant F4 (G2); on the prime side there is none
below the exponent barrier.
