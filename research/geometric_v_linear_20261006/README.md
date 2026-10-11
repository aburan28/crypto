# Linear (A, P) oracle over a geometric factor base

**Stage diagnostic.** `S`, end-to-end cost and speedup are **unset**. No work
below `2^61` is claimed.

**Class: accounting correction plus boundary.**
- It lowers this thread's own frontier.
- It closes the linearization class at three summands.
- It does not move any ratio to rho past the boundary.

## Answer

The frontier in
[`../linearization_reach_20260930`](../linearization_reach_20260930/README.md)
and [`../xl_bilinear_core_20261006`](../xl_bilinear_core_20261006/README.md)
was an artefact of using a **random** factor base `V`. With a **geometric**
`V = θ⟨1, g, …, g^{l−1}⟩`, a single linear solve resolves both summands of a pair
up to `l = ⌊(n+1)/3⌋`. That is 44 at `n = 131`, where random `V` stops at 14.

The `n = 131` gap moves as follows:

| | bits-only gap | charged gap, ω = 2 | charged gap, ω = 2.81 |
|---|---:|---:|---:|
| previous frontier (random `V`, `l = b = 14`) | 56.19 | 71.82 | 78.15 |
| **geometric `V`, subspace base, `m = 5`, `c = 28`** | **42.47** | **55.22** | **60.39** |

Both gaps are still positive: nothing here comes in under rho.

The binding constraint is now the linear-algebra cap `c ≤ 28–30`, not the
oracle. With a pair oracle the gap simplifies to `70.19 − c`. So no pair oracle,
at any per-trial cost, gets closer than about 40 bits.

A three-summand linear oracle, measured at `n = 131`, reaches only `l ≤ 5`.

## Derivation

Fix `x₁ = t` and put `A = X₂ + X₃` and `P = X₂X₃`. On `y² + xy = x³ + 1`:

    S₃(t, X₂, X₃) = (tA + P)² + tP + 1 = P² + t²A² + tP + 1.

This is linear over `F₂` in the bits of `(A, P)`, because squaring is
`F₂`-linear. `X₂` and `X₃` are the roots of `Z² + AZ + P`.

- `A ∈ V`.
- `P ∈ span(V·V)`. That span has dimension at most `2l − 1` for geometric `V`,
  and `l(l+1)/2` for random `V`.

So the unknown count `N` is `3l − 1` for geometric `V`, against `l(l+3)/2` for
random `V`.

The span bound was already measured in this repository and used inside a
Gröbner solve
([`../notes/koblitz-isogeny/subspace-structure-20261004.md`](../notes/koblitz-isogeny/subspace-structure-20261004.md)).
What is new here is the purely linear `(A, P)` solve and its effect on the chain
cost map. The derivation was proposed by an idea-generator subagent and checked
by hand before the protocol was frozen.

## Result 1: the law, the reach and the recovery (G1–G3 passed)

The confirmatory run used seed `20261008`, `n ∈ {19, 23, 29, 31}` and
`l = 4 … ⌊(n+1)/3⌋ + 2`, with both bases and 30 planted pairs per cell.
Trials were sized for about 24 expected successes, capped at `2^15` for
geometric and `2^13` for random. Every counted success is verified on the curve.

- **G1 passed.** The 8 geometric cells within reach that had at least 8
  successes deviate from `2l − n − 1` by −0.34 to +0.32 bits. At `n = 19`,
  where at least 3 such cells exist, the slope on `l` is 1.88.

  *Limitation:* the trial cap left few low-`l` cells with enough successes, so
  the slope test ran at `n = 19` only. The per-cell deviation test covers all
  four `n`.
- **G2 passed.**
  - Geometric `V` within reach has a mean rank defect of at most 0.86.
  - Past the reach, every cell of either base sits at generic rank: the defect
    is within +0.48 of `N − n`. The solution family is therefore `2^{N − n}`
    points, and `N` is what separates the bases. Geometric `V` grows it
    linearly (`3l − 1 − n`); random `V` grows it quadratically
    (`l(l+3)/2 − n`).
- **G3 passed.** All 22 geometric cells within reach recover 30 of 30 planted
  pairs.

The cells with at least 8 successes, plus geometric cells just past the reach:

| n | l | base | N | trials | successes | log₂ p̂ | law `2l−n−1` | mean defect | planted |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 19 | 4 | geometric | 11 | 32768 | 9 | −11.83 | −12 | 0.00 | 30/30 |
| 19 | 5 | geometric | 14 | 24576 | 22 | −10.13 | −10 | 0.03 | 30/30 |
| 19 | 5 | random | 20 | 8192 | 14 | −9.19 | −10 | 1.48 | 30/30 |
| 19 | 6 | geometric | 17 | 6144 | 23 | −8.06 | −8 | 0.24 | 30/30 |
| 19 | 6 | random | 24 | 6144 | 22 | −8.13 | −8 | 5.02 | 29/30 |
| 19 | 7 | geometric | 20 | 1536 | 33 | −5.54 | −6 | 1.48 | 30/30 |
| 19 | 8 | geometric | 23 | 384 | 29 | −3.73 | −4 | 4.05 | 30/30 |
| 23 | 7 | geometric | 20 | 24576 | 22 | −10.13 | −10 | 0.13 | 30/30 |
| 23 | 8 | geometric | 23 | 6144 | 19 | −8.34 | −8 | 0.86 | 30/30 |
| 23 | 9 | geometric | 26 | 1536 | 33 | −5.54 | −6 | 3.13 | 30/30 |
| 23 | 10 | geometric | 29 | 384 | 17 | −4.50 | −4 | 6.01 | 28/30 |
| 29 | 9 | geometric | 26 | 32768 | 10 | −11.68 | −12 | 0.12 | 30/30 |
| 29 | 10 | geometric | 29 | 24576 | 23 | −10.06 | −10 | 0.85 | 30/30 |
| 29 | 11 | geometric | 32 | 6144 | 19 | −8.34 | −8 | 3.12 | 30/30 |
| 29 | 12 | geometric | 35 | 1536 | 26 | −5.88 | −6 | 6.02 | 29/30 |
| 31 | 10 | geometric | 29 | 32768 | 8 | −12.00 | −12 | 0.24 | 30/30 |
| 31 | 11 | geometric | 32 | 24576 | 30 | −9.68 | −10 | 1.45 | 30/30 |
| 31 | 12 | geometric | 35 | 6144 | 27 | −7.83 | −8 | 4.06 | 30/30 |

The decomposition law is the same for both bases, as predicted: they find the
same decompositions. The difference is work per trial. Past the reach, a trial
must enumerate `2^{N − n}` candidates, for example `2^{5.0}` for random `V`
against `2^{0.2}` for geometric `V` at `n = 19, l = 6`.

## Result 2: n = 131 frontier (arithmetic, `gcost.py`)

The heuristics and budget convention are those of `costmap_a1.py`. A relation
needs `2^{n − 2c}` trials. Each trial is charged `ω · log₂(3c − 1)` bits, plus
the family size past the reach. The full tables are in `results/gcost.md`.

The best subspace cell is `m = 5`, `c = 28`: **42.47** bits-only and **55.22**
charged at ω = 2 (60.39 at ω = 2.81).

The orbit-union rows agree within 0.1 bit but are not headlined. An orbit-union
trial must place `X₂` and `X₃` in one Frobenius conjugate for `span(V·V)` to stay
small, and that case is not measured here.

Once a target decomposes with probability ≈ 1, the per-attempt budget is
`2^{RHO − c}`, where RHO = 60.8090 is the log₂ of matched rho. So the bits-only
gap of the pair oracle is exactly `n − RHO − c = 70.19 − c`. Linear algebra
(`m · 2^{2c} ≤ 2^{61}`) caps `c` near 30. **No pair oracle can close the gap.**

## Result 3: three summands do not linearize (K1–K4 passed, at n = 131)

`S₄ = G²` holds exactly on this curve, with `G = D² + D√σ₄ + σ₁σ₄ + σ₃` and
`D = σ₁ + σ₃`. It was verified as an identity before freezing: 20,000 of 20,000
random quadruples at each of `n = 13, 17, 23, 31` against the resultant, and
`G = 0` on all true 3-sums.

With `x₄ = x(R)`, `G` is linear in `e₁, e₂, e₃` and four product classes, which
gives `31l − 24` unknowns over geometric `V` (`k3.py`). The run used seed
`20261009`, `n = 131`, 40 planted triples per `l`:

| l | unknowns | full rank | recovered | rank |
|---:|---:|---:|---:|---|
| 2 | 38 | 40/40 | 40/40 | 38 |
| 3 | 69 | 40/40 | 40/40 | 69 |
| 4 | 100 | 40/40 | 40/40 | 100 |
| 5 | 131 | 13/40 | 13/40 | 128–131 (square system; random-matrix expectation ≈ 0.29) |
| 6 | 162 | 0/40 | 0/40 | 131 in 40/40 |

`G = 0` held on all 200 planted triples. The three-summand linear oracle stops
at `l = 5`, giving at most 10 net bits per trial against `2l ≤ 56` for the pair
oracle. **Linearization is closed at three summands, at cryptographic size.**

## Ideas considered this round

These came from an idea-generator subagent plus my own screening. Each faced one
kill test: resolve ≳ 56 extra bits per trial at polynomial cost, or change the
outer algorithm so that the gap formula no longer applies.

| idea | status | outcome |
|---|---|---|
| Geometric-`V` linear `(A, P)` oracle | **measured** (G1–G3) | frontier 56.19 → 42.47 bits-only; still fails |
| Three-summand symmetric linearization (`S₄ = G²`) | **measured at n = 131** (K1–K4) | reach `l ≤ 5`; closed |
| Large-prime variation (`X₃` in a larger space) | analytic | fails: reach drops once `X₃ ∉ V` (symmetry defect lost), and partial-relation cycles cost 2^{L/2}; ≈ 2^{117} vs 2^{97} per relation at `l = 30` |
| Subgroup factor base (closed under addition) | analytic | impossible: `#E = 4r` with `r` prime has no small subgroups |
| Batching the solve across trials | analytic | recovers at most the `ω·log₂ M` term (≈ 13 bits), not the 42-bit floor |
| `t` restricted to a subspace (batch over `t`) | analytic (subagent) | needs `4l < 5 − n`: impossible |
| Torsion / model / λ invariants | not run | on this curve they are pullbacks by endomorphisms `1 ± τ` and `[2]`, so ≤ 2 bits per summand (subagent derivation, unverified) |
| Graph-cycle IC (pair oracle, no linear algebra) | analytic (subagent) | removes the linear-algebra cap (`l → 44`) but costs about `2^{n−l} ≈ 2^{88}` tests with about `2^{43}` memory; gap about 26 bits; asymptotically `2^{2n/3}`, dominated by rho |
| Excess additive energy of geometric `V` | not run | long shot (prior ≈ 2%); a cheap test is specified in the subagent report |

## What would close the gap

Each trial must resolve `k` summands at polynomial cost with `(k−1)·c ≥ 70.19`.
With `c ≤ 28–30`, that means `k ≥ 4`.

- Linearization cannot supply it.
- Geometric `V` is already the best subspace for pairs: the subagent recalled
  (unchecked) an additive-combinatorics bound `dim(V·V) ≥ 2 dim V − 1`.
- Three summands linearize only to `l = 5` (Result 3).

**The open problem, stated precisely:** decide whether `R` is a sum of 4 points
from a factor base of dimension 24–30, at about `2^{12}` group-operation
equivalents per target, with a non-Macaulay method. Candidates named in the
subagent report are `F₂[z]` factorisation and Popov-form module reduction. None
is tested here.

## Files and reproduction

| file | role |
|---|---|
| `PROTOCOL.md`, `PROTOCOL-K3.md` | frozen before their runs: G1–G3 and K1–K4, with smoke tests disclosed |
| `glin.py` | the `(A, P)` oracle, geometric and random `V`, planted and random targets (reuses the frozen `lr.py`) |
| `run_confirm.py`, `score_confirm.py` | grid and scorer for G1–G3 |
| `k3.py` | the three-summand symmetric linear oracle at `n = 131` |
| `gcost.py` | the corrected `n = 131` cost map |
| `test_glin.py` | the `S₃` identity, the geometric span bound, planted recovery |
| `results/` | `confirm.json`, `confirm.log`, `score_confirm.json`, `k3.json`, `k3.log`, `gcost.json`, `gcost.md` |

```sh
python3 -m unittest test_glin -v
python3 run_confirm.py results 20261008 && python3 score_confirm.py
python3 k3.py 131 40 20261009 2,3,4,5,6 results/k3.json
python3 gcost.py results/gcost.json
```

The runs wrote into the session scratchpad and their outputs were copied into
`results/` unchanged. The first launch of the `k3.py` run started from the wrong
directory and executed nothing; it was relaunched with the frozen seed.

Pure Python 3, no dependencies.
