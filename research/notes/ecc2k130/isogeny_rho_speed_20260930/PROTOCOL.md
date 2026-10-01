# Is Pollard rho faster on an ECC2K-130 isogenous curve?

**Status:** preregistered 2026-09-30, before measured class-member rho
rows and before any Modal/RunPod GPU launch.
**Question:** among curves `F_{2^{131}}`-isogenous to ECC2K-130, is there a
vertex where automorphism-aware Pollard rho is cheaper than on the
challenge curve `E_0: y² + xy = x³ + 1`?
**Unit:** `S = charged group operations / √r` (AGENTS.md §2).
**Reference:** signed-Frobenius rho on `E_0`, automorphism order `A = 2·131 = 262`.

This thread is about **rho**, not Gröbner/FFD. Prior notes
(`RESEARCH_ISOGENY_DEGREE_SEARCH.md`, `RESEARCH_ISOGENY_CLASS_SEARCH.md`,
`LEAF_HORIZONTAL_SEVEN_BRIDGE_20260925.md`) already closed related IC
questions; they do not answer the rho comparison.

## 1. Hypothesis and falsification

**Null (derived, stated before measuring).** No reachable isogenous vertex
beats `E_0` on end-to-end rho iteration count. Descending degree-263 leaves
lose the geometric 2-Frobenius endomorphism on a fixed model, so their
eligible automorphism order for a rational-point walk collapses from
`A = 2n` to `A = 2` (negation only). Expected walk length therefore rises by
`√n = √131 ≈ 11.45×`. With class number `h(O_K) = h(−7) = 1`, there is no
second `End = O_K` isomorphism class over the algebraic closure that could
restore cheap Frobenius while staying in the same isogeny class.

**Success (would reopen).** A measured, verified rho on a curve isogenous
to the same `(n, r)` instance with

- `S_candidate / S_E0 < 0.90` in the same operation unit, same seed policy,
  same distinguished-point rule family, and
- every recovered scalar checked as `[k]P = Q`,

on at least four independent targets, with transport cost charged when the
public instance is moved.

**Abandon.** Every measured non-Koblitz class member sits near the
`A = 2` floor (`S ≈ √(π/4) ≈ 0.886`) while `E_0` sits near
`S ≈ √(π/(2·2n))`, and no GPU step-cost measurement closes an `√n` gap.

Inadmissible: changing `r`; counting wall clock as the primary metric;
running leaf rho with a forged Frobenius that maps to a different curve;
dropping setup/verification; comparing multi-target amortized tables to a
one-target `E_0` walk; treating an unimplemented leaf GPU kernel as a win.

## 2. Boundaries (derived)

### 2.1 Floor — automorphism orders

| curve | End conductor | cheap geometric `τ: (x,y)↦(x²,y²)` on the model? | eligible `A` for rho | expected `S = √(π/(2A))` |
| --- | ---: | --- | ---: | ---: |
| `E_0` (challenge) | 1 (`O_K`) | yes | `2n` | `√(π/(4n))` |
| horizontal 263-loop | 1 | returns to `j = 1` | `2n` | same as `E_0` |
| descending 263-leaf | 263 | no: squaring leaves the model | `2` | `√(π/4) ≈ 0.886` |

At `n = 131`: `S_E0 ≈ 0.07743`, `S_leaf ≈ 0.886`, ratio `√131 ≈ 11.45`.

Prime extension degree matters: the only subfields of `F_{2^{131}}` over
`F_2` are `F_2` and itself, so the only `F_2`-rational ordinary model in
the challenge class is the Koblitz curve itself. Intermediate Frobenius
orders cannot appear by subfield definition of `a₆`.

### 2.2 Reference

Published ECC2K-130 figure: automorphism-aware rho
`√(π r / (2·262)) ≈ 2^{60.81}` iterations. This thread uses the same `A`
and prices every charged addition the repository's counted walk already
charges.

### 2.3 Class size (search cost, not a rho floor)

`H(Δ) ≈ 2^{65.06}` isomorphism classes. Exhaustive enumeration already
exceeds one rho solve; only the 263-neighborhood is a practical walk.
That neighborhood's non-horizontal vertices are the descending leaves
above.

## 3. Measurement plan

### 3.1 CPU class screen (primary empirical gate)

Sizes: `n ∈ {17}` from the frozen census in
`experiments/koblitz_isogeny_cost_sweep.json` (273 members, `a₂ = 1`
preferred family). Do **not** rescan; reuse the frozen `a₆` list.

Arms, same planted targets per member:

| arm | who | walk |
| --- | --- | --- |
| `sf` | Koblitz member `a₆ = 1` only | signed-Frobenius, `A = 2n` |
| `neg` | Koblitz + stratified sample of non-Koblitz members | negation, `A = 2` |

Sample: `a₆ = 1`, plus 8 non-Koblitz members at indices
`⌊i·(M−1)/8⌋` for `i = 1..8` in the sorted member list excluding 1.
Seeds: `2026093001 + 17·1000 + a₆ mod 997`. Three repetitions per
`(member, arm)`. Cap `max_steps = 200 · √(π r / 2)`.

Report one table: rows = members; columns =
`A`, `S`, `S_walk`, `steps/expected`, `verified`, `wall_ns` (secondary).
Classify each row advance / engineering / relabelling / accounting per
AGENTS.md §3. The prediction is **accounting**: non-Koblitz rows move to
the `A = 2` floor; ratio to the floor stays flat.

### 3.2 GPU (Modal and/or RunPod)

The campaign GPU stack under `ecc2k130/` is **Koblitz-only** (generated
`eccF131` / ONB walk). It cannot run a descending leaf. Scope:

1. Probe Modal and RunPod credentials in this environment; retain the
   receipt under `gpu/`.
2. If credentials work, run a short Koblitz GPU throughput smoke on the
   existing kernel as a **host practicality** note only. It does not
   compare isogenous curves.
3. A leaf GPU rho requires a new generic binary kernel. That is out of
   this protocol's budget. The automorphism floor already prices leaf
   rho `√n` worse in operations before any kernel constant.

Do not claim a GPU end-to-end speedup from a Koblitz-only smoke.

### 3.3 What this does not claim

- No ECC2K-130 discrete log.
- No IC / Gröbner crossover.
- No statement that every possible endomorphism-accelerated walk on a
  leaf is impossible: only that the **cheap geometric 2-Frobenius on a
  fixed leaf model** is absent, and that the public matched reference
  uses exactly that automorphism.

## 4. Deliverables

- This protocol (frozen before outcomes).
- `derive_boundary.py` — automorphism floors and `n = 131` table.
- `examples/isogeny_rho_compare.rs` + `evidence/n17-*.json` — measured
  screen.
- `gpu/ACCESS.md` + probe JSON — Modal/RunPod access outcome.
- `RESULT.md` — one table, one unit, verdict.

## 5. Cost accounting

CPU screen: single host, `RAYON_NUM_THREADS` recorded, release build.
Operation counts are primary; wall clock is secondary and not pooled
across contended runs. GPU: record SKU, driver, and that leaf curves were
not executed.
