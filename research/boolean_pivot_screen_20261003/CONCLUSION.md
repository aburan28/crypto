# Sibling restriction defeats full pivot-trace replay on the frozen discovery

The source-pinned native discovery completed all **80 public Boolean branch
pairs**. Every product row matched the independent set-parity oracle, the
compressed matrix corpus and all 19 manifest members replayed, and every pair
had the same labelled row skeleton. The preregistered feasibility gate is
**REJECTED: 0/64 nontrivial eligible pairs at n>=12** meet its combined
low-difference-rank and full-trace condition. No pair at any size even retains
the same pivot-column sequence. The four fresh holdout seeds were not run.

| Original variables n | Pairs | Base rank | Difference rank | Degree<=2 capacity | Difference/capacity | Difference/base | Same pivot columns | Same full trace |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 16 | 54–56 | 25–28 | 29 | 0.862–0.966 | 0.455–0.519 | 0 | 0 |
| 12 | 16 | 144 | 65–66 | 67 | 0.970–0.985 | 0.451–0.458 | 0 | 0 |
| 16 | 16 | 256 | 119–120 | 121 | 0.983–0.992 | 0.465–0.469 | 0 | 0 |
| 20 | 16 | 400 | 187–190 | 191 | 0.979–0.995 | 0.468–0.475 | 0 | 0 |
| 24 | 16 | 576 | 275–276 | 277 | 0.993–0.996 | 0.477–0.479 | 0 | 0 |

The table is copied from `discovery_01/result.json`. The ratio ranges are
descriptive exact rank quotients, not timing ratios. For each quadratic original
generator f and branch variable x, the Boolean difference between its x=1 and
x=0 restrictions is the affine derivative D_x f. The degree-3 Macaulay row
difference is `t * D_x f` for a multiplier t of degree at most one, so it lies
in the degree-at-most-two monomial space on n-1 variables. That space has
dimension `1 + (n-1) + C(n-1,2)`, shown above. The measured difference ranks
nearly saturate this allowed space. Low derivative degree therefore does **not**
give a low-rank matrix update in these fixtures.

The exact old pivot schedule cannot be replayed on the opposite branch: all
80 pairs have different pivot columns, before requiring the same source-row
choices. The native replay guard refused every pair; the independent fresh
reduction supplied the reference trace and rank. This rejects the specific
full-trace reuse mechanism under the frozen screen. It does not reject a graded
block algorithm: the cubic columns of the two matrices are identical because
their difference has degree at most two. Reusing that invariant high-degree
block, while recomputing and verifying the low-degree tail, is a distinct
mathematical hypothesis requiring a new protocol and matched complete-cost
measurement.

The retained [bundle](discovery_01/manifest.json) has manifest SHA-256
`d0a8e19fa107543e813bcc0d68f94b25b437478edb575ab2f8971c59b0737d00`,
result SHA-256 `5ecac8e7b7a56b2b7bb9aa0f79ebfd6726f13d74885f5736a82116e2819eb22e`,
and raw gzip SHA-256 `eeb1a711b7fbc4890525b78efc07c8c79e41e3ffc6c308140b5f39f76db955ec`.
It records source commit `25202d814`, compiler/host/binary hashes, 80 full
packed matrix pairs, ranks, pivot schedules, exact XOR-word counts and
payload-byte counts. At n=24 one labelled branch matrix contains 155,648
column-plus-row payload bytes; this excludes allocation metadata and is not
peak RSS. A changed structural result was rejected after its manifest hash was
updated, and an attempted holdout invocation rejected the failed discovery
before producing output. The original bundle was not modified.

Five optimized Rust tests pass, including a separate sparse-rank reducer on
paired matrices and differences, exact derivative checks, degree-drop fallback,
and pivot-disappearance rejection. The bundle's native verifier regenerated
every matrix and checked every file hash. No performance comparison was run.
Complete solver cost, relation yield, calibrated operations and rho ratio
remain null. There is no curve input, key-related result or cryptanalytic
breakthrough.
