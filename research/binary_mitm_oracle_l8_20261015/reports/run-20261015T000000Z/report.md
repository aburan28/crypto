# Run report: first binary-regime sub-2^(2ℓ) decomposition oracle at ℓ=8

**Run id:** `binary_mitm_oracle_l8_20261015` (example
`examples/binary_mitm_oracle_l8.rs`, `cargo run --release`).
**Date:** 2026-10-15. **Host:** arm64 macOS (Apple silicon, 14 logical
cores), unpinned, shared host — wall time is a practicality note; the
gate metric is the per-target median ratio of the two arms on the same
targets under the same contention.

## Question

Ledger regime A (binary non-Koblitz), stage `decomposition`: can an
oracle decide, in fewer than `2^{2ℓ}` operations per target, whether a
target decomposes over a dimension-`ℓ` subspace factor base, at the
fixed rung `ℓ = 8`, on ≤ 64 targets, with median cost at most half of
pairs-and-solve and zero disagreements against the group check?

## Setup (frozen)

- Field `F_{2^25}`, sparse irreducible modulus from
  `koblitz_index_calculus::find_irreducible_sparse(25)` (Rabin-checked
  in the builder). `n = 25` is odd so `y` lifts via the half-trace;
  `3ℓ = 24 ≤ 25` keeps the found/refuted mix genuine.
- Curve `y² + xy = x³ + a·x² + b` with **`a = z ∉ F₂`** (genuinely
  non-Koblitz) and `b = 1`.
- Factor base: all curve points with abscissa in the low-order subspace
  `⟨1, z, …, z⁷⟩` (dim `ℓ = 8`) — **256 points** (128 abscissae × 2 signs).
- Reference arm: `SubspaceOracle::low_order(8, 1, gf)` (Gaudry
  pairs-and-solve over the same subspace, same `b`; `Θ(2^{2ℓ})` pairs
  per target, 32,896 when refuting).
- Candidate arm: **group-law meet-in-the-middle pair table** — all
  `Pᵢ+Pⱼ` (`i ≤ j`) precomputed once (factor-base stage), query walks the
  `|F| = 256` base points, computes `R − P₃` and looks it up, for both
  signs of the target.
- Targets: 48 natural (random on-curve points, deterministic LCG seed
  `0x9E37_79B9_7F4A_7C15`) + 16 planted-satisfiable (sums of three base
  points).

## Measured (recorded in `results.json`)

| Quantity | Value |
|---|---|
| pair table entries | 32,896 distinct sums |
| table construction wall | 21.1 ms (run 1) |
| retained bytes (logged estimate) | ~1 MiB (32,896 × 32 B) |
| verdict mix | 20 found / 44 refuted (of 64) |
| disagreements MITM vs pairs-and-solve | 0 |
| witnesses re-added in the group | 20/20 verified (asserted per target) |
| median MITM per target | 0.120–0.155 ms (three runs) |
| median pairs-and-solve per target | 11.7–42.8 ms (same runs) |
| median ratio | 94×–277× (quiet host: 276.7×) |
| gate `median ≤ ½ pairs-and-solve` | **PASS** (all runs) |
| max counted MITM point-adds per target | 1,068 (= 2·2·|F| − slack; both signs) vs 32,896 pairs |

The counted-operation comparison is the decisive one: 1,068 point
additions and hash lookups against 32,896 quartic solves is a factor
`2^ℓ / 2` below the `2^{2ℓ}` enumeration frontier in the oracle's own
native unit, at the `2^{2ℓ}`-entry table cost charged to the
factor-base stage — exactly the trade the ledger's required
factor-base fields (`materialized`, `construction_wall_ms`,
`retained_bytes`) describe.

## Controls

- `#[cfg(test)]` exhaustive control at `n = 15`, `ℓ = 5`: the MITM
  decisions agree with the complete enumeration of all base-point triple
  sums (both signs) and with pairs-and-solve on 200 random targets
  (`cargo test --release --example binary_mitm_oracle_l8`).
- Every reported witness is re-added in the group and asserted equal to
  the target; refutations are exhaustive over the complete table.

## Claim boundary

Stage `decomposition` record only: constant-factor engineering at the
oracle level (`Θ(2^ℓ)` per-target at `Θ(2^{2ℓ})` table memory). No
key-recovery, no end-to-end ECDLP, no `vs_rho`, no asymptotic sub-ρ
claim — relation collection stays `Θ(2^n)` overall. FFD is logged as
`unmeasured` with reason (group-law oracle; no polynomial system is
built or solved), which the row's gate permits explicitly.

## Graph

`mitm_vs_pairs.svg` / `.png`: per-target wall times of both arms
(sorted, log scale) with the half-of-pairs-median gate line.
