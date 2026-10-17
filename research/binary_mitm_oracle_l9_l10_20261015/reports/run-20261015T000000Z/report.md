# Run report: binary MITM oracle ladder, ℓ = 9 and ℓ = 10

**Run ids:** `binary_mitm_oracle_l9_l10_20261015/l9n27` and
`.../l10n31` (example `examples/binary_mitm_oracle_l8.rs` with rung
arguments, `cargo run --release --example binary_mitm_oracle_l8 <n> <ℓ>
48 16 research/binary_mitm_oracle_l9_l10_20261015`).
**Date:** 2026-10-15. **Host:** arm64 macOS (Apple M4 Pro, 14 logical
cores), unpinned, shared with other sessions' builds — wall time is a
practicality note; the gate metric is the per-target median ratio of the
two arms on the same targets under the same contention.

## Question

The promoted binary `decomposition` row's next target (set with the
ℓ = 8 landing, PR #1919): does the group-law MITM pair-table design
extend to ℓ = 9 and ℓ = 10 on fresh frozen instances, holding zero
disagreements and the sub-`2^{2ℓ}` counted-operation property, with the
table staying cheap?

## Setup (frozen per rung)

- ℓ = 9: `F_{2^27}`, base 519 points over `⟨1, z, …, z⁸⟩`; pair table
  134,594 distinct sums, ~4.1 MiB, 601.6 ms build.
- ℓ = 10: `F_{2^31}`, base 1,053 points over `⟨1, z, …, z⁹⟩`; pair table
  554,304 distinct sums, ~16.9 MiB, 4,826.6 ms build.
- Same non-Koblitz curve `y² + xy = x³ + z·x² + 1` (a = z ∉ F₂, b = 1);
  same deterministic target seed `0x9E37_79B9_7F4A_7C15`; 48 natural + 16
  planted targets per rung; reference arm `SubspaceOracle::low_order(ℓ, 1)`.

## Measured (recorded in each `results.json`)

| ℓ | n | verdict mix | disagreements | max MITM ops | median MITM | median pairs | ratio | gate |
|---|---|---|---|---|---|---|---|---|
| 9 | 27 | 24 found / 40 refuted | 0 | 2,076 point-adds | 0.293 ms | 724.9 ms | 2,473.9× | PASS |
| 10 | 31 | 20 found / 44 refuted | 0 | 4,212 point-adds | 0.735 ms | 3,119.2 ms | 4,244.2× | PASS |

Counted units: 2,076 and 4,212 point additions + lookups per target
against `2^ℓ(2^ℓ+1)/2 = 134,912` and `525,312` quartic solves — the
`2^ℓ`-fold gap holds and widens in ratio as predicted by the design.
Every witness is re-added in the group and asserted equal to the target;
all MITM verdicts agree with pairs-and-solve (0 disagreements on both
rungs, 128 targets total).

## Claim boundary

Stage `decomposition` records only. Constant-factor engineering: the
`2^{2ℓ}`-entry table replaces the `2^{2ℓ}` per-target scan, so relation
collection stays `Θ(2^n)` overall. No key recovery, no end-to-end ECDLP,
no `vs_rho`, no asymptotic sub-ρ claim. FFD logged `unmeasured` with
reason (group-law oracle; no polynomial system), as at ℓ = 8.

## Graph

`mitm_ladder.svg` / `.png`: per-target wall times of both arms at both
rungs (sorted, log scale) with the half-of-pairs gate lines.
