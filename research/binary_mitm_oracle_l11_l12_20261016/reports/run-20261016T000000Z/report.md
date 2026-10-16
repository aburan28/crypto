# Run report: binary MITM oracle ladder, ℓ = 11 and ℓ = 12

**Run ids:** `binary_mitm_oracle_l11_l12_20261016/l11n33` and
`.../l12n37` (`cargo run --release --example binary_mitm_oracle_l8 <n> <ℓ> 48 16
research/binary_mitm_oracle_l11_l12_20261016`).
**Date:** 2026-10-16. **Host:** arm64 macOS (Apple M4 Pro), unpinned, shared
with other sessions' builds — wall time is a practicality note; the gate
metric is the per-target median ratio of the two arms on identical targets.

## Question

The promoted binary `decomposition` row's next target (set with the
ℓ = 9/ℓ = 10 ladder, PR #1944): extend to ℓ = 11 and ℓ = 12 — the table
reaching `2^{2ℓ} = 2^{22}` entries and the pairs-and-solve reference arm
reaching its historically usable ceiling (`ℓ ≈ 12`) — holding zero
disagreements and the counted sub-`2^{2ℓ}` property.

## Setup (frozen per rung)

- ℓ = 11: `F_{2^33}`, base 2,017 points over `⟨1, z, …, z^{10}⟩`; pair table
  2,033,876 distinct sums, 65.1 MB, 2.69 s build.
- ℓ = 12: `F_{2^37}`, base 4,023 points over `⟨1, z, …, z^{11}⟩`; pair table
  8,092,016 distinct sums, 258.9 MB, 14.6 s build.
- Same non-Koblitz curve `y² + xy = x³ + z·x² + 1` (a = z ∉ F₂, b = 1);
  same deterministic seed; 48 natural + 16 planted targets per rung;
  reference arm `SubspaceOracle::low_order(ℓ, 1)`.

## Measured (recorded in each `results.json`)

| ℓ | n | verdict mix | disagreements | max MITM ops | quartic solves | median MITM | median pairs | ratio | gate |
|---|---|---|---|---|---|---|---|---|---|
| 11 | 33 | 21 found / 43 refuted | 0 | 8,068 | 2,098,176 | 1.746 ms | 1,620.5 ms | 928.0× | PASS |
| 12 | 37 | 24 found / 40 refuted | 0 | 16,092 | 8,390,656 | 8.566 ms | 18,735.0 ms | 2,187.2× | PASS |

ℓ = 12 is the largest rung at which pairs-and-solve was historically
usable ("usable wall-clock to ~ℓ = 12", the pre-MITM ledger state); the
MITM arm decides the same 64 targets there for max 16,092 counted point
additions — a `2^ℓ / 2`-ish factor below the 8.4M quartic solves of the
reference arm — with the table charged to the factor-base stage (259 MB,
14.6 s, still cheap). Every witness is re-added in the group and asserted
equal to the target; MITM verdicts agree with pairs-and-solve on all 128
targets (0 disagreements).

## Claim boundary

Stage `decomposition` records only. Constant-factor engineering: the
`2^{2ℓ}`-entry table replaces the `2^{2ℓ}` per-target scan, so relation
collection stays `Θ(2^n)` overall. No key recovery, no end-to-end ECDLP,
no `vs_rho`, no asymptotic sub-ρ claim. FFD logged `unmeasured` with
reason (group-law oracle; no polynomial system), gate-permitted as at
ℓ = 8–10.

## Graph

`mitm_ladder_hi.svg` / `.png`: per-target wall times of both arms at both
rungs (sorted, log scale) with the half-of-pairs gate lines.
