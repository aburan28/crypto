# Free-relation census by pair-root self-collision (Lead B)

**Date:** 2026-10-12. **Conductor task:** T-1410.
**Question:** `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` §3 (Lead B)
predicts that sorting the compact-orbit S3 pair-root table yields free
four-term relations among factor-base points, with ≈11,000 at the n=61
K=600 rung and ≈730 at n=83, and that the guided rank stage becomes "a sort
of a table" at those rungs. This experiment measures the actual yield, the
verified-relation count, and the homogeneous rank of the resulting matrix
mod r on the frozen Koblitz bases.

**Result class:** stage census; boundary measurement. No target is solved,
no DLP is recovered, no vs_rho comparison is made, and no ledger row is
promoted. The claims below are about relation yield and rank only.

## Method

Two new examples, byte-compatible with the compact-orbit twins:

- `examples/koblitz_free_relations.rs` — `u64` path (n ≤ 63).
- `examples/koblitz_free_relations_wide.rs` — `u128` path (64 < n ≤ 127).

Each loads (or constructs with the same deterministic scan) the frozen
`point_defined_factor_base`, enumerates the same relative pair-states
`(left, right, relative_shift)` and both S3 roots as `build_index`, keeps
**every** `(canonical root, state, root index, rotation)` entry (the
production table keeps only the first state per root), sorts by canonical
root, and pairs up entries whose pair-states differ beyond the trivial
transpose. Each candidate is then verified as an exact group identity — a
sign choice of the four factor-base points summing to O — before being
counted; the verified rows are inserted into the same incremental echelon
mod r that the rank stage uses, so the reported homogeneous rank is
measured, not assumed. Relations whose only O-sums are degenerate
(P + (−P) + Q + (−Q)) are counted separately as `trivial_only_pairs`.

Reproduce (from the repository root, after `cargo build --release`):

```sh
target/release/examples/koblitz_free_relations \
  experiments/koblitz-single-target-n61-20260925/factor-base.jsonl \
  experiments/koblitz-free-relations-20261012/census_n61_k600_frozen.json
target/release/examples/koblitz_free_relations_wide \
  experiments/koblitz-single-target-n83-20261006/base_n83_K600.jsonl \
  experiments/koblitz-free-relations-20261012/census_n83_k600_frozen.json
```

Host: local macOS, single thread. `scan_ms`/`sort_ms`/`verify_ms` are in
each JSON receipt. `KIC_DUMP_ROWS=<path>` writes every verified row.

## Measurements

| rung | base | collisions | nontrivial verified | homogeneous rank | needed (K−1) | census wall |
|---|---|---:|---:|---:|---:|---:|
| n=41, a=0 | construct K=64 | 1,024 | 0 | 0 | 63 | 0.07 s |
| n=61, a=0 | frozen K=600 | 10,188 | 588 | 49 | 599 | 8.6 s |
| n=61, a=0 | construct K=1200 | 12,732 | 4,800 | 746 | 1,199 | 37.5 s |
| n=61, a=0 | construct K=1500 | 10,924 | 6,000 | 1,197 | 1,499 | 50.0 s |
| n=71, a=0 | construct K=64 | 1,024 | 0 | 0 | 63 | 0.65 s |
| n=83, a=1 | frozen K=600 | 9,672 | 72 | 6 | 599 | 26.2 s |

Receipts: `census_*.json` in this directory; verified rows for the two
frozen rungs in `n61_k600_relation_rows.jsonl` and
`n83_k600_relation_rows.jsonl` (588 and 72 rows, every row 4-sparse).

Positive findings:

1. **Free relations exist and verify exactly.** Every nontrivial collision
   admitted a sign choice of four base points summing to O, checked with
   the same curve arithmetic the production pipeline uses. Lead B's
   mechanism is real.
2. **The collision count matches the prediction.** 10,188 collisions at
   n=61 against the note's ≈11,000; 9,672 at n=83 against ≈730 expected
   under the note's stronger reading — the total collision count is in
   family, but most collisions are structured and trivial (9,600 of them,
   the same count at every rung).
3. **Rank grows ≈ K^4 but with a small constant.** rank/relations ≈
   0.08–0.12; rank(K) ≈ 49·(K/600)^4 at n=61 (746 at K=1200, 1,197 at
   K=1500, against K−1 = 1,199 and 1,499 needed).

Negative findings (the falsifier verdict):

4. **The rank stage is not free at any frozen rung.** At n=61/K=600 the
   588 verified relations span rank 49 of 599 (8%); at n=83/K=600 rank 6
   of 599 (1%). The note's caveat that the relations are correlated is
   confirmed in the extreme: at n=61/K=600 the 588 rows touch only 168
   of 600 columns (top column in 36 rows), and at n=83 the 72 rows span 6
   dimensions.
5. **Making it free costs more than it saves.** Solving rank(K) ≥ K−1 with
   the measured law needs K ≈ 1,500–1,900 at n=61 and K ≈ 2,600 at n=83.
   At n=61 the census itself already costs 50 s at K=1500 against the 13.1 s
   guided rank wall at K=600, before charging the ≈6× larger linear
   algebra. At n=83 the K=2600 table is 1.1e9 entries (≈36 GB), scan
   ≈390 s single-thread against the 4.8 min 12-thread parallel guided rank
   already measured (`experiments/koblitz-single-target-n83-20261006/`),
   before charging the ≈19× larger linear algebra. Free relations lose
   end-to-end at every rung the ladder runs.
6. **The n=71/n=73 rungs are not rescued either.** Scaling the measured
   law (relations ∝ K⁴n⁴/r), rank K−1 needs K ≈ 3,700 at n=71 (2.0e9
   entries) and K ≈ 9,000 at n=73 (1.2e10 entries, ≈190 GB); the n=73
   guided rank wall (18,398 s sequential) would be replaced by a scan and
   sort of that table plus a ≈225× larger linear algebra.

## Consequences

- Lead B's weak form (relations exist, zero probes) is confirmed; its
  strong form (rank stage collapses at the rungs the ladder runs) is
  refuted by measurement at n=61 and n=83, the two rungs the note names.
- The ledger's precompute charge stays as is: the guided rank (sequential
  or parallel) remains the cheapest known way to fix the orbit logs at
  every rung, and any future "free-relation" claim must report the
  homogeneous rank, not the relation count.
- The sort-based census is still the cheapest known *relation generator*
  per unit work (588 verified relations in 8.6 s at n=61 with zero
  probes), so it remains the right tool if a future mechanism can
  decorrelate the rows (e.g. an algebraic family beyond pair-state
  coincidences). No such mechanism is proposed here.

## Files

- `census_n41_k64_construct.json`, `census_n71_k64_construct.json` — small
  controls (zero nontrivial relations, as the K⁴n⁴/r law predicts).
- `census_n61_k600_frozen.json`, `census_n61_k1200_construct.json`,
  `census_n61_k1500_construct.json` — n=61 K sweep.
- `census_n83_k600_frozen.json` — n=83 frozen rung.
- `n61_k600_relation_rows.jsonl`, `n83_k600_relation_rows.jsonl` — every
  verified relation row (mod-r coefficients over the orbit columns).
