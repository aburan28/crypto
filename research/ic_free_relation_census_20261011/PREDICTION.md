# Free-relation census of the K = 600 pair-root tables: prediction (written before the fixture runs)

Lead B of [`RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`](../../RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)
§3: two pair-states with the same Frobenius-canonical x-coordinate give a 4-term relation
among the `K` orbit representatives, found with zero probes.

## The note's prediction and its falsifier (quoted)

- **Count:** `≈ n⁴K⁴/r`, which is ≈ 11,000 at n = 61 and ≈ 730 at n = 83.
- **Falsifier:** "≥ 500 verified relations at both rungs and full column rank at n = 61".

## A competing count, and why

The note bins `M ≈ K²n` pair-sum **classes** into `≈ r/2` x-values. But a class is a
Frobenius orbit, and the subgroup has `≈ r/(2n)` orbits of x-values. So the coincidence count
is `M²·n / r = K⁴n³/r`, smaller by a factor of `n`:

| fixture | `log₂ r` | note `n⁴K⁴/r` | classes `K⁴n³/r` | relations needed (rank `K − 1`) |
|:--|--:|--:|--:|--:|
| n = 61, a = 0 | 47.2 | ≈ 11,000 | ≈ 180 | 599 |
| n = 73 | 56.3 | ≈ 40 | ≈ 0.6 | 599 |
| n = 83, a = 1 | 52.9 | ≈ 730 | ≈ 9 | 599 |

## Calibration (synthetic bases, seen before this file was written)

[`examples/koblitz_free_relation_census.rs`](../../examples/koblitz_free_relation_census.rs)
`construct:<n>:<a>:<K>:<seed>` builds `K` seeded subgroup points. Records are in
`runs/calibration.jsonl`.

| n | K | `log₂ r` | nontrivial relations | `K⁴n³/r` | `n⁴K⁴/r` | rank |
|--:|--:|--:|--:|--:|--:|--:|
| 37 | 60 | 27.8 | 2,698 | 2,847 | 105,329 | 59 = K − 1 |
| 41 | 120 | 39.0 | 33 | 26 | 1,066 | 11 |
| 43 | 150 | 32.1 | 8,585 | 8,667 | 372,674 | 149 = K − 1 |

The n = 23, 29 and 31 bases are saturated, because the relation count is capped by the
states, so they test only the logic. Every collision verified as a group identity, and every
kept row verified as `Σ [c]P = O`.

## Prediction for the fixtures

| fixture | relations predicted | rank predicted | the note's falsifier |
|:--|:--|:--|:--|
| n = 61 | ≈ 180 (Poisson, so roughly 150–210) | ≤ 180 < 599 | fails at n = 61 |
| n = 73 | 0–3 | — | fails |
| n = 83 | ≈ 9 (roughly 3–17) | — | fails |

- **If the counts land near `n⁴K⁴/r` instead, this prediction is wrong.**
- Pair sums of these fixtures stay in the order-`r` subgroup, because the representatives
  are checked to be in the `λ`-eigenspace. So the subgroup bin count applies.
