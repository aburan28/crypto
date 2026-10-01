# Exact third-factor order screen on the archived degree-7 m3 bases

**Decision: do not spend the next batch experiment on the straightforward
exact-support third-factor reorder.** This retrospective screen reuses the
independently replayed [four-policy m3 panel](../m3_four_policy_20260930/RESULT.md).
It is a **probe-only diagnostic**, not a new held-out outcome, a complete
ECDLP comparison, or an ECC2K-130 method no-go. The original eight-point
source bases and their descendant-native pullbacks each have 36 unordered
pair-table entries and eight possible third factors. Transported and native
leaf partners give identical counts by the verified group map.

The unit below is a group addition for one third-factor residual probe. For
each archived base, [derive.py](derive.py) enumerates the exact set of subgroup
targets reachable by each third factor and optimizes their order over all
`8!` permutations. It evaluates the expected probe count on all 420 nonzero
subgroup targets, uniformly weighted; a miss also costs eight probes. An
independent subset dynamic program and exhaustive permutation enumeration
agree. The existing order is `0..7`. The straightforward exact-support
score uses 36 existing pair-table sums plus **288 additional group additions**
to evaluate eight thirds; the decision bound optimistically charges no set,
ranking, or memory cost. `L=1024` savings are the expectation times 1,024,
not a measured batch timing.

| Seed | Source policy (equal leaf partner) | Baseline additions / target | Optimal additions / target | Expected additions saved at L=1024 | Optimistic queries to repay 288-addition score |
|:--|:--|--:|--:|--:|--:|
| 2026093001 | original (transported) | 6.7548 | 6.7548 | 0 | none |
| 2026093001 | pullback (descendant-native) | 6.7524 | 6.7524 | 0 | none |
| 2026093002 | original (transported) | 6.7762 | 6.7667 | 9.75 | 30,240 |
| 2026093002 | pullback (descendant-native) | 6.7905 | 6.7667 | 24.38 | 12,096 |

Even if an oracle selected the best order separately for the already known A
and B orbit halves, its largest saving is 0.119 additions per target, needing
at least 2,420 targets to repay this score. This is a **post hoc upper bound**
for these exact bases and target distributions, not a candidate policy that
may look at held-out results. The fresh-Q batch size under consideration is
1,024, below all nonzero break-even counts here.

The same archive shows why new-base selection needs full cold accounting:
the observed source-base constructions cost 13,251 and 31,941 field
multiplications, while the four original scans through first full rank cost
4,114–12,826. An extra candidate base could exceed the entire observed rank
scan before it improves a witness. These are observed costs, **not** a lower
bound on every possible candidate generator or an impossibility result.

This screen does not include changes to the first-witness relation rows,
rank trajectory, or a cheaper approximate orderer. Those mechanisms remain
open. The next falsifiable direction is to improve **which useful points form
the factor base and whether their first-witness rows reach rank cheaply**,
while charging construction and selection in the one-target cold budget.
Any amortized batch claim must separately state its target count. A new
target-blind selector needs fresh seeds and disjoint target streams before
being evaluated. A Koblitz-family m31 exploratory comparison and the
repository's m83 confidence gate still precede n131 transfer. This note
changes no accepted method speedup or canonical scoreboard number.

The [machine-readable bound](BOUND.json) includes all 24 policy/universe
rows, source/result hashes, all exact probe totals, orders and break-even
counts. Reproduce from the repository root:

```sh
python3 research/notes/ecc2k130/m3_probe_order_bound_20260930/derive.py \
  --out /tmp/m3-probe-order-bound.json
cmp /tmp/m3-probe-order-bound.json \
  research/notes/ecc2k130/m3_probe_order_bound_20260930/BOUND.json
```
