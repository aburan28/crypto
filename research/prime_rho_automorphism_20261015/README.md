# Rho on the full automorphism group of the round-2 crater curves

**Date.** 2026-10-15. **Status.** Measured on public synthetic prime curves
(22–28-bit p). All records were verified and the audits passed. The wall-time
comparison was not admitted.

Round 2 of the conductor gap
([`../isogeny_conductor_gap_round2_20261011/`](../isogeny_conductor_gap_round2_20261011/README.md),
finding 4) matched index calculus against `rho.negation` only and noted that
"Rho gains √3 from the same automorphisms, which this repository does not yet
run on prime curves." This study runs that rho on the same twelve E2 crater
curves.

| file | what it is |
|:--|:--|
| [`PROTOCOL.md`](PROTOCOL.md) | question, arms and prediction, frozen before the run |
| `specs/` | `cubic.json`, `quartic.json`: E2's crater curves, targets and measurement settings unchanged |
| `sessions/` | ecbench sessions, audits, tables, saved comparisons |
| [`reports/run-20261015T2305Z/`](reports/run-20261015T2305Z/REPORT.pdf) | run report, graph (`ratios.svg`, `plot.py`) and PDF |
| `run_sessions.sh` | run and audit commands |

## Result

| comparison | mean-S ratio | 95% interval | predicted | sessions |
|:--|--:|:--|--:|:--|
| `rho.cubic` / `rho.negation` (j = 0, A = 6) | **0.611** | [0.505, 0.739] | 0.577 | `ECBS1h5e1b9912a46d` |
| `rho.quartic` / `rho.negation` (j = 1728, A = 4) | **0.791** | [0.673, 0.923] | 0.707 | `ECBS1hf8d1171784f1` |

Each row is 72 same-seed pairs over 6 curves × 4 targets × 3 rounds. The unit
is S (group-addition equivalents) and the intervals are cluster bootstrap.

1. **The automorphism rho is cheaper in S, and both predictions hold.** Each
   predicted value lies inside its interval. Pooled, the cubic walk needs
   about 1.64× fewer group-addition equivalents than the negation walk, and
   the quartic walk about 1.26× fewer.
2. **This is an operation-count diagnostic.** ecbench does not price the
   class-representative map (applying ζ₃ or ι and choosing a canonical
   element), so the ratio does not bound full cost. A priced comparison needs
   `price-units` support for prime curves.
3. **There is no wall-time result.** The host ran at L0 with about 82% foreign
   load. The descriptive wall medians (0.98 and 0.99) are dominated by
   per-solve overhead at these sizes.

## Consequence for round 2 (derived, not re-measured)

Round 2's crater rows compared index calculus with `rho.negation`. Against
the automorphism rho, every crater IC/rho ratio there grows by about the
factors above: about 1.6× on j = 0 curves and 1.3× on j = 1728 curves. Round
2's conclusion that no index-calculus cell is below its matched rho therefore
holds with a wider margin. This adjustment applies the pooled ratios to round
2's recorded S. It is not a joint re-run.

## Open

- A run on a quiet host (L2) for an admitted wall-time comparison.
- Pricing of the class-representative map on prime curves.
- A `benchmarks/` entry is needed once a wall-time or priced result exists.
  This S-diagnostic does not qualify.
