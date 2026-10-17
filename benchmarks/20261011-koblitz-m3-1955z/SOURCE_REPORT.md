# Run report: m = 3 symmetric-S₄ solvers (2026-10-11, UTC 19:55)

**What ran.** Two runs, both through `examples/koblitz_ic_pilot.rs`, built from the
branch at base `4b6e09b78` plus this PR's changes.

1. **Full solves.** `KOBLITZ_PILOT_M3=e3 KOBLITZ_PILOT_MAX_TRIALS=200000
   koblitz_ic_pilot <n> 1 <a6> 3 <l> 50 full [g1]`.
   - n = 17, l = 5 and n = 19, l = 6.
   - The first 4 curves of each existing presentation sample (a₂ = 1).
   - Each curve with the monomial V and with the geometric V seed 1.
   - 16 cells, run sequentially.
2. **Solver comparison.** `KOBLITZ_PILOT_M3=<e3|sym|chain> koblitz_ic_pilot 17 1 1 3 <l> 4 full`.
   - l = 3, 4, 5; 4 probes per point; per-cell timeout 300 s.

**Host.** Claude Code cloud container: 4 vCPU, x86_64, shared, no core isolation.
Wall-clock values are indicative only.

**Records.**
- `experiments/koblitz_m3/full_solves_e3.jsonl`: 16 lines, one full pilot row each.
- `experiments/koblitz_m3/solver_compare_n17.jsonl`: 9 lines, one of them a timeout
  record (`sym`, l = 5, no row produced in 300 s).

## Values (measured)

Solver comparison, seconds per decomposition call (`us_per_call`), n = 17:

| l | e₃ enumeration | symmetric Gröbner | S₃ chain |
|---|---|---|---|
| 3 | 0.000116 | 0.0106 | 1.61 |
| 4 | 0.00141 | 4.10 | 9.21 |
| 5 | 0.0216 | timeout (> 300 s for 4 probes) | 45.2 |

No probe produced a relation at these 4-probe points; the predicted yield is
≤ 0.034 per probe.

![per-call time](per_call_n17.svg)

Full solves: 16 of 16 completed, and in every cell `log == planted == bsgs_log` and
`verified = true`. Wall seconds per solve:

| | monomial V | geometric V (g1) |
|---|---|---|
| n = 17, l = 5 (4 curves) | 3.5, 5.3, 7.9, 3.9 | 3.9, 9.7, 11.0, 14.5 |
| n = 19, l = 6 (4 curves) | 93.6, 104.2, 73.7, 78.8 | 287.3, 260.9, 168.3, 124.8 |

![full solves](full_solves.svg)

Per-call time for the e₃ engine (measured): about 16 ms (monomial) and 28 ms
(geometric) at n = 17; about 0.19 s and 0.47 s at n = 19.

## Failures

One cell timed out: symmetric Gröbner at l = 5. That records a cost, not a mathematical
result. No other failures.

## Graph sources

`plot.py` reads the two JSONL files and draws the SVGs; it computes nothing new.
