# Run report: geometric-V cost diagnosis, m = 3 curve-vs-presentation (n = 19), PinX3 (2026-10-11, UTC 21:50 onward)

**Host.** Claude Code cloud container: 4 vCPU, x86_64 with PCLMULQDQ, shared, no core
isolation. Wall-clock values are indicative only.

**Code.** The branch on top of base `4b6e09b78` (this PR).

## Runs

1. **Geometric-cost diagnosis.**
   - Command: `KOBLITZ_WORDFIELD_SERIAL=<1|0> KOBLITZ_PILOT_M3=e3 koblitz_ic_pilot <n> 1 1 3 <l> 20 full [g1]`.
   - Cells: (n, l) ∈ {(17, 5), (19, 6)} × {monomial, geometric g1} × {serial, clmul}.
   - Records: `experiments/koblitz_m3/geom_cost_diag.jsonl`, 8 rows.
2. **m = 3 curve-vs-presentation**, pre-registered in
   `research/notes/koblitz-isogeny/m3-presentation-vs-curve-design-20261011.md`.
   - Command: `KOBLITZ_PRES_M3=e3 KOBLITZ_PRES_GEOMETRIC=1 KOBLITZ_PRES_FFD_TARGETS=0 koblitz_presentation_sweep 19 1 6 400 32 1,2,3,4`.
   - Cells: 128 of 128 completed in 1550 s.
   - Records: `experiments/koblitz_presentation/pres_19_1_l6_geom_m3e3.jsonl`.
   - Analysis: `experiments/koblitz_presentation/icc_n19_m3_geom.json`, from
     `examples/koblitz_presentation_icc.rs`.
3. **PinX3 vs e₃** on the same probes.
   - Command: `KOBLITZ_PILOT_M3=<e3|pin> koblitz_ic_pilot <n> 1 1 3 <l> 40 full g1`.
   - Records: `experiments/koblitz_m3/pin_vs_e3.jsonl`, 4 rows.
4. **PinX3 full solves.**
   - Command: `KOBLITZ_PILOT_M3=pin KOBLITZ_PILOT_MAX_TRIALS=400000 koblitz_ic_pilot <n> <a2> 1 3 <l> 50 full [g1]`.
   - Cells: (n, a₂, l) ∈ {(23, 1, 8), (31, 0, 10)} × {monomial, g1}.
   - Records: `experiments/koblitz_m3/full_solves_pin.jsonl`, 4 rows.

**Failures:** none. Every cell of runs 2 and 4 has `bsgs_ok` / `verified` true and 0
inconsistent relations.

## Values (measured)

### Run 1: ms per call

| n, l | V | serial | clmul | e₃ values/call | roots/call |
|---|---|---|---|---|---|
| 17, 5 | mono | 19.39 | 9.61 | 8192 | 1019.2 |
| 17, 5 | geom | 28.71 | 9.29 | 7790 | 966.0 |
| 19, 6 | mono | 211.24 | 106.47 | 65536 | 16409.6 |
| 19, 6 | geom | 503.19 | 122.77 | 65536 | 16345.6 |

### Run 3: ms per call, geometric V, 40 probes

| n, l | e₃ (clmul) | PinX3 | relations |
|---|---|---|---|
| 17, 5 | 9.702 | 0.052 | 2 / 40 for both |
| 19, 6 | 123.595 | 0.154 | 1 / 40 for both |

![per call](per_call.svg)

### Run 2: ICC_curve, 32 curves × 4 V

Pooled over all cells: 4180 relations from 51200 probes.

| metric | ICC | 95 % CI |
|---|---|---|
| yield per \|F\|³ | 0.098 | [0, 0.224] |
| ln µs per call | 0.011 | [0, 0.137] |
| work per call | 0.000 | [0, 0.051] |
| ln projected solve | 0.000 | [0, 0.069] |

![ICC](icc_n19_m3.svg)

### Run 4: PinX3 full solves

| n, a₂, l | V | seconds | probes | unknowns |
|---|---|---|---|---|
| 23, 1, 8 | mono | 0.5 | 418 | 118 |
| 23, 1, 8 | geom | 0.7 | 435 | 118 |
| 31, 0, 10 | mono | 33.5 | 5612 | 525 |
| 31, 0, 10 | geom | 45.9 | 5876 | 507 |

![full solves](full_solves_pin.svg)

`plot.py` reads the JSONL and JSON records and draws the SVGs; it computes nothing new.
