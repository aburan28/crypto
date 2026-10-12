# Run report: m = 3 curve-vs-presentation at n = 31 (PinX3), 2026-10-11 UTC 23:57

**Design.** Pre-registered in
`research/notes/koblitz-isogeny/m3-presentation-vs-curve-n31-design-20261012.md`.

**Command.** `KOBLITZ_PRES_M3=pin KOBLITZ_PRES_GEOMETRIC=1 KOBLITZ_PRES_FFD_TARGETS=0 koblitz_presentation_sweep 31 0 10 2000 32 1,2,3,4`.

**Code.** This PR's branch at `68349d4ea`, plus the sweep's `pin` mode.

**Host.** 4-vCPU shared cloud container; wall clock is indicative only.

**Records.**
- `experiments/koblitz_presentation/pres_31_0_l10_geom_m3pin.jsonl`: 128/128 cells,
  872 s.
- Analysis: `icc_n31_m3_geom.json` (pre-registered metrics) and
  `icc_n31_m3_geom_diag.json` (with `KOBLITZ_ICC_DIAG=1` diagnostics), both from
  `examples/koblitz_presentation_icc.rs`.

**Failures:** none. All 128 cells have `bsgs_ok` and 0 inconsistent relations.

## Values (measured)

20661 relations from 256000 probes; median 160 relations per cell and 16.6 ms per call.

| metric | ICC_curve | 95 % CI |
|---|---|---|
| yield per \|F\|³ | 0.000 | [0, 0.116] |
| ln µs per call | 0.048 | [0, 0.216] |
| work per call | 0.228 | [0.032, 0.403] |
| ln projected solve | 0.067 | [0, 0.219] |
| *diag:* \|F\| | 0.239 | [0.025, 0.437] |
| *diag:* pins per call | 0.230 | [0.035, 0.406] |
| *diag:* roots per pin | 0.000 | [0, 0.138] |
| *diag:* work per call ÷ \|F\| | 0.086 | [0, 0.297] |

![ICC](icc_n31_m3.svg)

`plot.py` draws from the JSON records only; it computes nothing new.
