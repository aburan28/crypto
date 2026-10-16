# N83 precomputed-table rho control (2026-10-16)

**Verdict:** under the decision rule frozen in [PROTOCOL.md](PROTOCOL.md), neither
N83 "IC beats rho" online reading survives a rho given matched precomputation
(Bernstein–Lange table of `G`-only distinguished points). In every time-matched and
memory-matched comparison the 95% interval of the median paired ratio
IC/precomputed-rho lies entirely above 1, so precomputed rho is faster online.
Report: [REPORT.pdf](REPORT.pdf) ([source](REPORT.typ)). Run reports:
[run 1](reports/run-20261016T0524Z/report.pdf), [run 2](reports/run-20261016T0909Z/report.pdf).

| claim | matched arm (run, runner) | IC online median ms | precomputed rho median ms | median IC/rho (95% CI) |
|---|---|---:|---:|---|
| K3300 Tier 1, 113× | (a′) table walks = IC setup wall 680 s, 1 thread (run 2, EPYC 7763) | 32.78 | 8.11 | **3.74** (2.65–11.05) |
| K3300 Tier 1, 113× | (b) table = IC max RSS 10.79 GB (run 2, EPYC 7763) | 32.78 | 1.17 | **40.61** (24.27–71.72) |
| K600 hosted, 12.29× | (a′) = IC setup wall 772 s (run 2, Xeon 6973P-C) | 563.9 | 5.17 | **81.52** (46.35–236.59) |
| K600 hosted, 12.29× | (a) 650 s ≤ IC setup wall 785 s (run 1, EPYC 7763) | 580.3 | 3.33 | **214.73** (47.94–483.84) |
| K600 hosted, 12.29× | (b) table = IC max RSS 5.75 GB (run 2, Xeon) | 563.9 | 1.19 | **395.53** (177.93–1134.45) |

16 public points each (the eval points of hosted run `gha-38597929540-1`); every
recovered scalar of every arm replayed with the generic library curve, all arms
agree, and all match the hosted run's answers (`all_scalars_verified_and_agreeing:
true` in all four `analysis.json`). Positive control, same runners: empty-table fold
rho / IC = 70.9× (19.2–92.4) for K3300 and 3.62× (0.93–4.48) for K600 (run 2).
Operation counts (§5a): no ranking admitted (K600 order flips inside the 0–0.1 GAE
band; K3300 IC emits no root-call count, so its count is a lower bound).

## Layout

- [PROTOCOL.md](PROTOCOL.md) — question, arms, decision rule, prediction (frozen before timing) and amendment 1 (before run 2).
- `targets.jsonl`, `targets.sha256`, `reference_scalars.txt` — inputs (answers read only by the analyzer).
- `configs/` — precomputed-rho configs; `run_job.sh` — thin per-job orchestration.
- `records/run-38855595549/`, `records/run-38878565324/` — raw records per job (`ic/`, `rho_plain/`, `bl_a/`, `bl_a2/`, `bl_b/`, `host/`), `analysis.json` from `examples/n83_precomputed_rho_analyze.rs`, CI test logs.
- `figures/` — graphs ([plot.py](plot.py) lays out `analysis.json` values only); `REPORT.typ`, `tables.typ`, `verdict.typ`, `notes.typ` — report source.
- Code: `examples/n83_precomputed_rho.rs`, `examples/n83_precomputed_rho_analyze.rs`, `StrongRhoG::dp_walks` in `src/cryptanalysis/koblitz_strong_rho.rs`; workflow `.github/workflows/n83-precomputed-rho-control.yml`.

## Deviations (also in the report)

1. One IC process per job takes all 16 targets (setup once, per-target online clock).
2. Run 1's k3300 job ran on an EPYC 9V45 where IC's setup wall (498 s) was below arm (a)'s 650 s; that arm is excluded from the verdict and amendment 1 added (a′), budgeted at the measured IC setup wall.
3. Run 1's in-job analysis failed (IC summary is on stdout); run 2's workflow re-analysed run 1's raw artifact with the fixed analyzer.
4. Run 2's first `tests` attempt failed two `jv_cover` tests on wall-clock timeouts in untouched code; run 1's identical library passed 4,146/4,146.

## Reproduce

Push to `research/n83-precomputed-rho-control` (or `workflow_dispatch`) runs
`.github/workflows/n83-precomputed-rho-control.yml`; locally:
`run_job.sh <k3300|k600> <ic-bin> <rho-bin> <bl-bin> <analyze-bin> <out-dir>` with
binaries built from the commits named in the protocol.
