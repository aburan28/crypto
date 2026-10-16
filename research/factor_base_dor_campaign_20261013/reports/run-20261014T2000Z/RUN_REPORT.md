# Run report — factor-base → D* campaign, round 3 (run-20261014T2000Z)

**Engine:** `sage/factor_base_dor.sage` (unchanged since round 1; selftest
clean in round 1, no code changes after).
**Orchestrator:** `run_round3.sh` (TMPDIR redirect kept; a second pass
`round3b.log` re-ran the interrupted cells).
**Host:** same as rounds 1–2 (macOS arm64, single sage process per cell,
no pinning; deterministic counts, no wall-clock claims).
**Raw records:** 7 symmetric-S4 cells (n=21 l=2, n=19 l=2, n=23 l=2 with
32 targets, n=23 l=3, n=29 l=3 with 32 targets, n=29 l=4) and the n=23 m=2
top dims (9, 10; seeds 1–2). 102 JSON cells total across the study.

**Failures (preserved):** `sym4_n21_l4` — PolyBoRi/sage hard crash
("Python will now terminate", `raw/sym4_n21_l4.log`), twice attempted; the
n=21 δ-ladder keeps its l=2 and l=3 points. Round-3 pass 1 was interrupted
mid-cell (documented in `round3.log`/`round3b.log`).

**Measured values:** D* mean/min/max (dstar-sample 6–8 per cell),
point-decomposable fraction, boolean-unsat fraction, nvars, equation
degree — recorded in each raw JSON, aggregated in `figures/data.json`.

**Derived:** the `D* = l + 1` diagonal across all heavy-δ cells (Figure 4's
guide line); the orbit-amplification factor
`measured point fraction / 2^(3l−n)` (2^15 at n=29 l=3; 2^12 at n=23 l=2).

**Missing data:** sym4 n=21 l=4 (segfault); no cells past n=29; the m=83
gate remains unmet.
