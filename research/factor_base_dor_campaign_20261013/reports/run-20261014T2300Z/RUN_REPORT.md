# Run report — factor-base → D* campaign, round 4 (run-20261014T2300Z)

**Engine:** `sage/factor_base_dor.sage` (unchanged since round 1).
**Orchestrator:** `run_round4.sh` (TMPDIR redirect kept).
**Host:** same as rounds 1–3 (macOS arm64, single sage process per cell, no
pinning; deterministic counts, no wall-clock claims).
**Raw records:** 13 cells: yield-precision symmetric-S4 cells
(l=2 at n=23, 29, 31, 37, 41 with 256 targets; l=3 at n=29, 37 with
128 targets; `--dstar-sample 4`, dmax 5) and the m=2 grid at n=29
(dims 4–9, standard/random, seeds 1–2, 16 targets). 121+ JSON cells in
the study; round 4 completed with zero failures.

**Measured values:** point-decomposable fraction at 1/256 and 1/128
resolution per cell (the round's primary observable), D* per dstar-sample,
nvars, boolean-unsat fraction, plus the m=2 grid aggregates.

**Derived:** the amplification identity
`A = pfrac · 2^n / 2^(3l) = 2^(n−13)` (random, l=2) and `2^(n−14)`
(geometric, l=2) — i.e. pfrac constant in n; the D\*=3 reading at
n=37, l=3 (below the round-3 l+1 diagonal); Figure 5's flat yield curves.

**Missing data:** standard bases yield zero point decompositions at every
cell (recorded, unexplained — a trace/parity obstruction is the suspected
cause, not measured); no cells past n=41; m=83 gate unmet.
