# Experiment log: root shape replay in the inherited-F4 oracle

This lists every run made for this study, in order, including the pilot and
the survey that chose it. It is append-only. Times are UTC, 2026-10-15.
Host: Claude Code cloud container, Xeon @ 2.10 GHz, 4 vCPU. The machine was
shared with the characteristic-2 ecbench sessions
(`research/gaudry_quartic_d4_char2_20261015`, PR #1883); counts are
deterministic, and wall time is descriptive only.

| # | when | status | what | command | outcome | data |
|--:|:--|:--|:--|:--|:--|:--|
| P1 | 05:45 | survey | Read-only survey of the prime-degree binary IC oracles and of where five solver levers from the `k = 4` work could apply | an Explore agent over `src/cryptanalysis/koblitz_*`, `inherited_f4.rs`, `binary_semaev*`, `ic_framework/`, `research/pdp_*`, `research/gf2_*` | Only cross-target shape replay has headroom; weighted grading, the chain eigen step, sparsest pivots (rejected 5× in `research/gf2_*_pivot_20260929`) and Frobenius reduction cannot move cost. Table oracles already lead the algebraic ones by orders of magnitude and still lose to rho | summarised in the PR |
| P2 | 06:15 | reading | The headroom from existing records: roots are 27 % of inherited F4 on the dominant rung, and cubic `m = 3` roots are ~40 % rank-deficient | `research/notes/ecc2k130/RESEARCH_INHERITED_F4.md` §3.3–3.4 | Ceiling ≈ 10 % of the oracle's work on such cells | — |
| P3 | 06:18–06:20 | pilot | Frozen ladder (6 cells), default against `KIC_F4_ROOT_REPLAY=1` | `groebner_stage_bench --label default\|replay --out …` (uncommitted build, the code of the commit below) | Total word_ops 45,517,794 against 44,815,209: **1.016×**. Per cell 0.952–1.171×. Verdicts identical on every cell; replay adds splits on 3 cells (e.g. 212 against 186 splits on `K_0/2^9`, `m = 2`) | `pilot/frozen-default/`, `pilot/frozen-replay/` |
| P4 | 06:25 | — | Unit tests | `cargo test --release --lib inherited_f4` | 11 passed, 1 ignored, at `32bf9d358` | — |
| P5 | 06:50 | — | **Protocol registered** | `PROTOCOL.md` | — | — |
| C1 | 06:50– | claim | `chain`, `chain-holdout`, `holdout` ladders × {default, replay} | `run.sh` | running | `raw/` |
