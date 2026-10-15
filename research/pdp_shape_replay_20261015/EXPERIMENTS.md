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
| P2 | ~06:10 | reading | The headroom from existing records: roots are 27 % of inherited F4 on the dominant rung, and cubic `m = 3` roots are ~40 % rank-deficient | `research/notes/ecc2k130/RESEARCH_INHERITED_F4.md` §3.3–3.4 | Ceiling ≈ 10 % of the oracle's work on such cells | — |
| P3 | 06:20 | pilot | Frozen ladder (6 cells), default against `KIC_F4_ROOT_REPLAY=1` | `groebner_stage_bench --label default\|replay --out …` (uncommitted build, the code of the commit below) | Total word_ops 45,517,794 against 44,815,209: **1.016×**. Per cell 0.952–1.171×. Verdicts identical on every cell; replay adds splits on 3 cells (e.g. 212 against 186 splits on `K_0/2^9`, `m = 2`) | `pilot/frozen-default/`, `pilot/frozen-replay/` |
| P4 | 06:33 | — | Unit tests | `cargo test --release --lib inherited_f4` | 11 passed, 1 ignored, at `32bf9d358` | — |
| P5 | 06:42 | — | **Protocol registered** | `PROTOCOL.md` | — | — |
| C1 | 06:42 (under a minute) | claim | `chain`, `chain-holdout`, `holdout` ladders × {default, replay} | `run.sh` (binary `be69c5ca…`) | default ÷ replay **0.970, 0.954, 0.983**: replay costs more on all 9 cells run (0.753–0.985×); verdicts identical on all 9; more splits on 6. Cofactor-inadmissible cells skipped by the harness (`chain` 3 of 5 listed, `chain-holdout` 4 of 11). R1–R3 held; decision: recorded null | `raw/`, `summary.json` |
| D1 | 06:47–06:49 | diagnostic (after the claim) | Why: roots recorded / replayed / rows dropped, per cell, `chain-holdout`, replay arm | counting build (global counters + one line per cell; word operations and verdicts identical to C1's replay arm, checked) | large rank-deficient roots rarely repeat their layout (`K_0/2^15 m=4`: 18 recorded, 6 replayed); repeating roots drop ~4 rows (`K_1/2^15 m=4`: 1,832 replays, 3.8 rows each) | `diag/chain-holdout-replay-counted.json`, `diag/chain-holdout-replay-counts.txt` |
| D2 | 06:49 | diagnostic | D1 rerun, to check the counters are deterministic | as D1 | identical lines | `diag/chain-holdout-replay-counts.txt` (from this run) |
| L2 | 07:10 | reading | Lead 2 (weighted/block grading of the symmetrised model): already measured | `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` "The block structure, measured" (2026-09-12) | 1.00× where the minimum setting decides, 0.21× at `n = 15` otherwise: a recorded negative. Not rerun | — |
| L3a | 07:20 | — | Counters added to `decompose_by_e3` (field multiplications, word XORs); `koblitz_sym_s4` tests 5/5 | `cargo test --release --lib koblitz_sym_s4` | passed | — |
| L3b | 07:25 | — | **Lead 3 protocol registered** | `PROTOCOL-lead3.md` | — | — |

Times are file and commit timestamps. An earlier revision of this log gave
estimated times for P2–D2, up to twenty minutes late; they were corrected
the same morning, before review.
