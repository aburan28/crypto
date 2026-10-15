# Experiment log: `k = 4` in characteristic 2

This is every run made for this study, in order, including pilots, scratch
runs, failures and the survey. It is append-only and updated in the same PR
as each run lands. Times are UTC, 2026-10-15. Host: Claude Code cloud
container, Xeon @ 2.10 GHz, 4 vCPU.

Statuses:
- **pilot**: before the protocol; it chose options or predictions, and no
  claim is made from it.
- **claim**: registered in `PROTOCOL.md` before it ran.
- **survey**: reading only, no run.

Ordinary unit-test runs during development are listed once, as E8.

| # | when | status | what | command / spec | code | outcome | data |
|--:|:--|:--|:--|:--|:--|:--|:--|
| E1 | ~04:00 | pilot | First build of `gaudry_quartic_bin`, `m = 5, 8, 10`, weighted (2,2,2,1) + sparsest + replay, 3 planted targets each | `examples/gaudry_quartic_bin --m M --seed 1 --targets 3 --planted --weighted --sparsest --replay` | uncommitted, before `5a5828371` | `m = 5`: 2 of 3 did not close (dimension 568), 1 solved at `C₄ = 1.10·10⁶`, exact. `m = 8, 10`: **panicked**: no irreducible `z⁴ + c₁z + c₀` (an affine linearised polynomial in characteristic 2). Fixed by adding a `z³` term | not kept (terminal only) |
| E2 | ~04:10 | pilot | The same after the fix, `m = 8, 10` | as E1 | uncommitted | 3/3 closed at each `m`, dimension 512, exact against MITM; `C₄` 3.6·10⁶ (`m = 8`), 6.9·10⁶ (`m = 10`); elimination only ~6·10⁵, against 4.9·10⁷ for the same matrix shape in odd characteristic | not kept |
| E3 | ~04:15 | pilot | Term count of the interpolated `H` (an ignored test, since removed) | `cargo test --release --lib shape_probe -- --ignored --nocapture` | uncommitted | `m = 5` and `8`: **21 of 105** monomials carry `H`, 20–21 terms per Weil equation. Every present monomial has even exponents in `s₂, s₃, π`, which led to the Frobenius form | output quoted in `PROTOCOL.md` and the code's docs |
| E4 | ~04:25 | pilot | Frobenius form (2,2,1,1), degree 11, `m = 5, 8, 10`, 3 planted | `… --frobenius --sparsest --replay` | uncommitted | `m = 8, 10`: 3/3 closed, dimension 64, 387 × 448, exact; `C₄` 1.9–2.1·10⁵. `m = 5`: the same 2 of 3 did not close (dimension 67) | not kept |
| E5 | ~04:30 | pilot | Is the `m = 5` failure `t_R ∈ F_{q²}`? 30 planted per `m = 5, 6, 7` and 40 uniform per `m` | as E4 with `--targets 30 --planted` / `--targets 40`, subfield diagnostic | uncommitted | Hypothesis **rejected**: `t_R ∉ F_{q²}` in every failure. Planted failures 7/30, 5/30, 1/30; uniform 0/40 at each `m` | `pilot/planted-m5-frobenius.jsonl` (the closed `m = 5` rows) |
| E6 | ~04:35 | pilot | Second hypothesis: every decomposition of a failing target contains a pair `{t, t + 1}` (a positive-dimensional family) | as E5 with MITM on failures | uncommitted | **Confirmed** on all 7 failures at `m = 5` (12/12 or 7/7 decompositions paired) | log lines quoted in `PROTOCOL.md` |
| E7 | ~04:40 | pilot | Cost only at `m = 12, 14, 15`, 3 planted, no MITM | `… --frobenius --sparsest --replay --no-mitm` | uncommitted | `C₄` 2.3·10⁵, 2.5·10⁵, 2.6·10⁵; only root finding grows with `m` | not kept |
| E8 | 04:00–05:20 | — | Unit tests during development | `cargo test --release --lib gaudry_quartic_bin` (and `ecbench`, `ic_framework`) | up to `5a5828371` | All passed at the commit: 17 `gaudry_quartic*`, 79 `ecbench`, 49 `ic_framework`. No failing test run occurred | — |
| E9 | ~04:50 | pilot | ecbench plumbing, `m = 5`, rho + IC Frobenius, `max_trials = 30` | `pilot/scratch-spec-m5.json` with 30 trials, session `ECBS1h3e9aeab53ee8` | uncommitted | rho verified; IC exhausted, as expected at 30 trials | **overwritten** by E10 (same output directory) |
| E10 | ~04:55 | pilot | Same, `max_trials = 20,000` | `pilot/scratch-spec-m5.json`, `ECBS1h69d8f27b06bc` | uncommitted | IC verified: **S = 930** against rho 1.23; 741 trials, 11 relations, 2 unsolved systems, 711 replays / 28 fallbacks | `pilot/scratch-session-m5/` |
| E11 | ~05:05 | pilot | Same at `m = 8` | `pilot/scratch-spec-m8.json`, `ECBS1h7555b4beea76` | uncommitted | IC verified: **S = 607** against rho 0.96; 15,518 trials, 60 relations (~259 trials per relation, not the 67 of `m = 5`) | `pilot/scratch-session-m8/` |
| E12 | ~05:08 | pilot | Same at `m = 10` | `pilot/scratch-spec-m10.json`, `ECBS1hec6c4f0cc48f` | uncommitted | IC **exhausted** at 20,000 trials with 100 relations (needs ~256). The registered specs use 5,000,000 | `pilot/scratch-session-m10/` |
| E13 | 05:29 | — | **Protocol registered** | `PROTOCOL.md`, `specs/`, `run_ladder.sh` | `666438375` | — | — |
| E14 | 05:30–05:45 | claim | Part A ladder and the `m = 8` ablation (22 configurations) | `RAYON_NUM_THREADS=3 run_ladder.sh` | `5a5828371`, example sha256 `25a20a67…` | All finished. 0 failed solves against MITM; planted unsolved 5 (`m = 5`), 3 (`m = 6`), 0 after; uniform unsolved 0 | `raw/*.jsonl`, `raw/*.log`, `raw/summary.json` |
| E15 | 05:30–08:00 | claim | Session B1 `ECBS1h03a1d62aba88`, `m = 6`, three forms + rho | `specs/SPEC-char2-forms-m6.json` | recorded commit `666438375`, ecbench sha256 `92610bdf…` | 10/10 verified. `S` plain 713,704 / 722,233, weighted 38,542 / 31,128, Frobenius 4,501 / 3,654, rho 1.017 / 1.314; identical trials across IC arms | `sessions/char2-forms-m6/` |
| E16 | 08:02–08:49 | claim | Session B2 `ECBS1h4d63e1551dc1`, `m = 6, 8, 10, 12`, Frobenius + rho | `specs/SPEC-char2-ladder.json` | recorded commit `f0a8f49d9` (data and docs only since B1), same binary | 24/24 verified. IC ÷ rho theory 5,080 / 4,124, 648 / 519, 80.6 / 82.9, 49.7 / 48.1 | `sessions/char2-ladder/` |
| E17 | 08:52–09:54 | claim | Session B3 `ECBS1hb9229e3a0f30`, `m = 14`, Frobenius + rho, one target | `specs/SPEC-char2-m14.json` | as E16 | 3/3 verified. IC `S` 5.370, rho 1.258; IC ÷ rho theory **6.06**; IC wall 2,903 s against 168 s (L1, descriptive) | `sessions/char2-m14/` |
| E18 | 05:45 | survey | Read-only survey of prime-degree binary IC oracles (for the second half of the task) | an Explore agent over `src/cryptanalysis/koblitz_*`, `binary_semaev*`, `ic_framework/`, `research/pdp_*` | — | Only cross-target shape replay has headroom there; weighted grading, the chain step, sparsest pivots (rejected 5× before) and Frobenius reduction cannot move cost. Root-finding table oracles already lead and still lose to rho | summarised in the PR |
| E19 | 09:47–11:41 | audit | Replay audit of B2 | `ecbench verify --replay-all` | — | OK: 24 replays reproduced, 0 problems | `audits/char2-ladder.*` |
| E20 | 11:41–11:46 | audit (killed) | Replay audit of B1, chained after E19 in one background shell | as E19 | — | **killed** 5 minutes in, when the shell hit the tool's 2-hour background limit; no output. Rerun as E21 | `audits/char2-forms-m6.verify.log` (overwritten by E21) |
| E21 | 11:46– | audit | Replay audits of B1 and B3, detached | as E19 | — | AUDIT_PENDING | `audits/char2-forms-m6.*`, `audits/char2-m14.*` |
| E22 | ~11:50 | analysis | Summary (intervals, paired ratios, trend fit, model), `ecbench compare` (6 comparisons saved) and `ecbench table` | `gaudry_quartic_bin_summary`, `ecbench compare --save`, `ecbench table` | — | slope −0.248, `n*` ≈ 2^68.1 (leave-one-out 2^66.7–2^71.6); model 2^63.1–2^69.0; B1–C2 held | `raw/summary.json`, `raw/table.md`, `sessions/*/comparisons/` |

## Part A results (claim runs; Part B is in README.md)

- **E14 ablation (`m = 8`, `C₄`, `F_q` multiplications):**
  - plain + sparsest + replay: 9.54·10⁷;
  - weighted: 3.60·10⁶ (**26.5×** below plain, 26.4–26.6);
  - Frobenius: 1.98·10⁵ (**18.1×** below weighted, 18.08–18.18; **481×** below plain overall).
  - Sparsest pivots give 1.33× in the Frobenius form; replay gives 1.0003×.
- **E14 ladder, Frobenius + sparsest + replay, replayed mean `C₄`:**
  - `m = 5`: 1.13·10⁵ (uniform) and 1.23·10⁵ (planted);
  - `m = 6`: 1.55–1.68·10⁵;
  - `m = 7`: 1.81–1.88·10⁵;
  - `m = 8`: 1.90–1.99·10⁵;
  - `m = 10`: 2.07–2.16·10⁵;
  - `m = 12`: 2.23–2.33·10⁵;
  - `m = 15`: 2.50·10⁵ (cost only).
  - Every closed solve with MITM agrees.
