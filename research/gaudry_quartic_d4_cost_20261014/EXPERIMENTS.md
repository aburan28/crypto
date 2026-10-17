# Experiment log: a cheaper `k = 4` decomposition over `F_{p⁴}`

This lists every run made for this study, in order: pilots, scratch runs,
failures, the claim runs and the audit. It is append-only. Times are UTC.
Host: Claude Code cloud container, Xeon @ 2.10 GHz, 4 vCPU.

Statuses:
- **pilot**: before `PROTOCOL.md`; it chose options and predictions, and
  no claim rests on it.
- **claim**: registered beforehand.

"`p = 269` planted" means `examples/gaudry_quartic_d4 --p 269 --seed 1
--planted` with the flags shown. Pilot rows are in `pilot/`.

| # | when (2026-10-14) | status | what | flags / spec | outcome | data |
|--:|:--|:--|:--|:--|:--|:--|
| O1 | ~18:40 | pilot | Weighted form (σ₄ = π², weights 2,2,2,1), `p = 269` planted, 3 targets | `--weighted` (degree 26) | 3,979 × 4,200, dimension 512, `C₄` 1.49·10⁹ (against 2.55·10¹⁰), exact 3/3 | not kept |
| O2 | ~18:45 | pilot | Weighted degrees 22–25 and 27, 2 planted each | `--weighted --macaulay D` | 22–25 **do not close** (dimension 502, 506, 510, 511). 27 closes at 1.72·10⁹. 26 is the minimum. Phase split at 26: elimination 32 %, normal forms 50 %, charpoly 9 %, eigenvectors 9 % | `pilot/d26.jsonl`, `pilot/d27.jsonl` |
| O3 | ~18:50 | pilot | Border count of the multiplication matrix (debug print) | `--weighted --debug`, plain `--debug` | Weighted `π`: **44 of 512** products leave the standard monomials. Plain `σ₁`: 260. This led to the chain step | not kept |
| O4 | ~19:00 | pilot | Chain step, `p = 269` planted, 4 targets | `--weighted --chains` | `C₄` 5.09·10⁸, normal forms 1.7·10⁷, det A 1.3·10⁷, exact 4/4 | `pilot/ch.jsonl` |
| O5 | ~19:10 | pilot | + sparsest pivots, 3 planted | `--weighted --chains --sparsest` | `C₄` 1.41·10⁸ (elimination 4.8·10⁸ → 1.2·10⁸), exact | not kept |
| O6 | ~19:20 | pilot | All 8 combinations of weighted/chains/sparsest, 2 planted each | see the ablation in `README.md` | Plain + sparsest alone: 2.55·10¹⁰ → 1.0·10⁹. Chains on the plain form: det A dominates (260 chains). All exact | not kept (re-run as the registered ablation) |
| O7 | ~19:30 | pilot | + shape replay, 4 planted | `… --replay` | `C₄` 6.9·10⁷ replayed (elimination 1.2·10⁸ → 4.9·10⁷), exact | not kept |
| O8 | ~19:35 | pilot | `p = 613, 1009` (interpolated det A), 2 planted; weighted dense, 1 target | `… --replay` and `--weighted --sparsest` | `C₄` 8.2·10⁷ replayed against dense 5.3·10⁸, the same solutions | not kept |
| O9 | ~19:40 | pilot | `p = 29` (ecbench curve), 300 uniform, replay **without top-up** | `--p 29 --find 16 --targets 300 --weighted --chains --sparsest --replay` | **42 % of replays fell back** (124/297); 3 did not close (dimension 509–511); exact 297/297 | `pilot/p29_new.jsonl`, `.log` |
| O10 | ~19:40 | pilot | `p = 29`, committed oracle, same 300 targets | `--p 29 --find 16 --targets 300` | All 300 closed; the same solution sets on the 297 both closed | `pilot/p29_plain.jsonl`, `.log` |
| O11 | ~19:50 | pilot | Why replays fell back (reason counter), 40 targets | as O9 | All 12 fallbacks: "pivots or dimension differ" (rank shortfall). This led to the top-up | not kept |
| O12 | ~19:55 | pilot | Replay **with top-up**, 300 uniform | as O9 | 298 replays, 1 fallback (leads differ); mean `C₄` 6.23·10⁷; exact 297/297 | `pilot/p29_new2.jsonl`, `.log` |
| O13 | ~20:00 | pilot | Regression: the plain path against the committed counts, `p = 269` planted, targets 0–2 | no flags | **Identical** in every field to `../gaudry_quartic_d4_20261013/raw/solves_planted.jsonl` | `pilot/regress.jsonl` |
| O14 | ~20:00 | pilot | ecbench plumbing: rho + the new oracle, 20 trials | `pilot/scratch-spec.json`, `ECBS1h0921bb628aee` | rho verified; IC exhausted at 20 trials (as intended); 1.23·10⁹ ops, 19 replays | `pilot/scratch-session/` |
| O15 | 18:30–20:00 | — | Unit tests | `cargo test --release --lib gaudry_quartic` | 12/12 passed at `2b49f38d1` | — |
| O16 | 20:06 | — | **Protocol registered** | `PROTOCOL.md`, spec | commit `2be9ad3d7` | — |
| O17 | 20:08–20:35 | claim | Part A ladder and ablation (22 configurations) | `run_ladder.sh` (3 threads) | All held (A1–A4) | `raw/` |
| O18 | 20:07 to 00:45 (10-15) | claim | Session `k4-p29-d4cost` (`ECBS1h4687f1d786c2`), 8 runs | `specs/SPEC-k4-p29-d4cost.json` | 8/8 verified; new/committed 398× and 406× | `sessions/k4-p29-d4cost/` |
| O19 | 10-15 ~00:45–02:17 | claim | Replay audit | `ecbench verify --replay-all` | OK, 8/8 replays reproduced | `audits/` |
| O20 | 10-15 ~00:50 | — | `ecbench compare` (3 comparisons); one, committed against rho A/A, was **run by mistake and deleted**; `ecbench table` | — | 2 comparisons kept | `sessions/…/comparisons/`, `raw/table.md` |

Not experiments, but logged: two shell commands killed themselves while
waiting (`pkill`/`pgrep` matched their own command line). The ladder they
waited on was not affected.
