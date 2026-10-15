# Root shape replay in the inherited-F4 oracle: a recorded negative (0.95–0.98× on every registered ladder)

2026-10-15. Protocol: [`PROTOCOL.md`](PROTOCOL.md), registered (commit
`d0b27a7f4`) before the claim ladders. Every run is in
[`EXPERIMENTS.md`](EXPERIMENTS.md). Report: [`REPORT.pdf`](REPORT.pdf).
**Class: engineering, a stage diagnostic.** This is the inherited-F4
decomposition stage only, not an end-to-end ECDLP measurement (AGENTS.md
§8).

## Question

This is the follow-up the prime-degree survey ranked first: the one solver
lever from the `k = 4` work (PRs #1863, #1883) that could apply to binary
index calculus over `F_{2^n}`, `n` prime. F4Remake's idea, for the inherited
engine's roots:
- record which Macaulay rows are dependent on earlier rows on the first
  root of a shape (column layout, row count, degree, multiplier policy);
- drop those rows from later roots of the same shape.

A replayed root spans a subspace of the from-scratch one, so its
consequences are sound and the splitting tree still decides every target.
Recording is charged.

## Result

`examples/groebner_stage_bench`, both arms on one binary
(`be69c5ca…`, code `32bf9d358`). Unit: 64-bit word operations in the
eliminations, specialisation included.

| ladder | cells run | default | replay | default ÷ replay | verdicts |
|:--|--:|--:|--:|--:|:--|
| frozen (**pilot**, before the protocol) | 6 | 45,517,794 | 44,815,209 | 1.016 | identical |
| `chain` (registered) | 3 | 7,103,844 | 7,322,495 | **0.970** | identical |
| `chain-holdout` (registered) | 4 | 43,501,345 | 45,579,107 | **0.954** | identical |
| `holdout` (registered) | 2 | 41,414,967 | 42,137,329 | **0.983** | identical |

**Replay costs more on every one of the 9 registered cells: 0.753–0.985×.**
It splits more on 6 of them (e.g. 3,512 against 3,489 splits on
`K_1/2^15`, `m = 4`). Verdicts, decomposed and exhausted counts are
identical everywhere. The harness skips cofactor-inadmissible cells before
measuring them: `chain` ran 3 of its 5 listed cells, `chain-holdout` 4 of
11.

![default ÷ replay by cell](figures/ratio_by_cell.svg)

## Why (diagnostic D1, after the claim runs)

This was a counting build: counters only, with word operations and verdicts
identical to the claim run.

| cell | roots recorded | roots replayed | rows dropped per replay | ratio |
|:--|--:|--:|--:|--:|
| `K_1/2^9`, `m = 4` | 11 | 17 | 5.1 | 0.929 |
| `K_0/2^15`, `m = 4` | 18 | 6 | 36.3 | 0.753 |
| `K_1/2^15`, `m = 4` | 81 | 1,832 | 3.8 | 0.985 |
| `K_0/2^13`, `m = 3`, fresh targets | 3 | 13 | 13.9 | 0.962 |

- The large, rank-deficient roots (where `RESEARCH_INHERITED_F4.md` §3.4 saw
  ~40 % deficiency) rarely repeat their exact column layout across targets,
  so they are recorded and not reused.
- The roots that do repeat lose only ~4 rows.
- The subspace root then reads a slightly weaker tail and splits more.

The survey's ceiling (≈ 10 % of the oracle on such cells) is not reached
on any cell.

## Predictions (PROTOCOL.md)

- **R1 (exact): held**, 9/9 cells.
- **R2 (ratio in `[0.95, 1.10]` per ladder): held**, at 0.970, 0.954 and
  0.983.
- **R3 (no cell gains over 1.25×): held**; none gains at all.
- **Decision:** below 1.05 on every ladder, so lead 1 is a **recorded
  null** (a loss). `KIC_F4_ROOT_REPLAY` stays off and is kept as a control;
  the oracle is unchanged.

## Scope

- The inherited-F4 stage only. The algebraic oracle it belongs to is already
  orders of magnitude behind the root-finding table oracle
  (`research/pdp_root_finder_20261012`: 4.1·10⁷ word XORs against 592
  multiplications a target at `(23, 8)`), and that oracle still loses to
  rho on every measured size. No scoreboard row moves.
- The prime-degree survey's other leads (a weighted/block grading of the
  full symmetrised model; `decompose_by_e3` as an ecbench engine) are
  separate studies.

## Reproduction

```sh
cargo build --release --example groebner_stage_bench --example pdp_shape_replay_summary
research/pdp_shape_replay_20261015/run.sh
target/release/examples/pdp_shape_replay_summary research/pdp_shape_replay_20261015
python3 research/pdp_shape_replay_20261015/figures/plot.py   # layout only
```

Environment: Claude Code cloud container, Xeon @ 2.10 GHz, 4 vCPU, shared
with the characteristic-2 ecbench sessions of #1883. Counts are
deterministic, and wall time is not used. Conductor task T-3.
