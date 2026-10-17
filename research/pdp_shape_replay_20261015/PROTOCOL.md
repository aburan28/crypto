# Protocol: cross-target root shape replay in the inherited-F4 oracle (prime-degree binary IC)

Registered 2026-10-15, before the claim runs below, after the one pilot
listed in `EXPERIMENTS.md` (P3, the frozen ladder). Code: commit `32bf9d358`
on branch `claude/pdp-shape-replay-20261015`, off `main` at `e15e29dbd`;
`examples/groebner_stage_bench` sha256 `be69c5cac3bdbfbd…`. The switch is
`KIC_F4_ROOT_REPLAY=1`; it is off by default, so existing counts are
unchanged. Class: **engineering** (a stage constant at most).

## Question

For Koblitz curves over `F_{2^n}`, `n` prime, the inherited-F4 decomposition
oracle reduces a degree-`D` Macaulay matrix at each target's root, then
specialises it down a splitting tree. `RESEARCH_INHERITED_F4.md` §3.3–3.4
puts roots at 27 % of the oracle on its dominant rung, with the cubic
`m = 3` roots about 40 % rank-deficient. The question: does recording a
root's dependent rows on one target, and dropping them on later targets of
the same shape (F4Remake), lower the stage's word operations without
changing any verdict?

## Arms

The same binary, the same ladders:
- default;
- `KIC_F4_ROOT_REPLAY=1`.

Replay drops rows, so a replayed root spans a subspace of the default's row
space. Its consequences are sound, and the tree still decides every target;
the cost of a weaker root is in the counts. Recording costs one incremental
basis per shape and is charged.

## Runs (`examples/groebner_stage_bench`)

Ladders the pilot did not use:
- `--ladder chain`: `m ≥ 3`, 5 cells, the rank-deficient cubic roots;
- `--ladder chain-holdout`: 11 cells, `m = 3, 4`;
- `--ladder holdout`: `m = 2`, 2 cells.

Each ladder runs once per arm. Counts are deterministic, so one run each
suffices. Outputs go to `raw/<ladder>-<arm>/stage.json`.

## Predictions (from the pilot's total of 1.016× on the frozen ladder)

- **R1 (exact).** Every cell's verdict digest, and its decomposed and
  exhausted counts, are identical between the arms.
- **R2.** Default ÷ replay total `word_ops` lies in `[0.95, 1.10]` on every
  ladder.
- **R3.** No cell gains more than 1.25×.

## Decision

- If R2's ratio is below 1.05 on every ladder, lead 1 of the prime-degree
  survey is a recorded null. The switch stays off and the oracle is
  unchanged.
- At 1.05 or above on `chain` and `chain-holdout`, it is a constant on the
  algebraic oracle. That oracle stays orders of magnitude behind the
  root-finding table oracle (`research/pdp_root_finder_20261012`), so no IC
  scoreboard row moves.
- A miss on R1 retracts the change.
