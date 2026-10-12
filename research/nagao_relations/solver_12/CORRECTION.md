# Correction: the d = 7 operation ratio in `summary.json` is not a comparison

**Raised by the author of the round, against its own result, before it was
written up.**  `raw.jsonl` is untouched; this supersedes the reading of one
derived column, per §2 of `AGENTS.md`.

## What is wrong

`summary.json` reports `pair_over_image`, the ratio of counted field API
operations, by summing the pair oracle's operations over an instance and
dividing by the operations `solver_10`'s `quadratic-image` run spent on the
same instance.  For `d = 6` that is a comparison of two **completed**
exhaustions and it stands.  For `d = 7` it is not:

| `d` | pair oracle | quadratic-image (solver_10, `enumerate`) | ratio as printed |
|--:|---|---|--:|
| 6 | complete, 8/8 | **complete, 8/8** | `0.494` / `0.496` — valid |
| 7 | complete, 8/8 | **timeout, 0/8** | `0.769` / `0.774` — **not a comparison** |

Every `d = 7` `quadratic-image` enumeration in `solver_10` hit the 120 s
all-phase budget and stopped early.  Its operation count is therefore the work
it managed before the deadline, not the work the instance costs.  Dividing a
complete run by a truncated one understates the true ratio by an unknown
factor, and `research/nagao_relations/scaling_23_29.md` already names this
exact error: *"comparing a full run with a timeout would change the amount of
work being solved"*.

## What survives

The `d = 6` result is unaffected and is the one the round's falsifier turns on:

> **Brute-force pair enumeration exhausts the same instance in `0.49×` the
> counted field operations of the `quadratic-image` solver** — about twice as
> cheap — while returning the identical solution set on all eight instances.

The cross-oracle equality check is also unaffected, because it compares
solution sets and was only applied where `solver_10` reported `complete`; the
`d = 7` instances have `cross_oracle: null` in `raw.jsonl` precisely because
no completed run existed to compare against.

## What replaces the d = 7 row

`solver_13` runs both oracles to **completion** at `d = 4, 5, 6, 7, 8` with
budgets sized per dimension, which is the measurement the `d = 7` row was
reaching for and could not make.  Until then, there is no valid `d = 7`
operation ratio in this thread and none should be quoted.

## Why this was not caught by the contract

The contract required cross-oracle equality *"on every `d = 6` instance where
`solver_10`'s `quadratic-image` run reported status `complete`"* — correctly
conditioned — but the `primary_metric` clause said only *"matched against
`solver_10`'s numbers on the same targets"* and did not repeat the condition.
The guard existed for correctness and was missing for cost.  A future round
comparing against stored numbers should filter on `status == 'complete'` in the
aggregation itself, not only in the equality check.
