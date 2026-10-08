# Amendment 1: rounds are independent seeds, not replays

Status: written after the m = 13 and m = 19 sessions ran, and **before**
the m = 23 and m = 31 sessions ran. The specs, sessions, sizes, caps and
decision rules are unchanged. This amendment fixes one factual error in
`PROTOCOL.md` and the aggregation that error left undefined.

## The error

`PROTOCOL.md` says "Counts are deterministic, so repeats only serve wall
time." That is false for this harness. Each measured round of a workload
carries its own `algorithm_seed`, so the IC arms walk a different relation
stream in each round, and rho walks a different path. Every round is an
independent run on the same public point. What *is* deterministic is a
single run given its seed: `ecbench verify --replay` reproduces it exactly.

Within a (workload, round), all three arms share one seed. That was checked
on both completed sessions: 24 of 24 pairs at m = 13 and 24 of 24 at
m = 19. So the F4 and F6 arms are paired per **(target, round)**.

## The aggregation, fixed now

- **Unit.** The paired unit is (target, round). For each pair the analysis
  checks that the F4 and F6 arms share the seed and issue the same query
  stream: the same trials and the same relations.
- **Per-target value.** `y(t)` is the mean of `log₂(W4/W6)` over target
  `t`'s measured rounds: three at m = 13, 19 and 23, one at m = 31. Every
  other per-target quantity (each arm's `S`, `S/S_rho`, `W`, `G6`,
  relations) is likewise the mean over that target's rounds. Per-size
  figures are medians over targets of these per-target values.
- **Fit and bootstrap.** These are unchanged in form: OLS of `y(t)` on `m`
  over targets, with targets resampled within size, so a target's rounds
  move together.
- **Sensitivity.** The same fit is also reported using each target's
  **first measured round only**, the one-value-per-target reading the
  original text implied. If the two readings lead to different decisions,
  that disagreement is reported. The primary reading decides.

## Disclosure

Before this was written, the analyzer's first draft (which assumed
identical rounds and took the first measured round) printed interim
per-size medians of `log₂(W4/W6)`: 0.377 at m = 13 and 0.434 at m = 19. It
also showed that the abandon rule did not fire: at m = 19, W6 < W4 on every
target. No fit, interval or decision was computed. The choice above follows
the protocol's own per-target wording and uses every run. It was made to
fix the error, not because of those two numbers.
