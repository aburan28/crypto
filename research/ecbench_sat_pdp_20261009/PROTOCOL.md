# Protocol: matched SAT and PDP choices inside one cold IC pipeline

Preregistered 2026-10-09, before the sessions in this directory.

## Question and boundary

Can ecbench select a SAT engine and distinct point-decomposition
procedures (PDPs) in a complete, one-target IC run, verify their answers,
and preserve the whole cost of each arm? This panel tests the current
plug points and accounting; it does **not** treat a solver-call reduction
or a partial `S` as an attack speedup.

The curve is `icv1-f2m17-tm101-00378d4e`, with subgroup order 65,587.
Its signed-Frobenius generic floor is
`S_floor = sqrt(pi / (2 * 34)) ≈ 0.215`. The same-session strong
single-target rho arm (`rho.signed_frobenius_strong`) is the measured
reference on every workload. All rows use `S = cold total GAE / sqrt(r)`;
the `lower_bound` and `unpriced` fields are part of the result. Memory,
native solver operations, instruction counts where available, and wall
time are separate resource coordinates. A row with unpriced work has
no fully charged `S` speedup.

## Frozen design

| dimension | value |
|:--|:--|
| target law | `hash_to_subgroup_v1` public points, one target per workload |
| factor base | `koblitz-orbit:divisor=0;1` for every IC arm |
| decomposition | `mitm:m=2`, `mitm-frobenius:m=2`, and `descent-algebraic:m=2` |
| algebraic engines | `buchberger-f2`, `sat-cdcl:sat_conflict_budget=20000,sat_xor_encoding=cnf`, and the same CDCL budget with `sat_xor_encoding=native` |
| relation search | `targets=walk`, `max_trials=20000`, `linalg=incremental-gauss` |
| pilot | two independent targets, one measured round, one warm-up, target seed 41, algorithm seed 9123401; per-run timeout 60 seconds |
| main panel | eight independent targets, five measured rounds, one warm-up, target seed 42017, algorithm seed 9123402; per-run timeout 120 seconds |
| host | Mac arm64 for L0 counted-work feasibility; wall time descriptive. A separate L2/L3 host and A/A interval are required for runtime conclusions. |

Specs: [`specs/pilot.json`](specs/pilot.json) and
[`specs/main.json`](specs/main.json). Every arm gets identical workload
ids and the same seed in a round. The rho control repeats the reference
method to measure an A/A spread; it cannot promote L0 wall figures.
The SAT engine returns one checked model per call, whereas Buchberger
can enumerate multiple solutions. That difference is recorded and any
hit-rate change is part of the complete pipeline outcome.

## Pilot stop and main acceptance

- Run `ecbench plan` first. The slug must be registered and every
  resolved method id must be distinct except the intentional rho A/A.
- The pilot succeeds as a feasibility check only if all arms complete
  on both public targets, every returned scalar checks as `[k]G = Q`,
  and no arm exceeds 60 seconds or its deterministic trial/solver budget.
  Retain every error, budget exit and timeout. If it fails, stop the
  main panel and report the exact failed arm and phase; do not relax
  a budget after seeing it.
- If the pilot succeeds, run the main panel with the frozen inputs.
  For any speed or correctness comparison, require all matched
  (workload, round) pairs complete and verify. Quote the ratio of total
  `S` sums to matched rho with a two-stage cluster interval, but label
  it a lower-bound ratio wherever work is unpriced.
- Never infer a complete-method gain from SAT conflicts, solver wall
  time, a PDP hit rate, or relation counts alone. Native work without
  a pinned same-curve conversion remains explicitly unpriced. A
  runtime claim additionally needs paired L2/L3 runs and an A/A
  interval on a dedicated host, with a 95% interval excluding one.

## Reproduction and decision

Use the committed Rust `ecbench` binary to run each spec into a fresh
directory, `verify --replay-all`, `table`, and matched `compare` against
`rho-strong`. Store the raw sessions, audit receipts and comparison
outputs in this directory, including unsuccessful runs. The SQLite
index is rebuilt from those files. A follow-on note will state whether
the current solver/PDP plugs are operational, which native units still
need calibration, and whether a complete comparison is admissible.
