# Budget-qualified counted SAT/PDP panel

Frozen 2026-10-09 after the `counted-pilot` audit and before this panel.
The pilot's two measured `ic-sat-cnf` runs verified but recorded 56
per-call `solver_budget_exceeded` events in total. They violated the
`PROTOCOL-2.md` gate, so its six-arm main panel is inadmissible. That
panel was accidentally started; the 208-record stale attempt and the
80-record cleanly interrupted attempt are retained separately and none of
their records enters this comparison.

This follow-on removes **only** the budget-hitting CNF arm. The remaining
five arms had verified answers, no trial or solver budget exit, no timeout,
and identical counted replays in the pilot. This selection is based on the
predeclared budget gate, not on their measured cost. No method parameter,
factor base, target law, solver budget, accounting rule or reference changes.

## Hypothesis and fixed inputs

The hypothesis is that `ic.pipeline_counted` with MITM, Frobenius MITM,
and native-XOR SAT can each complete a matched eight-target panel with
bit-exact replays and expose the complete counted and unpriced resource
vector. The full-cost comparison to the strong signed-Frobenius rho is
**undetermined** wherever `cost.lower_bound` is true. A low counted `S`
alone cannot establish an IC speedup.

Use [`specs/count-qualified-main.json`](specs/count-qualified-main.json):
the registered `icv1-f2m17-tm101-00378d4e` Koblitz curve, subgroup
`r = 65587`, eight independent `hash_to_subgroup_v1` public targets from
seed 42018, `koblitz-orbit:divisor=0;1`, five measured rounds after one
warm-up, algorithm seed 9123403, 120-second execution cap and
`max_trials=20000`. The strong rho arm and its A/A duplicate remain the
reference and control. The IC arms are `mitm:m=2`,
`mitm-frobenius:m=2`, and algebraic descent with
`sat-cdcl:sat_conflict_budget=20000,sat_xor_encoding=native`; all use
`incremental-gauss` and the same `walk` target source. The generic
collision floor is `sqrt(pi/68)` in `S`; the matched strong rho is the
measured reference. Every arm uses the same public workloads and seeds.

Stop on any unverified run, wall timeout, `max_trials` exhaustion,
per-call SAT conflict-budget event, or failed deterministic replay. Retain
the failure and do not silently rerun or omit that arm. If all complete,
audit every measured deterministic run, compare only matched
(workload, round) pairs by the ratio of summed cold `S`, and report the
two-stage bootstrap interval, per-arm completion, setup and online phases,
solver calls and conflicts, budget events, memory, isolation levels and
all `*_uncharged` units. Save the sealed session, audit receipt and raw
comparison files. The Mac gives L0 counted evidence; any wall result
requires a separate L2/L3 panel with an A/A floor. Scaling claims require
other curve sizes, including the m = 83 gate for ECC2K-130 transfer.
