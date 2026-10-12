# Follow-on protocol: replayable counted SAT/PDP panel

Frozen 2026-10-09 after the first pilot, before either session below.
The original [`PROTOCOL.md`](PROTOCOL.md) and its failed pilot remain
unchanged. That pilot had 24 verified records, four 60-second
`ic-buchberger` timeouts, and a failed audit: one SAT native-XOR replay
matched its integer work but differed in the low bits of `total_gae`.
Those results stopped the original main panel as preregistered.

This is a **new, narrower comparison**. It removes the timed-out
Buchberger arm; it does not increase its budget or recast its timeout.
All IC arms use the new `ic.pipeline_counted` method id, which re-prices
the relation phase from deterministic group and pinned native counts
instead of cancelling a measured solver-wall charge. Unpriced SAT
conflicts remain explicit lower-bound work. The old `ic.pipeline`
identity and pilot records are preserved for audit.

## Frozen inputs, reference and decision

- Curve: `icv1-f2m17-tm101-00378d4e`, `r = 65587`, signed-Frobenius
  floor `S_floor = sqrt(pi/68) ≈ 0.215`; matched strong rho is the
  measured reference. Each workload is one independent public target
  from `hash_to_subgroup_v1`.
- Factor base: `koblitz-orbit:divisor=0;1`. Relation source `walk`,
  `max_trials=20000`, linear algebra `incremental-gauss` on every IC
  arm. PDPs: `mitm:m=2`, `mitm-frobenius:m=2`, and algebraic descent
  with SAT CNF and SAT native-XOR, both at 20,000 conflicts per call.
- Pilot: [`specs/counted-pilot.json`](specs/counted-pilot.json), two
  targets, target seed 41, one warm-up and one measured round,
  algorithm seed 9123401, 60-second execution timeout. It repeats
  known inputs to test the new accounting and replay contract.
- Main: [`specs/counted-main.json`](specs/counted-main.json), eight
  independent holdout targets, target seed 42018, one warm-up and
  five measured rounds, algorithm seed 9123403, 120-second timeout.
- Reference and control: two arms of
  `rho.signed_frobenius_strong`, on the same targets and round seeds.
  On this Mac all runs are L0; `S` and integer counters are primary,
  wall time is descriptive. An L2/L3 host and A/A interval remain
  necessary for a runtime conclusion.

Run the main panel only if every pilot run verifies, every measured
deterministic run replays identically, and no execution hits the
trial, solver or wall budget. Retain all failures. In the main panel,
compare only complete matched pairs by ratio of total `S` sums with
a two-stage workload/round interval. Mark a ratio involving an IC
row `bounded` if any native work is unpriced; no end-to-end speedup
or cross-size scaling statement follows from such a ratio.

The result must report per-arm completion, cold `S`, setup and online
phases, solver calls/conflicts/budget exits, memory, isolation levels,
and every `*_uncharged` unit. Save source/spec hashes, session files,
full replay receipts, and comparison output. A method that fails its
own check is a failed row, never removed from the table.
