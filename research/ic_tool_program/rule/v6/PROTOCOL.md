# The rule's comparison at baseline v6

**Declared 2026-10-16, in the pull request that accepts v6 (R08), before
any run below.** IC_TOOL_PROGRAM.md asks for the rule's comparison at every
new baseline. v4 (R10), v5 (R06) and v6 (R08) were accepted the same day,
and the comparison runs once, at v6. That baseline carries both of the
programme's levers. Rules at v4 and v5 are not run.

**Everything rule v3 declared holds here, except the rows under "What
changes"** ([`../v3/PROTOCOL.md`](../v3/PROTOCOL.md), and through it rule
v2's). That covers:
- §23's sizes, rows, arms, retries and stop rule: six sizes, 64 targets
  and four A/A repeats a size, 408 rows;
- the claims and the reading, and the inadmissible list;
- the native tools (`icprog rule`, `isolated_bench`);
- `icprog`'s oracle as the replay's checker;
- the pin against §23's `T01-R1`, and each row's logarithms checked
  against §23's row;
- the reference check against the strong rho;
- the accounting, class and non-claims.

## Question

One unseen public point, the index calculus against Pollard rho on that
same point, online and cold, at §23's six sizes. **How far is rho still
ahead at v6?** At v3 rho was 5.2–21.0× faster cold (13.2×, 12.0× and 21.0×
at the three largest sizes). That remaining distance is what CLAUDE.md
asks every result to report.

## What changes from rule v3

- **The binary** is v6's, the one R08 measured: SHA-256
  `4ccb1951134ff31b5183fc00033933b9add2e3c0c52003f9041fd04bac1f4012`, built
  from local commit `3d4c5c0ef88057e8ece09dcc91e8fad4ff7e781a` (v4 plus
  `R06-scan-key/candidate-v4.patch` plus `R08-scan-fold/candidate-v5.patch`,
  tree `ab57f3fc`). The runner is given the commit (`--ic-commit`).
- **The arms' sources are still §23's list.** Between v3 and v6, main's
  drift and the two rounds changed seven of its nine files. The candidate
  IDs are new, as they should be.
- **The strong rho is v6's own.** `koblitz_rho_fixture` is built from the
  same tree, so the rho walk's per-step cost is measured afresh. R06 and
  R08 change only the index calculus's scan and build; rho's online
  interval moved by nothing reliable in either round.
- **What each size is read against:** §23's measurement under the key
  `s23_then_v6`. The README sets rule v3's figure beside each one.
- **The host** is the Emerald Rapids of R10, R06 and R08, kernel build
  `6.18.44-fc-v139`. Rule v3 ran on another host (`6.18.44-fc-v70`). Its
  manifest is `runs/host.json`.

## The pricing band (AGENTS.md §5a)

**Both arms are timed, not counted.** Each row runs the index calculus and
rho in one process on one isolated core, and the cold figure is wall time
over wall time. Every native operation either arm performs is inside its
time: rho's canonicalisations and table work, and the index calculus's
keys, filters and sorts. So no counter is left unpriced on either arm, and
the band at 0, 0.1 and 1 group-addition equivalent per native operation
gives one figure.

**Where counted units do appear, the band applies:**
- the rule's secondary reading, rho priced at the canonical step;
- the Bernstein–Lange model.

Those are reported as rule v3 reported them, as lower-bound readings, and
neither ranks the arms.

## Prediction

Not drawn from a run. If the rho arm's cost is as at v3, a cold ratio
should fall by about what the rounds gained in set-up at that size. At the
three largest sizes that is about 1.6–1.7× (R06 × R08, which also include
v4's ~1%), so about 8×, 7× and 12×. At the three smaller sizes it is less.
The host changed, and main's drift has moved rho's online interval before
(R07, v2 to v3: 1.11–1.82× faster), so this is a direction, not a band.

## Success and stop

The rule makes no acceptance decision. It is a measurement.

**Stopped** if:
- a row fails its check;
- a logarithm differs from §23's;
- the pin fails;
- the reference check finds the comparison's rho slower a step than the
  strong walk. That makes the rho arm inadmissible, and the comparison is
  reported as not run.

## Cost, as an estimate

v6's set-up is 1.6–1.7× faster than v4's at the top sizes, so the 408 rows
take under an hour, and the reference check's 48 walks take minutes.

## Commands

    icprog rule manifest --comparison v6 --ic <v6 ic> --ic-commit 3d4c5c0ef88057e8ece09dcc91e8fad4ff7e781a --isolate <isolated_bench>
    icprog rule all --comparison v6 --ic <v6 ic> --ic-commit 3d4c5c0ef88057e8ece09dcc91e8fad4ff7e781a --isolate <isolated_bench>
    icprog rule reference --comparison v6 --fixture <koblitz_rho_fixture> --fixture-commit 3d4c5c0ef88057e8ece09dcc91e8fad4ff7e781a --isolate <isolated_bench>
    icprog rule claims --comparison v6
    icprog rule analyse --comparison v6 > research/ic_tool_program/rule/v6/analysis.json

## Amendment 1 (2026-10-16, after the rows): the claims and the reference walk

The 408 rows ran as declared, 15:40–16:56 UTC. While diagnosing what
follows, one figure was seen: the cold ratio at `2^36.6`, 3.78. Two tools
then failed for
reasons that lie outside the measurement. Both are fixed without
re-timing any row.

**The claims.** `icprog rule claims` wrote claims that the current
`vs_rho` checker rejected: 0 of 408 passed.
- No row had a validation error, and every row passed independent
  validation.
- Each claim lacked the same five stage fields: `rho_online_phase_ms`,
  `ic_online_interval`, `rho_online_interval`, `operation_accounting` and
  `record_class`.
- The ledger (`docs/ic/boundary_targets.json`) began requiring them on
  2026-10-13 (`e2b865e4`), a week after rule v3 ran. The rule's claim
  writer was never updated.

`icprog` now writes all five fields from what each run record already
holds:
- **Rho's phases:** `rho_solve` and `recovery_check`, in milliseconds.
- **The intervals:** each arm's recorded start and stop events.
- **The operation accounting:** each arm's native counter, kept distinct
  (the descent's trials and rho's walk steps), with no operation ratio.
- **The record class:** `controlled` when all five of the index
  calculus's exclusive target phases have a cost, `exploratory`
  otherwise.

With the rewritten claims, 408 of 408 rows pass, all of them
`verified_answer_controlled_wall`. The rows and their timings are
unchanged.

**The reference walk.** `koblitz_rho_fixture strong`, built from v6's
tree, refused `n = 47`, `57` and `59` with
`assert!(matches!(n, 7 | 11 | 13 | 17 | 19 | 23 | 37 | 41 | 53 | 61))`.
- The gate is not in v3's fixture, whose walk the rule v3 reference ran
  at all six sizes.
- It came back with a union-resolve merge (`fe07f6088`) from a branch
  line (`c72806c87`).
- No code after it depends on `n`. The single-word walk is v3's, and
  wider rungs return before the gate.

The gate now reads `2 <= n <= 63`. The walk is otherwise byte for byte
v6's: local commit `e095054a` is v6's tree plus that one change, and the
fixture's SHA-256 is
`5259f49cdf188adb4fc16749ec49f40f8fdd3edae3853285720aed390f82b5a4`.

**The reference check re-runs at all six sizes** with that one build:
- The refused attempts stay on record under
  `runs/reference-refused-size-gate/`.
- In those attempts the three sizes the gate let through were admissible:
  rho at 0.65–0.75 of the strong walk's time a step. They are not pooled
  with the re-run.
- If the re-run finds the comparison's rho slower a step than the strong
  walk at any size, that size is reported as not run, as declared.
