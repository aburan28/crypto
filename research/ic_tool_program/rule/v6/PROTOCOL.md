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
