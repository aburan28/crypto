# The genus-3 IC-2LP study that 04b follows up is not in the repository

- **Status:** `open`
- **Kind:** access/tool limit: the source artefacts are absent, or held outside this repository
- **Affected work:** [`research/ic_improve_20261014/04b-genus3-2lp-rebalance/`](../research/ic_improve_20261014/04b-genus3-2lp-rebalance/README.md) (Conductor T-1)
- **Next action owner:** the requester; whoever holds `04-genus3-2lp`
- **Last updated:** 2026-10-15

## Problem

04b cannot pair its arms against the original IC-2LP arm, on the same
instances, seeds and counted unit, because that arm's study, protocol, rules,
tool and records are not in this repository.

The request named these artefacts:

- `research/ic_improve_20261014/04-genus3-2lp/README.md` and its `PROTOCOL.md`;
- `research/ic_improve_20261014/RULES.md`;
- `examples/ic_improve_04_genus3_2lp.rs`;
- `examples/ic_improve_04_analyze.rs`.

It also quoted the original's findings:

- local exponents of 1.25 to 1.46 between 2^10 and 2^20;
- linear algebra at q^1.52;
- relation phase at q^1.30;
- row weight from 7 to 25;
- F = ceil(2 q^{2/3}), tuned at 2^14.

At the requester's choice, 04b rebuilt the baseline instead. Its arms are
paired against that rebuild, not against the original.

## Evidence and uncertainty

- **Observed.** The artefacts were searched for on 2026-10-15 and found nowhere:
  - on `main` (at `0dc9ef4a7`);
  - on the working branch;
  - on every remote branch named with 20261014 or 20261015, and `codex/ic-improvement-round3`;
  - in GitHub code search for `ic_improve_04`;
  - in pull-request searches for `ic_improve_20261014` and genus-3 large-prime/GHS work;
  - among the 30 most recent pull requests.
- **Derived or proved:** not applicable.
- **Suspected explanation:** unknown. Candidates are an unpushed branch in
  another session, another repository, or a mistaken path.
- **Missing capability:** access to the original records and tool.
- **What this evidence does not establish.**
  - It does not show that the original numbers are wrong.
  - 04b's rebuilt baseline differs in constants and probably in design detail:
    its row weight goes from 4.9 to 14.4 rather than 7 to 25, and its
    linear-algebra exponent is 1.447 rather than 1.52.
  - So 04b's verdicts are about the rebuild only.

## What would unblock it

- [ ] The original tool and records are committed, or linked from another
  repository, at a pinned revision.
- [ ] The two baselines are compared on matched sizes in one counted unit:
  local exponents, phase exponents and row weights. Any design difference is
  named. In particular: the large-prime elimination order, the stop margin,
  how the factor base is ordered, and the Lanczos variant.

Meeting both would show whether 04b's answer carries over to the original
baseline. If the designs differ, A1 and A2 would have to be re-run on the
original tool.

## Next discriminating step

- **Action:** obtain the original 04 artefacts from the requester.
- **Inputs and controls:**
  - 04b's committed curves (seed 20261016) and `analysis.json`;
  - the original's instances, if published.
- **Execution limits and stop condition:** the original's own §-rules (60
  minutes per job, 3 GB).
- **Outcomes and next actions:**
  - If the baselines agree in exponents, 04b's verdicts stand for both.
  - If they disagree, run A1 and A2 on the original tool.
- **Recheck trigger:** the artefacts appear.

## Resolution and history

- 2026-10-15: opened with 04b. The rebuilt baseline departs from 4/3:
  E = 1.383, 99% interval [1.382, 1.384]. Neither requested arm removes the
  departure.
