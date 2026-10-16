# SQISign statistical assertion prevents a green full root gate

- **Status:** investigating
- **Kind:** observed validation failure; cause unknown
- **Affected work:** [fixed-modulus source delivery](../research/isogeny_world_frontier_20261016/replay-optimization-v1/custody-approved-20261016/README.md), TASK-20261016-de577e, crypto PR1991
- **Next action owner:** independent source/resource reviewer for attribution; SQISign maintainer for any separate follow-up
- **Last updated:** 2026-10-16

## Problem

The complete merged-source release library suite is not green because its SQISign statistical assertion failed at source 3ce9694072773badcedc815ddbb953339888e433. The original full-suite requirement remains open. Producer evidence delivery is distinct from source admission, merge eligibility and scientific approval.

The existing test draws a random signing instance over its 37-vertex graph, verifies 4000 tampered messages and requires the observed acceptance rate strictly between 1/(37*4) and 4/37. This report changes neither that test nor its threshold.

## Evidence and uncertainty

- **Observed:** the serial full suite completed 4161 passed / 1 failed / 117 ignored / 0 filtered, exit 101. The assertion observed 27 accepted out of 4000, rate 0.00675. Exact commands, timestamps, source/executable identity, output, failure stack location and OS receipts are in [root-gate.json](../research/isogeny_world_frontier_20261016/replay-optimization-v1/custody-approved-20261016/root-gate.json) and its hash-checked root-gate-3ce9694.tar.gz.
- **Observed:** one prospectively bounded whole-module retake on the same executable passed 14 / 0 / 0, with 4265 filtered and exit 0. It is archived beside the full failure and does not replace it. No adaptive repeats occurred.
- **Observed:** SQISign source, its SHA3 source dependency and Cargo.lock are unchanged between frozen source015fa7bb504e and current3ce, and between main parentd22253655d and current3ce. Exact git comparisons are retained in diagnostic.json inside the archive. The SQISign module SHA256 is 7caef45c3cff41e67f431da9179fb189c9a143e86f37c7c6ef253b6a1428be6f. Current root executable SHA256 is d0df72f83a47794a2099b873ea1a047c46b7fe08762bb26200b4d3e82a1d3535.
- **Observed:** earlier source015 full serial suite passed 4148 / 0 / 117, including this assertion. Its historical receipts remain frozen. Earlier parallel JV failures remain separately retained.
- **Derived:** 27/4000 = 0.00675 is below 1/148. This arithmetic explains the assertion outcome.
- **Suspected explanation:** unknown. A later focused pass and unchanged module do not establish the cause, distributional correctness or parallel stability.
- **Missing capability:** a controlled explanation of the failure, including the actual endpoint law and relevant random instance; the failing run did not retain its random seed.
- **Limits:** source admission remains false. These are same-host engineering tests, not independent review or an algorithm performance result. No isogeny board score changes.

## What would unblock it

- [ ] Establish the failing assertion's cause using a separately reviewed bounded diagnostic that preserves random instance identities and checks its endpoint law against actual native verification.
- [ ] Obtain and retain a complete passing root suite on the final proposed source, with required tests unchanged or any correction independently justified and reviewed.

A passing focused retake alone satisfies neither condition. Any SQISign implementation or statistical-test correction belongs in separately scoped work; this fixed-modulus change does not patch around the failure. Ordinary branch protections and required checks remain in force.

## Next discriminating step

- **Action:** independent source/resource review checks the asserted separation between the fixed-modulus implementation and this validation failure. A maintainer can then design a seeded diagnostic for the native SQISign test in a separate task.
- **Inputs and controls:** exact current and frozen source identities, both preserved executions, unchanged module/dependency comparisons, matching executable SHA256; retain the entire failed sample when a diagnostic exists.
- **Execution limits and stop condition:** no additional execution is authorized by this problem record. This checkpoint completed exactly one bounded whole-module retake and stopped.
- **Outcomes:** a specific cause permits a separately reviewed fix; variation without a cause leaves this record open; infrastructure failure is recorded separately. A focused pass never relabels the full failed run.
- **Recheck trigger:** independent review, a separately approved diagnostic or a new complete root gate.

## Resolution and history

Unresolved. 2026-10-16: preserved the full failure, one focused pass and complete approved remote custody. The producer delivery can be archived honestly while source admission remains closed.

The shared problem-index reservation was released by T-3323. A separate T-3364 reservation permits the new index row; no held scope was overridden.
