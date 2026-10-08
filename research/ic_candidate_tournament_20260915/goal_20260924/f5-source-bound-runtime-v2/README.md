# F5 v2: complete source-bound development solve

The single invocation of the [frozen protocol](PROTOCOL.md) is closed.
[PR #1060](https://github.com/aburan28/crypto/pull/1060) accepted the controller at
`8b49542e7e7aa50f65edc8831faabd0b57032acb`. [PR #1064](https://github.com/aburan28/crypto/pull/1064)
committed the exact [preexecution receipt](PREEXECUTION.json) at
`f1c5e41bebe73a008d8861a2bdee23f2695fd93a` before the only invocation; the
[dispatch observation](DISPATCH-OBSERVATION.json) followed at
`51d13a1187b04eb92090ce649f96c743e642a5e2`. The historical pending states in
PROTOCOL and PREEXECUTION are preserved. [TERMINAL.json](TERMINAL.json) is current.
No retry, resume or budget extension is permitted.

The physical macOS ARM64 development job used exposed public point
`[52411,72106]`, algorithm seed 2026093032, one worker, degree-three bounded
MatrixF5 with an 8192-node PDP budget, dense final relation LA, 512 ordinary/target
attempt caps and 7200/7170-second controller/native watchdogs. Input seed list
was empty; fixture provenance was null. No scalar was supplied. Native source,
binary, stdin, interpreter and complete Python pre/post gates match the frozen
invocation. Both controller and native process returned zero, without timeout.
The [native interface prerequisite](interface-control/README.md) remains a
separate parser-only control.

## Mathematical admission

| Candidate / workload / run | Actual usable base / folded columns | Result | Native online ns | Rho / IC |
| --- | --- | --- | --- | --- |
| `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h84b504d19844` / `12f20a62c94c` / R0 | 62 / 29 | Verified one-target IC; rank 29/29 | 10,512,454,542 | Unknown |

The independent checker reconstructed the exact geometric base, cofactor
projection, sign/Frobenius quotient coefficients, scalar-prime matrix, full
rank, all 29 column logs and the target relation. It verified 61 relations in
216 ordinary queries, zero duplicates and 32 rows beyond rank. The target
needed three attempts: two independently proved group-infeasible queries,
then a valid decomposition giving scalar 24886; `[24886]G=Q` replay passed.
All unsuccessful attempts remain charged. This is a complete bounded IC path
on this instance; the engine is bounded Macaulay elimination with an F5 row
criterion, not a full incremental Gröbner-basis implementation.

The ordinary census contains 155 `proved_unsat` and 61 verified witnesses.
Exact postexecution three-sum enumeration found exactly 61 feasible ordinary
queries and zero feasible misses on this stream. The measured witness fraction
is 61/216. Its stored nominal Wilson interval is [0.226581, 0.345838]; collection
stops by rank, so this is a descriptive interval without a fixed-sample coverage
guarantee. No planted decomposition or oracle witness entered the native job.
This census does not establish completeness or general F5 soundness beyond the
observed queries.

## Retained audit rejection and explicit negative proof

The original v2 auditor rejected the report with `negative claim needs an
independent proof adapter`: its generic query gate only independently checks
two-summand tiny-field negatives. That gate and every historical evaluator are
unchanged. The first v3 proof prototype also rejected: it wrongly required the
geometric PDP input points to be in the prime subgroup. These two postexecution
audit failures are retained inside the archive, with no retroactive complete
execution-source attestation claimed for those failed audit commands.

The reviewed versioned proof adapter checks the **63 geometric PDP points**;
cofactor projection yields the separate **62 usable IC points**. It keeps torsion
and coset points, repeated summand indices and identity pair sums. It constructs
all pairs and checks Q-P against the complete pair set for each base point.
Every triple appears in this check, so absence proves group infeasibility of
that recorded query. All 157 negative claims pass: 155 ordinary and two target.
Controls reject forged negatives on both an actual ordinary witness and the
actual target witness, oversized domains, and removal of repeated/torsion cases.
The adapter is bounded to three summands, field degree at most 17 and at most
128 geometric points; unsupported domains fail closed.

The new independent analysis uses `audit_f5_runtime_v3.py` and versioned exact
query-law/stage/admission wrappers. It changes no executed candidate source,
algorithm setting, invocation identity, query or cost. The archive separately
freezes a complete postexecution auditor context. Replay imports that retained
context under isolated, site-disabled Python and checks every loaded local
module before/after analysis. It executes no native worker.

## Timing and scope

The native online interval starts at first target-dependent native IC work after
reusable preparation and ends after general-group scalar replay. External
Python admission is excluded. Integer phase costs close exactly:

| Exclusive online phase | ns |
| --- | ---: |
| Target query generation | 17,000 |
| Target PDP, all three attempts | 10,512,420,582 |
| Target relation check | 168 |
| Target descent/recovery arithmetic | 10,292 |
| Scalar replay | 6,500 |
| Total | 10,512,454,542 |

The supplementary native cold ledger closes at 918,648,281,209 ns. Whole
controller wall is 920,906,322,416 ns; controller/source startup and ending gates
are operational diagnostics outside the native arithmetic ledger. Native
per-child RSS is 75,186,176 bytes, not a measured process-tree peak. Raw resource
seconds are retained unchanged and published as roundtrip decimal strings.
This host is uncalibrated, the point is exposed, and no incumbent or rho arm was
registered. Speedup, normalized S/floor, fresh-target qualification and promotion
remain null/false. The F5 and SAT development results are not a paired ranking.
The three sealed confirmation sets and every historical consumed registration
remain closed. The active goal still requires the reviewed exposure/reference/
calibration adapter and a new frozen fresh paired comparison.

## Durable evidence

[results-20260930/evidence.tar.gz](results-20260930/evidence.tar.gz) has SHA-256
`d62ff4b0e5848b7e36474d1a6af2be7fef1ee060d6cc05ae1e02c86c2957c1b2`,
38,314,958 bytes. It retains registration, execution, complete source/native/
interpreter assets, raw stdout/stderr/resource receipts, the two audit rejections,
mathematical admission and separate auditor source context. The external
invocation SHA-256 is
`9409cb6c40548816153f5b9978e8fa0637b5aaa2801472575e7cc55a9f61d56d`.
[receipt.json](results-20260930/receipt.json) inventories member bytes, hashes and
modes. [AUDIT.json](results-20260930/AUDIT.json) contains all query/rank/log/cost
proofs. [TRANSPORT-REPLAY.json](results-20260930/TRANSPORT-REPLAY.json) and
[TRANSPORT-SOURCE-GATES.json](results-20260930/TRANSPORT-SOURCE-GATES.json) retain
passing isolated archive/source/mathematical replay.

[FINAL-FOCUSED.log](results-20260930/FINAL-FOCUSED.log) records ten passing new
controls. [FINAL-FULL.log](results-20260930/FINAL-FULL.log) records 380 regressions,
377 passing and three declared macOS skips (strict frozen Linux statistics and
two Linux affinity/instruction controls), in 200.979 seconds. Skill validation
passes. Exact-head hosted CI remains a separate acceptance gate.

Replay from this repository to a new empty destination:

```sh
python3.12 research/ic_candidate_tournament_20260915/publish_f5_runtime_v3.py replay \
  --bundle research/ic_candidate_tournament_20260915/goal_20260924/f5-source-bound-runtime-v2/results-20260930 \
  --expected-execution-sha256 9409cb6c40548816153f5b9978e8fa0637b5aaa2801472575e7cc55a9f61d56d \
  --out /tmp/ic-f5-v2-transport-NEW
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 \
  -p test_f5_v2_full_development.py -v
```
