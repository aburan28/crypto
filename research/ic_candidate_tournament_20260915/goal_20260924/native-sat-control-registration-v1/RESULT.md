# Native prepared SAT control completes and independently replays

The sole invocation preregistered in [PR #1274](https://github.com/aburan28/crypto/pull/1274)
is **consumed and closed**. The original frozen native checker admits a complete
source-bound disclosed-input control: scalar `24886` was recovered for
`[52411,72106]` after two proved negative queries and one source-verified point
witness. Do not retry, resume, extend or execute an archive restoration.

PR #1274 was merged at `2026-10-03T03:52:44Z`, commit
`f2a4698f4e913de5e6c4468f639e45a35b8356b0`. All ten applicable checks on accepted
head `9d469d658da8039f2ccf60e24dc23e60fa143da4` passed before dispatch; the last
integration check completed after the owner merged the registration. Linux and
macOS custody receipts are byte-identical to the original committed receipt,
SHA-256 `5de0de203d089df32633cd88163b77114502de030f286f3dc6f14dd9f9abc3a5`.
Execution began only after that final gate passed.

## What executed

The [claim](result-v1/data/consumed.json) was created before launch and fixes
registration `b098ebbe135abd4706e62d6286f513fe1c4f07fad61f7250487bd2d86ad317e4`.
The [producer](result-v1/data/execution/producer.json),
[terminal receipt](result-v1/data/execution/terminal.json), all exporter/CMS
launches, ANF/CNF/Magma inputs, complete solver stdout/stderr, source-model
digests, failed progress records and PID ledgers remain unchanged in
[the retained data tree](result-v1/data).
There were three exporter invocations and three one-thread CryptoMiniSat
invocations. No further scientific invocation or ordinary relation query ran.

| Trial (zero based) | Coefficients `(a,b)` | Public query point | Result | Witness |
| --- | --- | --- | --- | --- |
| 0 | `(60182,25978)` | `[40991,73355]` | Source UNSAT; independently proved no geometric three-sum | None |
| 1 | `(4091,5477)` | `[73003,104622]` | Source UNSAT; independently proved no geometric three-sum | None |
| 2 | `(12679,5261)` | `[59775,2910]` | Valid source model and group relation; scalar recovered | Geometric indices `[53,1,48]` |

The frozen [independent audit](result-v1/data/independent-audit.json) has SHA-256
`3570722341b8ca2cd8794d51829cac1a1f3b8c1ff02a07a6af9377a36d2cbab7`.
It reports `PASS_NATIVE_SOURCE_BOUND_CONTROL_AUDIT`, source-bound execution
admitted, target complete, scalar verified, and zero native children executed
by the auditor. It independently reconstructs the query law, preparation,
source model and actual group witness, cofactor/orbit recovery and `[d]G = Q`.
All source and child-drain gates pass. The query loop stopped at first
success within the frozen eight-query, one-million-conflict-per-query limits.

The successful 51-bit source model has SHA-256
`37d86697adf0f9dc51a527784039f44aa4b6c607599c8a4f3c9dba8972c9f973`;
the full expanded model has SHA-256
`d7271754b64b0ff72211cf53b2cabe41629721948cd3f93108c3ef7514ce408b`.
The retained source equations and independent group replay connect that model
to the actual decomposition; satisfiability alone is not the admission gate.

## Timing and scope

These are uncalibrated physical macOS ARM64 control observations on an Apple
M4 Pro, 14 logical CPUs, 48 GiB memory, macOS 26.6 / Darwin 25.6.0. The
[host context](result-v1/data/host-context.txt) records the original facts;
the busy lock provides serialization, not L2 timing isolation. One solver
thread was used; the memory ceiling and memory peak are unknown.

| Exclusive target-dependent phase | Original nanoseconds |
| --- | ---: |
| Query generation | 11,000 |
| PDP, including both unsuccessful attempts and all native child work | 56,585,062,375 |
| Relation/source-model verification | 396,917 |
| Target descent | 333 |
| Producer scalar recovery check | 9,542 |
| Sum / attempted and completed control interval | **56,585,480,167** |

The outer controller interval is `62,921,422,708 ns`. The five exclusive phases
sum exactly to the producer interval, and that interval fits inside the outer
controller interval. Independent frozen replay took `6,258,817,750 ns` and is
**outside** the producer interval. Source/loading startup and the external
checker are separate. This timing boundary cannot establish the primary
independently verified online speedup. Operation calibration, A/A noise,
normalized operations, matched rho/incumbent costs and speedup remain unknown.

The point is already disclosed; this control has no fresh target. It imports
the existing preparation: 63 geometric points, 62 usable points, 29 folded
columns and rank 29. Its original 149 ordinary records contain 37 verified
relations, 106 exact negatives and six budget-inconclusive attempts, with eight
dependent relations and no duplicates. All are independently reconstructed
again, retaining the historical Python controller provenance. There are zero
new ordinary queries, so this invocation supplies no new natural-yield estimate.
Fresh qualification, headline eligibility and promotion are false.

## Retained publication and portable replay

[manifest.json](result-v1/manifest.json) binds every raw data file;
[seal.json](result-v1/seal.json) fixes its canonical digest.
The native publisher uses create-only copies and never edits the original
claim, producer, audit or source files. Original operator stdout/stderr and
host context remain present. The complete executable/source/dependency capsule
stays in the immutable before-execution [registration publication](publication).

The portable native replay is a **postexecution publication context**. It
retains the exact original frozen audit by its independently recorded digest;
it neither replaces that auditor nor retrospectively changes its source seal.
It rechecks data hashes, original registration/claim/terminal identities,
preparation, every target query, model/witness/scalar, phase sum and exact
role/helper argv, environment and output hashes. It executes no exporter,
solver or archived Mac binary. The [original local publication replay](result-v1/publication-replay.json)
passes and agrees with the unchanged frozen audit.

```sh
target/debug/icprog sat-control-replay-publication \
  --publication research/ic_candidate_tournament_20260915/goal_20260924/native-sat-control-registration-v1/result-v1 \
  --out /tmp/native-sat-consumed-result-replay-new.json
```

Controls reject changed solver thread policy even after the publication
inventory is resealed, altered CNF, an omitted attempt and attempted fresh
promotion. CI runs the replay and controls on Linux x86-64 and macOS ARM64.
Those passes establish portable data/math replay; actual native solver hardware
admission remains the original macOS ARM64 build only.
The [first broad publication harness failure](first-publication-harness-failure.txt)
records a test namespace collision from compiling the transport module twice.
Publication now imports the existing module; no scientific invocation, raw
result or original frozen audit was changed or repeated.

The full F4/F5-plus-SAT goal remains active. This closes actual native SAT
disclosed-input execution. Native F4/F5 controller/preparation verification,
new ordinary-yield evidence, complete exposure exclusions, source-bound strong
references and a new ecbench paired-target protocol remain required. There is
no fresh tournament result or global-fastest claim, and all three historical
confirmation sets stay closed.
