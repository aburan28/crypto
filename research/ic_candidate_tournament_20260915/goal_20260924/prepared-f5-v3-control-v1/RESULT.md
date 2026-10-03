# Complete disclosed prepared F5 source control — closed

The one permitted execution recovered scalar **24886** for the disclosed public
point `[52411,72106]`. The original source-bound runtime completed, and its own
preexecution-frozen transport independently admitted
`ADMITTED_COMPLETE_PREPARED_F5_CONTROL`. This is a correctness/source milestone,
not solver-family qualification or a speedup. The registration is consumed and
closed. Do not repeat it, including from a restored archive.

The preregistration [PR #1159](https://github.com/aburan28/crypto/pull/1159)
merged at `45b2b9adc3405309ba07d005d9d141f89c4a2afb`, after both applicable
checks succeeded on `b3fe186b22b4ef092d782865018d05161dc975ad`:
[IC integration](https://github.com/aburan28/crypto/actions/runs/36890344161)
and [names](https://github.com/aburan28/crypto/actions/runs/36890344091).
The [protocol](PROTOCOL.md), [panel](panel.json), [original registration](REGISTRATION.json)
and external execution seal are unchanged. The historical before-execution
ledger stays historical; [TERMINAL.json](TERMINAL.json) and the
[exclusive original claim](registration/execution-claim.json) record closure.

## Verified behavior and limits

| This one disclosed control | Recorded outcome |
| --- | --- |
| Candidate | `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h5cc365241aec` |
| Workload | `8b3d9065d872` |
| Run | `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h5cc365241aecW8b3d9065d872R2026100104` |
| Curve / target | Frozen n17 binary curve / `[52411,72106]` |
| Actual usable / geometric / folded base | 62 / 63 / 29 |
| Target attempts | 3: two independently proved negatives, one witness |
| Recovery | Scalar 24886, independent group replay passed |
| Raw online wall, all target attempts | 11,955,969,917 ns, uncalibrated development diagnostic |
| Matched rho / incumbent online wall | Unknown; neither arm ran |
| Online speedup / S / floor ratio | Unknown |
| Fresh-target qualification / promotion | False / false |

The algorithm seed was deliberately selected from the known interface fixture.
Neither its 1/3 witness count nor this one timing estimates natural query yield
or expected target performance. It uses a bounded Macaulay engine with an F5
row criterion, not a full incremental Gröbner-basis implementation. Observed
MatrixF5 dispatch, witness indices `[0,4,50]`, exact target query law, the imported
factor logs and final scalar are independently replayed. The two negative target
points `[117235,33918]` and `[89585,21540]` were checked against the exact finite
three-sum group oracle; this is a bounded n17 absence check, not a general
algebraic soundness theorem. No known scalar was supplied to the worker.

The reusable preparation was checked before online start. It has rank 29/29,
61 verified relation rows from 216 ordinary preparation queries. Those historical
natural-query outcomes and failed attempts remain in their accepted preparation
history. This job executes **zero new ordinary queries**, no new relation-matrix
solve and no new factor-base construction. Its final `LAgauss` identifier refers
to the imported preparation's relation matrix; internal GF(2) Macaulay reduction
belongs to PDP. This control cannot replace a new natural-yield study.

## Timing boundary

Native online start follows input loading, preparation/log verification and the
fixed target-independent symbolic template. Stop follows target scalar recovery
and general-group replay. All failed target attempts are included. Independent
Python admission and archive restoration are outside this interval. The physical
macOS ARM64 process used one worker through the repository busy wrapper; there
was no calibrated matched environment or fixed memory ceiling.

| Exclusive target phase | Native ns |
| --- | ---: |
| Query generation | 11,293 |
| PDP encoding and solving, including failures | 11,955,937,581 |
| Relation verification | 251 |
| Target descent/recovery algebra | 11,875 |
| Scalar recovery check | 8,917 |
| **Online sum** | **11,955,969,917** |

The raw native process wall is 11,999,138,125 ns; the controller wall is
14,289,997,333 ns. These have different boundaries and are not added to the
online phases. Native peak RSS was 51,281,920 bytes, a process diagnostic;
it is not an allocation measurement or a memory resource ceiling. The native
and transport before/after source gates agree. [Original native output](controls/pipeline.stdout),
[raw meter](controls/pipeline.metrics.json), [admission](controls/admission.json)
and [transport](controls/transport.json) retain the evidence.

## Durable bytes and independent restoration

The Git-tracked [archive](publication/evidence.tar.gz) retains all 269 original
execution and audit members, complete runtime/assets, stdout/stderr, claim and
source gates. It is **25,476,819 bytes**, SHA-256
`12ea04c458036c0ca56c588d151ae317ada2cbc661124344be472bbecc94804b`.
[receipt.json](publication/receipt.json) inventories every member's bytes, hash
and mode. [pack_evidence.py](pack_evidence.py) captures/restores bytes only and
executes no solver. Exact restoration passed, then the restored execution's own
frozen helper produced a byte-identical mathematical admission and
`PASS_FROZEN_PREPARED_TRANSPORT`; [restoration receipt](controls/restore.stdout)
and [restored transport](controls/restored-transport.json) preserve both checks.
The original interpreter executable and stdlib are hash-bound, not packaged as a
portable installation. Full source/interpreter replay requires that interpreter;
archive integrity can be checked on other hosts without executing a solver.

```sh
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-v3-control-v1/pack_evidence.py restore \
  --bundle research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-v3-control-v1/publication \
  --out /absolute/new/restoration
python3.12 /absolute/new/restoration/execution/extracted/prepared_f5_runtime_v3.py audit \
  --execution /absolute/new/restoration/execution \
  --expected-execution-sha256 5c0a1473b7d2637a7e97666db68262d038c628f937f3b0ea0338fd83f949f98d \
  --out /absolute/new/audit
```

Audit performs no native solver attempt. Never invoke `execute` against these
consumed records. All earlier registrations and three confirmation sets remain
closed. The full goal stays active: a separate SAT runtime/budget gate, historical
exposure union, strong matched reference, comparison-host build, calibration,
resources/order and fresh paired qualification remain outstanding. This result
is classified **accounting/correctness**, with comparative ratios unknown.
