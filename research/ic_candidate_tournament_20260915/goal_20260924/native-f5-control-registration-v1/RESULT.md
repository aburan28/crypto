# Native F5 disclosed-input control completes and independently replays

The sole invocation of registration
`5a50e01c02fa7e14f227dea413cb9e569f73d68d684a68541e49fd357660f060`
is **consumed and closed**. The original preexecution-frozen native checker
reports `PASS_NATIVE_SOURCE_BOUND_F5_CONTROL_AUDIT`: source-bound execution
admitted, target complete, and scalar `24886` independently verified for the
already disclosed synthetic n17 point `[52411,72106]`.
Never retry, resume, extend or dispatch an archive restoration.

## Acceptance and the actual execution

[PR #1294](https://github.com/aburan28/crypto/pull/1294) accepted the complete
source/custody implementation and actual preregistration on `main` at
`2026-10-03T21:54:57Z`, merge commit
`367920ad035168fefa0965ddd9a957e7e4e664bb`, accepted head
`fd64db667d5be5c0e05aef9cdc55c2c462e415c0`. All 17 applicable checks completed
successfully before merge; the two advisory checks were neutral and there were
no review threads. [The exact original gate response](result-v1/data/source-acceptance.json)
is retained. The prior implementation
[PR #1284](https://github.com/aburan28/crypto/pull/1284) was already accepted.
Dispatch waited for the complete gate.

The immutable compiled source remains
`3eb99b4d2d934059834174a4d75789bf410b916a`; its executed source, actual worker
example and compiled resources are byte-identical in the accepted combined
head. Child/worker SHA-256 is
`c4a2f160faeb57a2451171a2f105c05c77d57abad76e13eb8112f57ac47c20b0`.
The original frozen dispatcher/checker is
`5055542e85dc805c6a0c8ee9e09548d53f8bf2a84892e15e02695d580221e8b1`.
The full before-execution capsule remains in the accepted
[registration publication](publication); no source was resealed after dispatch.

The [original claim](result-v1/data/consumed.json) was created before launch.
One native worker consumed the exact frozen 7,223-byte mathematical job with
the cleared, explicit one-thread/CPU/no-cache environment. No stdin delivery
error, watchdog timeout, output cap or source gate failure occurred.
[Raw stdout, empty stderr, receipt, terminal and PID ledger](result-v1/data/execution)
retain the complete invocation and confirmed worker drain. No new ordinary
query ran, and the auditor launched no native child.

| Trial, zero based | Query coefficients `(a,b)` | Public query point | Original result | Witness |
| --- | --- | --- | --- | --- |
| 0 | `(57707,37266)` | `[20087,50635]` | Proved no geometric three-sum | None |
| 1 | `(43768,33418)` | `[39701,75542]` | Proved no geometric three-sum | None |
| 2 | `(45669,63073)` | `[111691,70749]` | Verified decomposition and recovered scalar | Geometric indices `[26,36,46]` |

The original independent audit SHA-256 is
`994ab5bf1f3153eb6d22425eb15647a152cc1ebd14b68347147c279479a9ed69`.
It checks full frozen source/input custody, the consumed claim, exact argv and
environment, delivered input/output bytes, build identity and complete PID
ledger. It independently reconstructs the seeded query law, both finite-group
negatives including torsion/repeated/identity-pair cases, witness group sum,
cofactor/orbit transport and `[24886]G = Q`. The inner mathematics-only receipt
does not itself admit runtime; the outer frozen source/transport audit does.

## Observed timing and its limits

These values are uncalibrated physical macOS ARM64 diagnostics on the frozen
Apple M4 Pro host. Busy-lock serialization does not establish L2 timing
isolation, an A/A noise bound or calibrated operations. Memory peak is unknown.
This disclosed point and query seed do not establish fresh performance.

| Exclusive target-dependent producer phase | Original nanoseconds |
| --- | ---: |
| Query generation | 36,127 |
| PDP, including both negative attempts and success | 22,420,441,207 |
| Target relation verification | 417 |
| Target descent | 9,583 |
| Producer scalar recovery check | 7,208 |
| Sum / complete producer diagnostic | **22,420,494,542** |

The recorded outer online stopwatch is `22,420,494,042 ns`, 500 ns below the
exclusive phase sum; the charged producer diagnostic is the exact phase sum,
not a claim that these are the same observed stopwatch. The outer worker
interval is `22,451,179,708 ns` and includes startup/preparation validation.
The original independent external audit took `10,480,790,791 ns` and is outside
the producer interval. These records therefore cannot establish the primary
independently verified online speedup. Matched rho/incumbent costs, normalized
operations and speedup remain unknown. Do not rank this diagnostic against the
different SAT or historical F5 controls.

Preparation is historical: 216 ordinary attempts, 61 witnessed relations,
155 independently proved geometric negatives, 32 dependencies, no duplicates,
rank 29 and all 29 column logs verified. Geometry has 63 points; there are 62
usable points before folding and 29 columns. Original Python preparation
provenance remains explicit. There are **zero new ordinary queries**, so this
control supplies no new natural-yield estimate.

## Publication and portable replay

[The result manifest](result-v1/manifest.json) and [seal](result-v1/seal.json)
bind every copied raw file. Publication uses create-only writes and retains
the original live claim, terminal and frozen audit by their recorded digests.
It neither reruns the worker nor replaces the original admission checker.
Portable replay checks the exact old registration/inputs/claim/terminal,
transport policy and data, and independently reconstructs all three attempts
and preparation again. Its context is explicitly postexecution; it executes
no archived code or Mac binaries. The original raw operator JSON is retained
in the terminal/audit files; no separate captured operator stream is claimed.

```sh
target/debug/icprog f5-control-replay-result \
  --publication research/ic_candidate_tournament_20260915/goal_20260924/native-f5-control-registration-v1/result-v1 \
  --out /tmp/native-f5-consumed-result-replay-NEW.json
```

Mutation controls must reject a replaced original audit/claim, changed output,
thread policy, extra PID, altered acceptance response and fresh promotion even
when a copied publication inventory is resealed. Result replay on another
platform establishes portable data/math checks, not native Mac-binary execution
on that platform. Actual CI results remain pending until this PR's exact head
passes; local results are retained separately.

The full goal remains active. Both native disclosed-input target controls now
have original source-bound admission. New native natural-yield evidence for
the exact F5 and CryptoMiniSat pipelines, failed-attempt accounting, complete
historical exposure exclusions, strong references, IC1 identities and a new
frozen fresh paired one-target comparison remain required. The three historical
confirmation sets remain closed. There is no fresh winner or globally-fastest
claim.
