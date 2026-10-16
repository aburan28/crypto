# Native baseline at the requested prime-field size

Status: **partial**. The native repository collector was executed at
`p = 2305843009213693951 = 2^61 - 1`, with three inverse pairs and twelve
deterministic seeds. Each run exhausted 100,000 target trials with zero
relations and no recovered scalar. These are capped stage diagnostics.
The screenshot's 1,000/2,000-candidate ranking and supplied degree-7 map
were not located in the inspected default-branch code and were not reproduced.

| Workload | Runs | Target trials per run | Relations | Verified scalar |
| --- | ---: | ---: | ---: | --- |
| 20-bit positive control, 16 inverse pairs | 1 | 62,676 | 33 | yes |
| Requested prime field, 3 inverse pairs | 12 | 100,000 each | 0 each | no, trial cap |

The prime-field instance is `y² = x³ - 3x + 111`, with prime group and
subgroup order `2305843010818082053`, cofactor one, and generator
`(1, 1796364306280342617)`. Its exact raw ICV1 identity and point encodings
are in the manifest. Registration of this new model remains pending; it
must not be treated as a registered leaderboard candidate.
The three stored representatives are emitted in the manifest together with
the synthetic known-answer target. The collector never receives the target
scalar. This new instance is not a transport of the screenshot's instance.

The source library is unchanged from `de0c9391480155ebdd5a8495d2c79e4508d4795d`.
The driver calls `prime_fast::find_a3_curve`, `FastFactorBase::build`,
`collect_relations` and `solve_relations`. It checks subgroup membership,
replays every relation with ordinary point operations, and independently
checks any recovered scalar with `[k]G = Q`. The control replayed 33 relations
and recovered `596767` after 62,676 target trials. No relations existed to
replay in the large-field runs. A same-host, same-binary deterministic replay
is retained separately; it is not independent validation.

## Frozen protocol and reproduction

The [protocol](protocol.json) was committed before execution. Commands:

```sh
cargo build --release --no-default-features --example prime61_factor_base
timeout 300 target/release/examples/prime61_factor_base 20 100000 1 16
timeout 300 target/release/examples/prime61_factor_base 61 100000 12 3
```

`sh research/prime61_factor_base_20261016/reproduce.sh NEW_DIRECTORY` builds
the driver, repeats both workloads, replays the requested-field stream and
compares every JSONL byte with the saved records. That wrapper passed in
this session. It refuses an existing output directory.

The positional arguments are field bit length, per-seed target cap, seed count,
and inverse-pair count. At 61 bits the driver asserts the exact requested
modulus. The deterministic seed interval is `2026101600..2026101611`.
The known-answer target is one fixed point per curve; twelve seeds are
independent collector streams on that point, not twelve different targets.

Raw records, stderr, source/binary SHA-256s, source revision, and the ecbench
host capsule are in [results](results/). The host is Linux x86_64 on AMD EPYC
9V74, with five visible logical CPUs and a 16-GiB container memory limit.
No wall timing or performance isolation is claimed. No operation-cost or
IC/rho ranking is possible: target trials are not total operations, no rho
reference was run, and none of the 61-bit attempts solved. No security,
asymptotic, or nonexistence conclusion follows from zero relations.

## Requirements and remaining work

| Requirement | Status | Evidence or gap |
| --- | --- | --- |
| Use code from `aburan28/crypto` | verified for baseline | native `prime_fast` backend and source hash |
| Execute at either requested field size | verified for capped baseline | exact prime asserted; all twelve records retained |
| Positive control and deterministic replay | verified | control scalar; byte comparison receipt |
| Reproduce screenshot's 1,000/2,000-candidate search | blocked | exact driver, candidate policy and degree-7 map not located |
| Recover and verify the large-field scalar | incomplete | all twelve runs hit their frozen trial cap |
| Rank candidates by full operation cost | not attempted | no candidate search, full-cost counters or completed large-field solves |
| Register new model and regenerate its browser entry | pending | raw exact parameters retained; no leaderboard promotion |

The source lookup is the open [problem](../../problems/screenshot-factor-base-driver-unlocated.md).
Locating its path, branch, PR, or original files is the next discriminating step.
Once identified, retain its exact candidate/base/map semantics, validate a
large-field subgroup and map, and preregister the requested candidate pools
and paired solver runs. This baseline does not discharge those requirements.

## Visuals and applicability

Repository validation: `cargo test --release --lib -j 2` passed with
4,154 passed, zero failed and 117 ignored. The initial sparse-checkout build
failed on an absent tracked identity-certificate fixture; it was materialized
and the full suite then passed. Both logs are retained in `results/`.
The example test target compiled successfully (zero unit tests); behavioral
validation is the verified CLI control and complete byte replay above.
`rustfmt --check`, `git diff --check`, shell syntax, all recorded SHA-256s and
the native backlog-card validator passed. Conductor is not installed, and
the referenced `.conductor/WORKFLOW.md` is not tracked at the source pin.

[Run report](reports/run-20261016T0720Z/report.typ),
[PDF](reports/run-20261016T0720Z/report.pdf), and
[execution diagram](reports/run-20261016T0720Z/execution.svg) describe the records.
A diagram is used because the control and requested-field runs have different
curves and base sizes and do not provide comparable cost measurements.

The canonical scoreboard, progress timeline, leaderboard and lab browser were
checked for applicability. No S, full-pipeline cost, IC/rho ratio, accepted
candidate, ecbench session or registered model is added here, so those
end-to-end graphs are unchanged. No old finding is superseded. This is a
new bounded baseline observation, not a measured improvement.
