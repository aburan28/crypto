# Further degree-83 rho tuning

The fresh held-out panel measured reference Fold256/DP6/escape-off against the
previously tuned Fold64/DP6/escape-off control. Both use production u64 scalar
arithmetic in the same executable. Geometric mean control/selected whole-process
cost ratio: **1.069651413981749**, target-block 95% interval
**[1.0431783781923283,1.0892379529034781]**. Arithmetic mean process costs were
**2.987007664533333 s** and **2.780545233400000 s** (ratio 1.0742524986298727).
All five target means improved; each mean contains three paired fresh seeds.
The pooled 15-run medians reverse that ordering slightly: control 2.833963098 s,
selected 2.850071739 s. The 15 paired ratios have median 1.0757736750707487 and
minimum 0.9943479875332358; an individual launch can regress. These descriptive
summaries mix target/seed cases and do not replace the frozen five-target
estimate. They and all raw cases are retained rather than selecting one summary.

This is a measured engineering runtime comparison on one CPU/host and five
fresh synthetic targets. The point estimate exceeds the largest observed A/A
ratio, 1.0465408131500493; the lower interval edge overlaps that noise envelope.
Five target clusters and a shared VM limit generalization. The interval is
conditional on frozen training selection and target-orbit exclusions, not a
statement that 256 lanes is globally optimal or that every seed improves.

## What was actually tried

The native campaign froze 139 launches before timing: ten identical A/A controls,
18 arithmetic-kernel runs, 81 policy runs, and 30 held-out runs. Every launch
completed and verified. Four campaign reservations observed zero contention.
There were no process timeouts or global solve-cap failures. Internal per-lane
capped-walk restarts do occur: their counters remain in the producer outcomes,
and their work is included in the charged process/online costs. Completed and
uncensored refers to the overall solve, not to absence of these lane restarts.
The actual child wall/CPU totals were 423.47495148600007 s / 422.98710100000005 s;
they exclude builds, parent orchestration and independent external verification.

Exact Barrett reduction and conditional modular addition passed seven worker
tests, independent native arithmetic checks, checked Lean mathematical lemmas,
and nine matched complete solves with equal recorded answers, counters and
charges. They did **not** reduce measured training cost: 2.923128010666667 s
versus reference 2.9030551186666664 s. The frozen selection retained reference
arithmetic. The written proof and negative timing remain part of the evidence.
Recorded-counter agreement is not a complete intermediate walk trace.

The reference-kernel policy sweep tested lanes 64/128/256 and DP6/8/10 with
fruitless escape fixed off. All nine policies completed all nine training cases.
Fold256/DP6 was selected by whole-process mean, 2.680834025111111 s, against
Fold64/DP6 at 2.8917663447777775 s. Selection was immutable before holdout.
The full losing-policy table, every target and seed, and costs are retained.

## Identity, measurement and correctness

The exact legacy curve is `icv1-f2m83-t6151469093347-9eba6cd5`, a=1,b=1,
field z^83+z^7+z^4+z^2+1, r=8569786107849059, h=1128547018. It differs
from the primary a=0 ECC2K confidence-gate fixture. The field backend and
signed Frobenius quotient (166 phases) are unchanged. The generic reference
scale sqrt(r/166) is unchanged; calibrated operation-equivalent S remains null.

Native C++ regenerated eight fresh point-only fixtures and verified exclusion
from the original point, all K600 orbits, and eight previous optimization points.
The validation sidecar never entered child argv or environment. A separate
native checker replayed all **139 answers**, checked the full stage/fixture
counts and all nine arithmetic semantic pairs, with **zero invalid rows**.
Its receipt SHA-256 is
539751e5dc58c0390b2d36df52e016d674eac51c61d770ffb5ffb6e2597e9988.
The separate native file audit passed 2,439 recorded hash bindings across 455
distinct bound files and 1,036 snapshots, including all 139 execution, producer
and observation records, all clean environments/resources and four frozen
matrices. Its receipt SHA-256 is
8511411695b7de2a2fd693b64122c42a75dc89bf5b9b694ea8c60af710eb5bd1.
External verification took 0.886574336 s wall, 0.833018 s user and 0.051324 s
system with peak 109188 KiB RSS, on reserved CPU3 without observed contention.
This combined fixture-and-answer check is outside both solver process intervals.

Whole-process costs include startup, parsing, curve/normal-basis/jump setup,
walk, collision recovery, built-in scalar verification and output. Online costs
are secondary: geometric ratio 1.0752029765269837, 95% interval
[1.0503821796055077,1.0942886571132544]. Each solver uses one OS thread on CPU3,
12 GiB address-space limit, a 60-second cap, and the same clean environment.
The 5 ms exit-poll interval includes scheduling overhead and is not a strict
upper bound. Physical-host exclusivity and independent-host replication are not
established. Native operation counters are retained; primitive costs are not
calibrated, so this does not promote an operation-based boundary ledger row.

No new IC run or IC/rho ratio was made. The previous matched cross-method study
and all previous PDF pages remain historical evidence. Ratios from different
target panels are not multiplied and new rho times are not paired with old IC
times. The scoreboard, IC leaderboard and progress timeline were checked and
left unchanged because this rho-only engineering comparison adds no admitted
cross-method ratio, calibrated S or new curve identity.

## Evidence and reproduction

- [Frozen protocol](PROTOCOL.md), [machine protocol](provenance/protocol.json),
  [baseline and scope](BASELINE.md), [selected config](selected_config.json).
- [Native statistics](records/panel/statistics.json),
  [native descriptive medians/minima](analysis/descriptive_summary.json),
  [A/A noise](records/panel/aa/noise_summary.json),
  [arithmetic selection](records/panel/arithmetic/selection.json),
  [all policy scores](records/panel/policy/selection.json),
  [all held-out rows](records/panel/holdout/records.jsonl).
- [Worker validation](checks/worker_validation.md),
  [full library log](checks/full_lib_prebench.log),
  [runner tests](checks/native_runner/test.log),
  [native answer receipt](verification/final_receipt.json).
- [Native file-provenance receipt](verification/file_provenance_receipt.json),
  with its explicit compiled-source coverage and admission limits.
- [Fixture generation and verifier](verification/README.md),
  [arithmetic proof](audit/ARITHMETIC_PROOF.md),
  [Lean proof map and exact scope](audit/PROOF_MAP.md),
  [protocol audit](audit/PROTOCOL_AUDIT.md).
- [Native controller](native_runner/README.md), raw process directories under
  records/panel, and [isolation records](records/isolation).

Build and run the example on one point-only input with the selected config:

```sh
cargo build --release --locked --example n83_rho_speed_eval
target/release/examples/n83_rho_speed_eval point.jsonl WALK_SEED research/n83_speed_20261012/selected_config.json NEW_RECEIPT.json
```

Use the native runner and isolation commands in native_runner/README.md for a
new complete campaign; output paths are create-new. Original absolute workspace
paths remain in frozen receipts. Reproduction relocates inputs through a new
protocol rather than rewriting old records. Source and binary pins are in
provenance/frozen_sources.sha256 and provenance/frozen_binaries.sha256.

The implementation/protocol was committed as 9da9ff2c and preflight pins as
adddc355 before the first timed child. Base main was c62555b7 (merged PR1722).
Conductor was attempted and unavailable; ownership and scope are in BASELINE.md.
The missing /usr/bin/time wrapper failed before any verification child started;
its error is retained in checks/replay_attempt1.stderr.log. Direct native
monitoring then completed successfully. Earlier compile/proof attempts remain.

## Requirement status

| Requirement | Status | Evidence / remaining boundary |
| --- | --- | --- |
| Faster exact degree-83 rho on fresh workloads | Verified on this host/panel | 1.069651 ratio, five-target interval and all seed rows |
| Fair comparison against already tuned rho | Verified | Same executable/backend, frozen current64 versus selected256 |
| Preserve negative approaches and all cases | Verified | Barrett and all nine policy rows retained; no failures discarded |
| Exactness and independent answer replay | Verified | Seven worker tests, native arithmetic and 139 scalar replays |
| Mathematical proof and formalization | Checked mathematical lemmas | Four Lean declarations; machine-word bridges remain written/native |
| Full local validation | Verified | 3992 library passed,112 ignored; seven worker and22 controller passed |
| Comprehensive PDF and visual record | Delivery artifact | New additive report, earlier60pages preserved, benchmark copies |
| General optimality, calibrated S, L2, independent host | Not established | Outside the admitted scope of this runtime comparison |

This study and its new reports are versioned additions. Previous evidence is
unchanged. MANIFEST.sha256 binds the final artifacts; the file-provenance receipt
binds the original recorded sources, binaries, commands, inputs and outputs.
The additional native descriptive module passes five tests. The pinned Lean
build, axiom check and fresh kernel replay pass locally; the new CI workflow
reproduces those checks on future changes. These are local validation results;
remote CI and landing status belong to the follow-up PR.

Delivery resumed on Oct14 and integrated main275beb93f at2ea7fd320. The final
integrated release library passes4089 tests with117 ignored; the worker passes
all seven tests. [Integration receipts](checks/integrated_20261014/README.md)
and [delivery provenance](DELIVERY_PROVENANCE.md) separate this validation from
the unchanged Oct12 measurements. Byte-for-byte measured source copies are in
[the archive](archives/measured_sources/README.md). The [70-page comprehensive
PDF](reports/run-20261012T1543Z/unified_curve_basis_report_faster_rho.pdf) and its
historical60-page input are now durable repository artifacts too.
