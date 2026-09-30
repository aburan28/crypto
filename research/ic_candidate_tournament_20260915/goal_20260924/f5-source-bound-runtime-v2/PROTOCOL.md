# Source-bound F4/F5 development pipeline, version two

Status: implementation and preregistration, **not dispatched**. This is a new
bounded development control, independent of every terminal historical pilot,
paired arm and F5 correctness registration. Those attempts remain closed.
No complete v2 family admission, measured cost or speedup exists yet.

Hypothesis: the accepted native bounded MatrixF5 pipeline can collect a full
relation matrix, solve its logs and recover one supplied n17 development target
within a new frozen envelope while the complete Python controller, interpreter,
native source/build/binary and submitted JSON job are bound before execution.
The adapter also accepts the existing bounded MatrixF4 path, but this protocol
registers **only one F5 arm**. These engines use bounded Macaulay elimination;
F5 adds a row criterion. They are not full incremental Gröbner-basis packages.

## Exact instance, input and envelope

The unchanged native worker is from accepted source commit
`765c3c5f19032bd852163805f257c56babef2040`, complete source manifest
`c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7`,
worker SHA-256
`94caf3d67e57dde09488763ec19e791aca54bbb35f1968c73e42d37a210b436a`.
Its retained source/build/dependency bytes are independently checked before
registration and again in execution/audit. Staging assets never runs the worker
or generates a new query. The complete newly accepted controller snapshot will
be explicitly recorded in the preexecution receipt; an evolving checkout is
not executable evidence.

`panel.json` fixes exposed point `[52411,72106]` on
`EC1N17Ckb1hbbe2b5b6b1e6`: field degree 17, modulus `z^17+z^3+1`, kb1,
prime subgroup order 65587, cofactor two, generator `[43693,23339]` and admitted
Frobenius eigenvalue 17184. No scalar is submitted. Standard-subspace dimension
six gives 63 geometric points, 62 usable distinct nonidentity points before
folding and 29 sign/Frobenius columns. No isogeny is used; unknown conductors
remain null. This kb1 toy control does not establish ECC2K-130 fidelity.

The new algorithm seed is `2026093032`, frozen before native queries or any
exact oracle on this stream. The job fixes three summands, degree-three
MatrixF5, dense final relation LA, 8192 reduction-node budget per PDP, 100000
conflict limit in the native configuration, batch size eight and 512 maximum
ordinary trials. The unchanged target policy also allows up to 512 sampled
`aG+bQ` attempts. Actual dispatch and every attempted batch must reconstruct
from the source-bound report. Incomplete/unsupported/failed branches remain
charged; exact feasibility is assessed only after the native producer ends.

The host envelope is physical macOS ARM64, one worker, one target, unbounded
memory and a 7200-second controller watchdog. The native worker has a
7170-second cap, reserving 30 seconds for preflight and ending source gates.
Every child inherits the controller's owned process group; all listed native
thread pools are fixed at one. This adapter admits no native plugin that forks
or detaches, and no cross-platform execution of the retained macOS binary.
Portable CPU subprocess controls and transport replay do not validate a Linux
F4/F5 build. Per-child RSS is diagnostic, not total process-tree peak memory.

## Reference, accounting and acceptance

The correctness reference independently reconstructs field/group arithmetic,
the usable factor base, Frobenius quotient coefficients, scalar-prime matrix
and rank, all column logs, actual target descent relation and `[k]G=Q`.
Postexecution exact three-sum enumeration retains feasible misses and ordinary
status mix; it supplies no witnesses to the measured worker. Natural yield is
reported with rate uncertainty and the rank-stop selection boundary disclosed.
The native embedded build identity, effective configuration and observed
collector/solver paths must match the declared method.

Success requires complete pre/post Python/interpreter/native/source gates,
full 29-column rank from verified ordinary relations, verified column logs,
one relation-derived target scalar, exact online and exclusive cold ledger
closure and identical independent source/math replay after archive transport.
A scalar alone or the passing guided production controls do not meet this gate.
Candidate, workload and run IDs remain unresolved until registration freezes
the full executed surface. The new canonical candidate includes the controller
as well as the native method; native-only admission is subsidiary diagnostic
evidence and must not silently substitute its old candidate ID.

One-target online time starts at the native first target-dependent computation
after reusable preparation and ends after its independent general-group scalar
replay. All unsuccessful target attempts and native checking remain charged.
The external Python audit is recorded separately and remains outside this
native interval. Startup and reusable preparation stay outside online time;
cold accounting remains supplementary. Failed targets have null verified
online cost. Native timeout preserves bytes and source gates where available,
but incomplete native output has unknown rank, solved-target count and costs.

This is a development correctness/source control on an exposed point. Native
wall readings are uncalibrated diagnostics: no incumbent or rho arm is
registered, and no competitive speedup, normalized S/floor ratio, fresh-target
qualification, family promotion, global optimum or larger-field extrapolation
is admissible. Passing full source/math admission is a prerequisite for the
new exposure/reference/calibration adapter and fresh paired tournament.

## One-shot dispatch, retention and closeout

Merge the reviewed implementation and this protocol after applicable exact-head
checks pass. In a new result PR, freeze reusable native assets and all accepted
Python/package/interpreter sources before deriving the canonical records. Record
the accepted controller commit and full external execution SHA-256 in a
committed preexecution receipt **before the only native invocation**. Verify
static encoder feasibility against the exact accepted source and panel; that
does not generate new query witnesses. Dispatch only the externally hashed
read-only registration. No native build/identity command is run outside the
source-bound controller as an unregistered feasibility probe.

Stop at the first completed native solve, the frozen attempt caps, or either
registered watchdog. No resume, second dispatch, restart, post-outcome seed
change, budget extension or retry-until-success is allowed. A negative result
keeps the family gate open and requires a distinct future protocol, not a rerun.
Publish the registration, all sources/assets, sealed stdin job, raw stdout and
stderr, native process records/gates, independent admission and natural-yield
audit. Use `publish_f5_runtime_v2.py` for terminal successful-controller records,
including native timeouts; use `publish_sat_runtime_failure_v3.py` for failed or
partial controllers without inventing mathematical admission. Commit a durable
content-hashed archive and transport receipt, then update status and scoreboard.

Typical commands from the accepted repository, using new output directories:

```sh
python3.12 research/ic_candidate_tournament_20260915/f5_runtime_inputs_v2.py --out NEW_ASSETS
python3.12 research/ic_candidate_tournament_20260915/f5_runtime_registration_v2.py \
  --repository ACCEPTED_REPOSITORY --assets NEW_ASSETS \
  --panel research/ic_candidate_tournament_20260915/goal_20260924/f5-source-bound-runtime-v2/panel.json \
  --out NEW_REGISTRATION
# Commit the preexecution receipt before continuing.
python3.12 research/ic_candidate_tournament_20260915/f5_runtime_pipeline_v2.py \
  --registration NEW_REGISTRATION --expected-execution-sha256 FROZEN_SHA256 --out NEW_EXECUTION
python3.12 research/ic_candidate_tournament_20260915/audit_f5_runtime_v2.py \
  --execution NEW_EXECUTION --expected-spec NEW_REGISTRATION/execution.json --out NEW_AUDIT.json
python3.12 research/ic_candidate_tournament_20260915/publish_f5_runtime_v2.py publish \
  --registration NEW_REGISTRATION --execution NEW_EXECUTION --audit NEW_AUDIT.json \
  --expected-execution-sha256 FROZEN_SHA256 --out NEW_BUNDLE
python3.12 research/ic_candidate_tournament_20260915/publish_f5_runtime_v2.py replay \
  --bundle NEW_BUNDLE --expected-execution-sha256 FROZEN_SHA256 --out NEW_TRANSPORT
```

This control consumes none of the three closed confirmation rounds. The active
goal remains complete F4/F5 and SAT family admission followed by a new frozen,
fresh paired incumbent/rho comparison; a merged plan is not its completion.

## Native interface and publication prerequisites

The consumed v1 invocation in PR #1058 failed before execute_job because it
submitted `[null]` to a `Vec<u64>` input. Version two submits an empty seed list;
output fixture provenance remains `[null]`. The v1 registration and source are
unchanged and closed. Its raw evidence and unknown mathematical costs remain
separate. The fresh algorithm seed is fixed prospectively above; no oracle has
examined its ordinary-query stream.

Before a measured v2 invocation, complete the separately frozen
[actual native Job/Config deserialization control](interface-control/PROTOCOL.md),
retain and independently replay its source/build/input/output evidence, and
merge its passing controls with this adapter. That control never calls
execute_job or produces an IC query. Freeze and commit the external measured
invocation receipt against that accepted snapshot before the single IC run.

Version two publishes raw finite resource seconds as roundtrip decimal strings
and hashes the original unchanged metrics file. Native integer clocks and
scientific phase costs stay integer nanoseconds. These resource diagnostics
are uncalibrated and cannot substantiate a speedup. Float-free publication and
transport controls are required before measurement; no resource number is
rounded, silently discarded or converted into a scientific operation count.
