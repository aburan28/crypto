# Draft for a fresh n17 external-SAT one-target transport

Status: target design draft, not a registered target experiment or dispatch
authority. The original MatrixF5 panel reached rank 29/29; the first
100,000-conflict CryptoMiniSat panel reached 22/29, and the separately frozen
[one-million-conflict panel](../native-sat-million-registration-v1/result-v1/RESULT.md)
reached 29/29 with all logs independently replayed. Both SAT registrations
are consumed and closed. The three historical confirmation sets remain closed. The disclosed
[transport controls](controls-v1/RESULT.md) passed exporter parity on all
three points; accepted-CMS stdin parity passed two of the three frozen
60-second controls and a separately labeled longer diagnostic for the third.
The first frozen control failure remains. A distinct post-buffer-marker CMS
binary was built, and its [nine-role disclosed transport result](result-v1/RESULT.md)
passed an independent data-only audit. The target-free SAT worker source is
staged, but no one-use SAT target capsule, fresh public-point card or target
solve has been admitted.

## Question and bounded scope

Can the accepted external CryptoMiniSat engine be part of a complete,
source-bound, single-public-point IC solver whose online interval begins only
after every target-independent process launch, executable load, curve/base/log
preparation, and native role setup? The target is the exact educational n17
Koblitz representation `EC1N17Ckb1hbbe2b5b6b1e6` used by the native
`ecbench` candidate/workload records, not a challenge target. The registry's
older `EC1N17Ce1hdfbf24105ef5` label uses a different encoding record for
the same field and curve; it is not the new paired workload identity.
The target input is one public point with no known scalar supplied to the
worker. Its factor logs must come from an independently audited, complete,
natural CryptoMiniSat ordinary panel. A separate F5 preparation cannot be
silently substituted for the SAT candidate.

## Why the existing adapter is not sufficient

`PreparedN17Target::solve_external_sat` already has the seeded `[a]G+[b]Q`
law, bounded attempts, full-point model lifting, independent scalar replay,
and the five-phase online clock. `SatTargetBackend` has start/query/complete
callbacks; the durable target journal supports the SAT row shape. The current
ordinary native backend launches `exporter` and `cms` per query. The accepted
exporter CLI receives the point as command-line arguments, and the ordinary
backend gives the accepted CMS CLI a CNF path. Reusing that backend as a primary one-target timing
would include process startup in the target PDP phase and fail the declared
online interval. It also lacks a separate SAT target capsule, original frozen
worker, native role audit, and one-use target claim.

## Preinitialization feasibility gate

The accepted CMS source archive has SHA-256
`467b1c3d00a7d6e893332b4d8b42c6326301974d22885aa745f1f926da050323`.
Its `src/main.cpp` enters `readInStandardInput` when no positional CNF filename
is supplied, prints a flushed reader line at verbosity 1, then parses stdin.
That is a source observation, not yet a binary-compatibility result. On
disclosed synthetic n17 points, test the exact pinned binary against this
path before adopting it for the target. The printed line precedes stdin
stream setup and `DimacsParser` construction, so the compatibility probe's reader line is
**not** proof that every target-independent initialization step has finished.
For a primary online interval, add a source-bound readiness marker after
reader/parser initialization and before the first read, prove byte- and
model-equivalence to the accepted binary on disclosed controls, or supply an
equally auditable post-initialization handshake. A new marked binary has its
own source and binary hashes in the target candidate; never call it the
unchanged accepted executable. Source review found that the earlier
[one-line sketch](cms-stdin-ready.patch) emits READY before
`StreamBuffer` allocates its 148,576-byte buffer; that constructor also calls
`fread`. The stricter [post-buffer patch](cms-stdin-ready-postbuffer.patch)
constructs the buffer without reading, emits READY only after that allocation,
then primes the stream and parses the target CNF. Its file-input default keeps
the original eager-read path. The distinct marked binary passed the
[fixed three-way disclosed control](result-v1/RESULT.md); the exact binary and
exporter still need a new one-use target capsule and its own audit.
Start **one stdin-piped CMS child for
each of the `max_queries` permitted attempts** in distinct process groups
with only target-independent flags. Wait for every child to complete the
post-initialization handshake before opening the online interval. A child exits after one
CNF and cannot serve a later attempt. Record every binary hash, argv,
environment, PID, readiness interval and resource use. The frozen resource
envelope includes the memory and descriptors of the whole idle pool. No child
may receive the public point, derived query, ANF or CNF before the online
clock opens. Source review and transport receipts must prove this boundary;
stopwatch subtraction cannot.

After the clock opens, deliver each source-verified CNF through its held stdin
pipe with a bounded nonblocking write, close it, collect the model/status,
enforce the per-attempt deadline and drain the group. An unused child receives
no CNF and is cancelled and drained with a retained receipt when the target
run ends. Compare source/model outcomes with the accepted cold-file CLI on
the same disclosed points, including SAT, geometric negative,
budget-inconclusive and timeout controls. A stdin-path failure is a retained
feasibility result; do not route around it in a measured run. If the pinned
binary cannot consume a prestarted stdin stream reliably, design and freeze
a separate in-process CryptoMiniSat library server with equivalent source and
model checks. Its startup must also precede the online interval.
The three disclosed points in the closed
[native SAT control](../native-sat-control-registration-v1/RESULT.md) give two
proved geometric negatives (`[40991,73355]`, `[73003,104622]`) and one
source-verified witness (`[59775,2910]`) for exact transport parity. Use new
low-budget and watchdog controls for inconclusive and timeout paths; do not
run or mutate the consumed control registration itself.
The [Rust prestarted-stdin probe](../../../../examples/ic_sat_stdin_probe.rs) is a
disclosed control for one retained CNF/ANF pair. It pins the CMS binary hash,
runs a cold-file and prestarted-stdin path with the same conflict budget, checks
the manifest's public point and all three export digests, retains both outputs,
checks each SAT assignment against the retained ANF and CNF and a full-point
lift, independently enumerates the three-sum class, and reports
status parity. Its [disclosed input plan](CONTROL_PLAN.md) fixes the three
audited points and exact source-file hashes. The
[retained result](controls-v1/RESULT.md) records two frozen passes, one frozen
failure and its separate post-hoc diagnostic. A passing probe is transport evidence,
not a source-bound complete target solve or a performance measurement.

The accepted exporter source constructs the standard factor-base predicate
and geometric base before it *uses* the explicit target, but its existing CLI
still reads target arguments at launch. The new source-bound
[one-request exporter mode](../../../../examples/koblitz_pdp_export.rs)
accepts **no target arguments or target environment values** and constructs
the field, predicate, geometric base and static basis encoding before a
flushed READY marker. Prestart one exporter child per permitted target attempt
as part of the frozen resource envelope. Each child then accepts exactly one
bounded point and blind instance ID from stdin after the online clock starts,
exports the source files and exits. Its release binary passed all three
disclosed export-parity controls; it is not the accepted ordinary-exporter
binary. Preserve the accepted ANF, CNF-XOR,
manifest and model semantics. Cross-check ANF, CNF-XOR and Magma export bytes
exactly; compare the deterministic manifest fields and source-instance
identity exactly while excluding timing and the new `transport` metadata,
which must differ between a cold CLI and a prepared process. The prepared
manifest explicitly warns that its internal whole-process timer includes
idle time before target delivery; it cannot be used as an online measurement.
Check independent full-point lifting against
the accepted exporter for a frozen disclosed corpus. Do not import known target
logs or derive `[a]G+[b]Q` before the clock. Record all exporter process
startup, request/response, source bytes, hashes, watchdog and drain receipts.
The `PreparedChild` transport in
`src/bin/prepared_sat_worker/native.rs` implements pinned-binary startup,
one exact flushed readiness marker, a durable no-input READY receipt,
deferred bounded stdin delivery and unused-child cancellation. The
[`ic_exporter_prestart_probe`](../../../../examples/ic_exporter_prestart_probe.rs)
passed its three disclosed-input controls, including malformed/oversize
rejection and unused-child drain. These controls may not be treated as a
target-worker admission without a new source-frozen target registration and
independent audit.

## Target implementation and admission

Use a new SAT target scope/worker and a new scientific registration; do not
weaken the F5 target scope or a consumed old controller. Freeze the complete
source, vendored dependencies, accepted/rebuilt native assets, target-free
config, preparation binding, resource envelope, startup policy, algorithm
seed and max attempts. Publish the source and registration with external seals
before generating a public point. Each of the four strict source descriptors
must name its role (`cms`, `f5`, `incumbent` or `rho`), the same canonical
`kb1` curve ID, sealed registration, source manifest and full archive digest,
and declare a target-free, unexecuted build. Unknown fields are rejected so a
descriptor cannot carry a target or scalar. These declarations require an
independent publication audit; the point-card generator does not prove them.
Then publish a scalar-free point card using `ecbench`'s
`hash_to_subgroup_v1` public-target law. It carries the native workload ID,
full workload digest, seed, index, derivation counter and exact point; the SAT
worker replays that workload before consuming the card. This is the same
one-point workload representation accepted by the incumbent IC and rho arms.
The `ecbench_workload_id` is a point/fixture cross-check here; a later campaign
must mint its own canonical workload ID from the frozen resource and online
accounting policy. The card alone does not qualify an IC or rho run.
The existing native F5 worker can accept a point in a newly frozen post-card
target config, provided its complete source/build archive is published before
the card and the final registration proves byte-for-byte source and executable
identity with that publication. Its old one-use registration remains closed;
the old auditor's `fresh_paired_qualification: false` still needs a separate
chronology and same-point campaign audit. Generate the incumbent and rho
`ecbench` workloads from the card's seed and check their workload digest and
point before either measured solve.
The card's source-publication hashes and creation order must be audited
independently. Bind that exact card to the
one-use claim before dispatch. The original frozen checker must verify the
published archive and card, the preparation's original source-bound audit and
recovered factor logs, native child
identities/inputs/outputs, every attempt, source ANF/CNF/model, geometric
negatives, final scalar and the five exclusive online phases. A bad model,
no full-point lift, budget exhaustion, timeout, transport error, or failed
durability check remains a distinct retained row. No failed target is a
verified recovery.

Preparation and native role startup are cold/setup costs outside the online
interval and reported separately. The interval opens immediately before the
first target-dependent query computation and closes after the independent
`[k]G=Q` replay. It includes target query generation, request transfer,
ANF/CNF construction, CMS solve, source/model/relation checks, descent and
recovery. The five phases must sum exactly to online wall time. Include every
failed attempt. An unentered phase is explicitly zero with a reason; a missing
phase is null and blocks an online total. No inference from the ordinary panel
or from a mock backend can establish this runtime admission.

| Event | Accounting boundary |
| --- | --- |
| Exporter and CMS pool launch, READY handshakes, reusable base and factor-log setup | Cold/setup; finish before online opens, with original receipts |
| Seeded `[a]G+[b]Q` construction and target validation | `target_query` |
| Exporter request/response, equation and CNF construction, bounded stdin transfer, CMS search, including failed attempts | `target_PDP` |
| Independent source/model verification and full-point lift | `target_relation_check` |
| Cofactor/log recovery from a verified witness | `target_descent` |
| Independent `[k]G=Q` scalar replay | `target_recovery_check`; stop online at its completion |
| Cancel and drain unused, never-fed CMS children; frozen audit | Post-online custody costs, separately retained; these children must have done no target work |

The first native run is a disclosed correctness control, not fresh comparison.
Only after that original frozen audit passes, freeze an exposure census,
canonical EC1/IC1/workload/run manifests, one new public point and the same
resource envelope for SAT, F5, incumbent IC and strong rho. Run the four arms
on that point in a predeclared order with failures retained; each has an
independent scalar replay. Use native `ecbench` for cross-method measurement,
and obtain a host-level isolation/noise receipt before promoting a CPU wall
speedup. No stage winner is a global winner.

The current `ecbench` `ic.pipeline` method runs its own in-process oracle and
does not attest an external one-use SAT capsule. Add a source-bound external
method/import boundary or an equivalent audited native integration before
putting SAT and rho into one `ecbench` claim. The imported target row must be
the exact original frozen audit, with its five online phases, child receipts,
public point, candidate/workload IDs and resource envelope; never relabel it
as an `ic.pipeline` execution. `ecbench` verification must reject a changed
or unaudited source row and keep the one-use registration closed, so replay is
data-only. If that adapter is absent, report the paired runs separately with
speedup unknown under the `ecbench` claim rules.

Keep SAT-native and F5-internal work unpriced in a counted `S` until one
predeclared, common operation boundary measures them; a group-operation lower
bound is not a complete cost. Report reusable preparation and process-start
costs separately from the measured online interval, and leave any unsupported
cold total null rather than summing intervals from incompatible run conditions.
