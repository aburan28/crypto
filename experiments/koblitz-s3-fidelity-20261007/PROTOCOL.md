# Paired S3 inversion: isolated two-curve measurement protocol

## Question and frozen decision rule

Does sharing one field inversion between the two S3 target queries shorten the
**verified one-target online interval**, and is that shortening different on
the two supported Koblitz instances below? Compare each target only with its
own baseline on the same curve. The primary interaction is the difference
between the curves' mean paired log speedups. This is an implementation
comparison; an IC-versus-rho claim additionally needs a matched strong-rho
solve of each public point under the primary one-target contract.

| Curve (registered ICV1 slug) | Model | Subgroup order | Role |
| --- | --- | ---: | --- |
| `icv1-f2m41-tm2308219-7f48b14a` | `y²+xy=x³+1` over the recorded GF(2⁴¹) polynomial basis | 549756390943 | First curve |
| `icv1-f2m53-tm56619371-dac20a85` | `y²+xy=x³+1` over the recorded GF(2⁵³) polynomial basis | 21044858204113 | Second curve |

Both use `KoblitzCurve::new(0,n)`, 244 signed Frobenius orbit columns, the
same low-x base construction, the same target-span stop rule, relation seed
20260928, and 14 relation/index workers plus the main thread. The only
within-curve algorithmic change is the paired S3 inversion in PR #1518.
The curves are **different mathematical instances**: they have separate curve,
factor-base, candidate and workload identities. No comparison of their raw
solve times alone answers the question above. The implementation currently
does not admit `KoblitzCurve::new(1,41)` or `(1,53)`; that is an implementation
constraint, not a claim that those curves do not exist.

The source snapshots in
[`../koblitz-s3-pair-query-20261007-v4/`](../koblitz-s3-pair-query-20261007-v4/)
are the baseline and candidate. Historical N53 Mac candidate IDs and timings
remain frozen there. For this Linux panel, `candidate_id: null` for all four
curve/variant combinations until the Linux binaries, exact base digests,
source hashes, compiler, flags and every manifest field are frozen. Compute
new `IC1...` hashes from sorted, compact JSON with no run fields; use each
target's canonical one-target workload hash and run IDs
`<candidate-id>W<workload-id>R<repeat>`. The N41 feasibility run found 20,008
usable base points and 244 folded columns; N53's frozen record has 25,864 and
244. Re-enumerate and independently verify both counts on the measurement
build before assigning IDs. A nominal column count is never an `fb` count.

### Hypotheses and acceptance

For each independently generated target, let `A` be the baseline and `B` the
candidate. Run `ABBA` or `BAAB` in fresh processes. The target's log speedup
is `d = log(geomean(A_online_ms) / geomean(B_online_ms))`; its two repeats
estimate host noise but count as **one** independent target. The preregistered
primary test is whether `mean(d_N41) - mean(d_N53)` differs from zero, using
a two-sided 95% Welch interval over independent targets. Report both
curve-specific geometric mean ratios and their intervals regardless of the
interaction result. A separate positive speedup claim for each curve needs
its own multiplicity-adjusted interval, with Holm correction across the two
curves. Do not infer an interaction because one curve's interval excludes 1
and the other's does not.

Every counted target must have four completed online intervals, identical
public points per within-curve pair, verified scalars, matching semantic
witnesses, complete exclusive phase costs, and an independent point replay.
Every failure, timeout, OOM, and rejected fidelity sample stays in the raw
ledger. A target with an algorithmic failure cannot become a win by selecting
its successful repetitions. Report completion rate and keep the aggregate
verified speedup unknown if the planned panel is not fully verified. Apply
noise exclusions using the host counters and the fixed criteria below, before
inspecting which variant ran faster.

## Host and fidelity gate

Use a **physical Linux host** with an auditable, exclusive cgroup v2 isolated
CPU partition covering at least 15 execution CPUs on one NUMA node, separate
housekeeping CPUs, whole SMT sibling coverage, memory binding, IRQ routing,
fixed performance frequency, no CPU quota, and no foreign runnable tasks on
the reserved CPUs. Record physical CPU model, microcode, core/SMT/NUMA
topology, cpuset partition readback, actual affinity of every worker thread,
frequency/thermal readings, `/proc/interrupts`, CPU steal, runqueue delay,
CPU/memory/I/O PSI deltas, cgroup throttle/OOM counters, memory peak and
`memory.numa_stat`, kernel command line, binary/source/workload hashes,
the exact execution order, and raw failures. Use the repository's ICMS
L0–L3 checks as the per-run diagnostic grade and the stronger host-level
receipt required by the cryptanalysis isolated benchmark rules as the claim
gate. The physical-host receipt is mandatory even when ICMS says L2 or L3.

The `isolab` strict policy in this repository already implements the
host-level cpuset/NUMA/IRQ/PSI/perf receipt. Its strict tier A and passing
noise checks are the intended execution route. First run the read-only host
preflight. A RunPod Pod may provide useful correctness and exploratory
profiling, but its container-visible affinity is insufficient evidence of
host-wide exclusivity. If its strict preflight fails, retain the rejected
receipt and label any subsequent Pod timing **exploratory**; do not turn an
ICMS grade or a custom numeric score into a controlled-speedup claim.

`perf stat` on the isolated execution CPUs records task-clock, cycles,
instructions, context switches, migrations, faults, branches and cache
events when supported. Report `IPC = instructions / cycles` with the event
counting scope and `time_running / time_enabled`; unsupported or multiplexed
events remain missing. IPC is a mechanism diagnostic, not a fidelity score or
a replacement for online wall time. Profile each variant in a separately
labeled pass: `perf` over a whole process also counts reusable base and index
setup, so those counts must not be described as target-online IPC. The main
ABBA/BAAB wall panel uses the same binaries without the optional profiling
pass. Measure the profiling overhead with an A/A control.

The benchmark executable owns its timer: target-dependent computation starts
after reusable base/index setup and ends after recovered-scalar point replay.
The five exclusive online costs are `target_query`, `target_pdp_charged`,
`target_relation_check`, `target_descent`, and `target_recovery_check`; their
sum must equal `target_online_after_reusable_setup` within timer rounding.
Record base/index setup and whole-process wall separately. A known-scalar
fixture is built outside the timed interval, and the solver receives only its
public point. No cross-target table is shared.

## Sampling and stopping before results are seen

1. Generate **12 pilot targets per curve** from one frozen SHA-256 seed law,
   modulo each curve's subgroup order; use separate previously unseen public
   points and one workload ID per point. Include a separate A/A baseline
   control on 12 points per curve. Balance `ABBA` and `BAAB`, as well as which
   curve runs first in each block. Use one warmup target outside the pilot.
   The inputs under [`pilot/`](pilot/) have already been generated; their
   scalar witnesses are retained only for independent post-run checking.
2. Preserve every pilot run. Require the A/A 95th percentile of
   `|log(A1/A2)|` to be below `log(1.05)` for a 10% interaction target, and
   require no fidelity failures among the pilot samples used for variance.
   If either gate fails, report the pilot and repair or change the host before
   starting a confirmatory panel. Never loosen the gate after seeing ratios.
3. Run `cargo run --release --bin s3_pair_power -- A.json B.json 1.10` on
   the **clean pilot's one-row-per-target pair summaries**. It calculates
   `ceil((z_0.975+z_0.8)^2 (s41²+s53²) / log(1.10)²)` independent targets
   per curve for 80% power at two-sided 5% alpha. Freeze the confirmatory
   size at `max(24, ceil(1.20 × calculated_n))` per curve, using new seeds;
   the pilot is excluded from the confirmation test. The 20% inflation is a
   prespecified allowance for pilot variance uncertainty. If that size is
   infeasible on the available host, report the cost and leave the
   confirmatory result pending rather than silently shrinking it.
4. Run the full fixed-size panel without effect-based early stopping. Compute
   the Welch interval on target-level `d` values, and a target-cluster
   bootstrap as a sensitivity check. Preserve ABBA/BAAB order and all raw
   host and solver records. If valid-sample replacement is necessary because
   of host noise, retain both attempts and generate replacement targets from
   the next prespecified seed; algorithmic failures are not replaceable.

The 12 N53 targets in the earlier **unisolated** holdout have sample SD
`0.209277` for paired log ratios. They imply about **38 independent targets
on one curve** to detect a 10% within-curve effect and **76 per curve** to
detect a 10% difference between curves, under equal variance. With the
preregistered 20% inflation, the provisional confirmatory size is **92 per
curve** (184 distinct targets, 736 ABBA/BAAB solver invocations), plus the
pilot and A/A controls. These are planning calculations, not measurements of
N41 variance, and the strict-host pilot supersedes them before confirmatory
target generation. At 5% interaction the same noisy SD gives about 289
targets per curve before inflation; at 15%, about 36.

## Execution and evidence

Build both source snapshots against one checked-out `crypto` dependency
revision, one `Cargo.lock`, Rust version, target architecture and release
profile. Record executable SHA-256 before and after the panel. Never reuse
the archived macOS binaries on Linux. Keep the baseline source snapshot
unchanged. Reserve the host before the first warmup and release it only after
the block ends; do no builds or audits during timed runs.

The manifest registers the immutable sources as the Cargo examples
`s3_pair_baseline_frozen` and `s3_pair_candidate_frozen`; build them together
with `s3_icms_adapter`. The adapter turns a single solver invocation into
ICMS's exact `online_one_target` window and refuses a mismatched target,
worker count, missing verification, or incomplete phase sum. Generate one
ICMS `fixed_file` spec per point and variant, with the file SHA-256, 15
declared threads/CPUs, and the Linux binary hashes. Its `command` argv is:

```text
s3_icms_adapter SOLVER 41|53 0 244 20260928 PUBLIC_POINT.json 14 METRICS.json
```

For ICMS, `METRICS.json` is the `{metrics_path}` command-adapter substitution.
Use the same public-point file in both variant specs. For a target session,
set `warmup: 0` and `repetitions: 1`, then list the two specs as `A B B A`
or `B A A B`; warmups use a separate public point. This retains the full
four-execution order in one session. ICMS's own spec/workload
IDs are additional measurement labels; retain the canonical EC1/IC1/W/R
identities required by the candidate catalog in the final rows. The ICMS
grade is diagnostic unless the separate host-level strict receipt passes.

For each run, retain candidate/workload/run IDs; curve and point encodings;
binary/source hashes; resource envelope; all five phase times; target and
rank S3 calls; attempts including misses/timeouts; verified relations, novel
rows and final rank; actual base points and folded columns; peak RSS and
NUMA placement; scalar and replay certificate; start/stop timestamps; host
fidelity receipt and `perf` status. Recompute semantic digests across arms,
then independently verify each target and its relation witnesses with the
checked Sage launcher `/Volumes/SSD990/cryptanalysis/sage`. Save
`/Volumes/SSD990/cryptanalysis/sage --runtime-info` before Sage validation;
keep that startup check outside online timing. Preserve the raw run even if
Sage rejects it.

No matched rho result is in this protocol's existing evidence. Should the
result be compared with rho, run strong single-target rho on every same
public point under the same resource envelope and report
`rho_online_ms / IC_online_ms` separately. Batch throughput, setup-inclusive
time and the S3 stage cost do not substitute for that primary comparison.

## Status at protocol freeze

| Requirement | Status | Evidence or remaining step |
| --- | --- | --- |
| N41 S3 path executable | Diagnostic pass | One local unisolated pair recovered scalar 123212651130; baseline 181.602292 ms, candidate 82.272583 ms; key relation fields matched. Independent Sage replay and host isolation pending. |
| N53 correctness across fresh targets | Verified historical control | PR #1518's twelve one-target pairs and checked Sage replay; its CPU ratios are exploratory. |
| Pilot public inputs | Frozen, unmeasured | Twelve SHA-256-derived points per curve under `pilot/`; `SHA256SUMS` checks the files. Exact workload IDs await the Linux candidate/curve manifests. |
| Physical-host receipt | Pending | Saved `runpod` SSH endpoint refused connection and no RunPod API key was present in this session at protocol preparation. |
| Strict two-curve pilot and A/A | Pending | Requires a reachable qualifying Linux host and frozen Linux binaries/targets. |
| Pilot-derived fixed sample size | Provisional only | 76 per curve from historical N53 variance; replace with the clean two-curve pilot. |
| Confirmatory curve interaction | Pending | No confirmatory targets have been generated or measured. |

The protocol is versioned before the larger panel. The local N41 diagnostic
preceded this freeze and is not part of the pilot or confirmation.
