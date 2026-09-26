# RTX PRO 6000 throughput measurement and workstream closeout

Date: 2026-09-26
Status: Planned; GPU measurements outstanding.

## Goal

Produce a reproducible report explaining how close the workload gets to its relevant hardware throughput ceiling, what limits it, and what uncertainty remains. Do not equate GPU utilization with peak arithmetic throughput.

## Starting point

The prior work reports that [PR #860](https://github.com/aburan28/crypto/pull/860), commit `2dd3dfd5`, was merged and all eight CI workflows passed, including CUDA compilation. This status comes from the completed conversation; it was not independently refreshed for this document.

The reported fixes cover consistent exits for detected short cycles, a default trail limit, and rejection of incompatible table data and checkpoints. GPU runtime correctness and performance remain unmeasured. The tests do not establish the absence of every possible cycle.

Operational prerequisite: table-walk v3 requires a fresh corpus. Retain the trail limit for residual cycles. Earlier performance figures must not automatically be attributed to the changed build.

## Action checklist

- [ ] **Pin the environment.** Record the exact GPU edition, power limit, sustained clocks, driver and CUDA versions, source commit, compiler configuration, and workload configuration.
- [ ] **Validate on hardware.** Complete the corrected build's GPU correctness checks and retain the results with the performance evidence.
- [ ] **Establish reproducible timing.** Separate warm-up from measurement, repeat runs, report variability, and distinguish kernel time from end-to-end time. Record thermal or power throttling and competing workloads.
- [ ] **Collect a hardware profile.** Capture instruction mix, execution-pipeline throughput, memory traffic, occupancy, and stalls. Use separate unprofiled timings where profiling overhead could distort results.
- [ ] **Document an analytical ceiling.** Base the estimate on the executed instruction mix and applicable hardware limits. State assumptions about clocks, instruction dependencies, and shared execution resources. Label undocumented limits as uncertain rather than inventing precision.
- [ ] **Document an empirical arithmetic ceiling.** Use an isolated, correctness-checked arithmetic benchmark on the same GPU. Explain differences between that benchmark and the full workload; a microbenchmark is a reference, not proof of an attainable application rate.
- [ ] **Write the measurement report.** Include raw evidence, configuration, achieved throughput, estimated fractions of the two ceilings, variability, and the supported bottleneck explanation.
- [ ] **Make the closeout decision.** Close the correctness-and-measurement workstream when the criteria below are met. Track any subsequent optimization as separate work.

## Interpretation rules

- Advertised floating-point TFLOPS and AI TOPS are not the relevant ceiling for binary-field arithmetic.
- Use the same work unit and timing scope when comparing achieved throughput with a ceiling.
- Distinguish theoretical upper bounds, empirical reference rates, and measured application performance.
- High occupancy or 100% GPU utilization alone does not establish arithmetic saturation.
- Kernel throughput and end-to-end throughput answer different questions; retain both.
- Do not claim that the corrected table walk is faster than another implementation without comparable measurements.

## Report fields

| Field | Result |
| --- | --- |
| Exact GPU edition and configuration | Pending |
| Source commit and software versions | Pending |
| GPU correctness result | Pending |
| Work unit and timing scope | Pending |
| Repeated kernel throughput and variability | Pending |
| End-to-end throughput and variability | Pending |
| Analytical ceiling and assumptions | Pending |
| Empirical arithmetic reference rate | Pending |
| Achieved fraction of each applicable ceiling | Pending |
| Dominant bottleneck and supporting counters | Pending |
| Uncertainty and measurement limitations | Pending |
| Raw evidence locations | Pending |

## Closeout criteria

- [ ] GPU correctness checks pass.
- [ ] Measurements are reproducible and their variability is reported.
- [ ] Raw results and the exact configuration are retained.
- [ ] Profiling supports the bottleneck conclusion.
- [ ] Any ceiling estimate includes its assumptions and limitations.
- [ ] Outstanding limitations and follow-up work are explicitly recorded.

An arbitrary claim of “100% peak” is not a closeout requirement. A defensible conclusion that the workload is limited by memory, dependencies, or another resource is a valid outcome.

## Reference

[NVIDIA Nsight Compute Profiling Guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/) describes instruction statistics, throughput, occupancy, and performance analysis.

## Tracking and evidence contract

This plan follows [PR #860](https://github.com/aburan28/crypto/pull/860).
The saved checklist is preserved above; its GPU measurement fields remain pending.
Merging this planning document does not close the measurement workstream.

Before measurement, commit the frozen workload, seeds, build flags, run order,
repetition count, warm-up policy, cost accounting, and stop conditions in a
follow-on experiment PR. Retain a clean baseline and an A/A noise check. Any
baseline/candidate comparison must use matched conditions and retain failures.
Keep profiling runs separate from unprofiled timing runs. Preserve raw reports,
commands, binary/source hashes, and the hardware manifest in that PR; large
artifacts need hashes and durable accessible locations.

Report kernel rates as stage diagnostics. Full-workload timing must state its
boundary and include setup, transfers, reporting, and verification where those
are in scope. Leave end-to-end DLP cost and speedup unset unless a complete,
verified, matched comparison supports them. Compare ceiling fractions only in
matching units and scopes; an empirical arithmetic reference is not a proven
upper bound on application throughput.

Use follow-on PRs for execution evidence and any resulting optimization. Update
this checklist and link those PRs when their results are accepted. Do not mark
hardware validation or workstream closeout complete merely because this plan
or a compile-only check is merged.
