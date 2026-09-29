# Standalone GPU known-scalar arithmetic diagnostic

Frozen before execution on 2026-09-29. Follow-up to crypto PR #913.
Repository baseline: aa38ca38eff9fcce5b51fcdb70910aea8343616e.

This experiment evaluates public, known scalar multiplications on synthetic
fixtures. It contains no discrete-log solver, unknown-scalar target, walk,
distinguished-point service, or campaign deployment.

## Hypotheses and controls

1. Reduced tau-NAF may reduce known-scalar GPU arithmetic cost relative to
   binary NAF using the same backend and fixtures.
2. A precomputed linear squaring map may reduce the cost of squaring relative
   to evaluating it as multiplication. Table creation/upload must be priced.

Use a 2x2 comparison: binary NAF / reduced tau-NAF, each with generic-multiply
squaring / linear-map squaring. Retain every row. No shipping GPU baseline is
claimed: this is a portable reference backend with one scalar per thread,
polynomial-basis arithmetic, and an identical Itoh-Tsujii inverse in all arms.
The linear square uses precomputed images of the polynomial basis vectors.

## Inputs and correctness

Reuse the 24 training and 24 independent holdout cases for each of m=83 and
m=131 from PR #913's immutable run-02.json. Verify its SHA-256 and outputs
against the archived Sage output digests. No challenge generator or target
is used. New GPU curve metadata carries EC1 aliases and full curve UIDs;
the training panel's first synthetic point fixes each record's generator.

First independently recompute outputs using Python integer/polynomial field
arithmetic. Compile the same arithmetic header for the CPU and compare both
squaring variants and both recodings against every output. Test zero, one,
infinity and the rational two-torsion point as well. On device, check every
output before accepting any measurement. A mismatch aborts the experiment.

The throughput workload repeats the 24 fixtures 1,024 times, deterministically
shuffled, to make 24,576 independent known-scalar evaluations per panel.
Repeating fixtures fills the GPU; it does not increase the number of unique
correctness cases. Record this explicitly. Block size is fixed at 128.

## Timing and accounting

Record compile/JIT, fixture preparation, CPU recoding, table construction,
upload, kernel events, result download, and verification separately. Also
time allocation/upload/kernel/download as a complete warm GPU invocation.
Report a cold arithmetic pipeline diagnostic including one-time costs; do
not conflate it with a rho or ECDLP pipeline. Reused recoding and buffers in
device-only measurements must be labelled amortized.

Warm each arm. Five A/A pairs for binary NAF with generic squaring, followed
by seven interleaved rounds, reversing order on alternate rounds. Calibrate
a common repetition count from the baseline to aim at >=50 ms per timed
sample, capped at 16 launches. Report actual durations and counts; if capped
below the target, flag short samples rather than claiming precise rates.

Record the exact GPU, compute capability, CUDA/CuPy/Python versions, CPU and
memory/NUMA visibility, launch geometry, compiler options, register counts,
clocks/power/temperature when visible, hashes and before/after load. Pin host
CPU and attempt local NUMA placement; report failure or unavailable metadata.

Screening target: candidate/control median paired kernel cost <=0.90 on both
degrees and both panels, with improvement exceeding the largest A/A deviation.
This is a stage screen, not a confidence-interval-qualified performance claim.
Full rho runtime, walk iterations/s, and end-to-end DLP speedup stay null.

## Bounds and stops

Use one ephemeral Modal RTX-PRO-6000 function, max_containers=1, no retries,
600-second function timeout and 540-second benchmark process timeout. No
persistent deployment, volume, endpoint, or worker fleet. GPU unavailable or
authentication/network failure means blocked, not a fabricated measurement.
Do not switch to another account or route around an authorization rejection.
Keep failed receipts. No automatic rerun to obtain a favorable measurement.

RunPod may execute the same one-shot container on an existing authorized GPU
pod, but is a portability option, not a second experiment required for closure.
Container exit alone does not necessarily stop RunPod billing.

At the published Modal rate checked 2026-09-29, RTX PRO 6000 GPU time is
$0.000842/second: 600 seconds would be $0.5052 for the GPU component alone.
CPU, RAM, image build/startup, and other billable time are separate; this is
not a hard total-dollar cap. Source: https://modal.com/pricing .

## Decision

Success requires verified device execution, not merely a CUDA source file or
CPU emulation. A blocked cloud run may close the packaging/host-validation
deliverable, but leaves GPU results explicitly pending. Preserve archived
CPU receipts unchanged.
