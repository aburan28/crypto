# Paired device test of squaring-loop unrolling

Frozen before device execution; follows PR #962. This is a standalone
known-scalar arithmetic diagnostic on the immutable PR #913 synthetic inputs.
No discrete-log solver or challenge workload is included.

## Fixed comparison

Select `--study square-unroll`. Run four arms: binary NAF and reduced tau-NAF,
each with linear squaring using default compiler unrolling or the opt-in
`LINEAR_SQUARE_NOUNROLL=1` pragma. The comparator for each candidate is the
same recoding with default unrolling. Cache compiled modules separately by
field degree, squaring implementation, and unrolling option. Record compiler
options and kernel attributes for every arm.

Use the same 24 training plus 24 holdout fixtures at each of m=83 and m=131,
and the same field/curve identities and verified Sage output digests. First
verify all 384 outputs across the four arms. Any mismatch stops the study.
Smoke mode stops after correctness and produces no timing summary.

## Timing and decision

Keep 24,576 evaluations per launch and 128 threads/block. Calibrate one common
repetition count using the first control, aiming for 50 ms, capped at 16.
Collect five A/A pairs separately for each recoding's default-unroll control.
Then run seven rounds, reversing the four-arm order on alternate rounds.
Checkpoint each verified invocation, each A/A pair, and each complete round.

Report candidate/control median paired kernel cost separately for both
recodings. A panel passes the screen only when that ratio is <=0.90, the
improvement exceeds its own control's maximum A/A deviation, and all control
and candidate samples reach 50 ms. An overall recoding screen requires both
degrees and both training/holdout panels to finish and pass. Preserve every
sample and unsuccessful outcome. This is a screening result, not a qualified
runtime speedup claim; no confidence interval or production-walk conclusion
is inferred. Register counts and spill bytes remain compiler diagnostics.

The original four-arm recoding/squaring study remains available as
`--study original`. Its historical receipts are immutable.

## Isolation and bounded execution

Follow current AGENTS.md section 10: timed runs go through the repository's
`tools/isolated_bench.py run`, with a full visible SMT sibling group reserved
and at least one other CPU left free. Keep its default quiet-machine and
contention thresholds. Preserve the isolation receipt; a contended run,
unavailable receipt, or unresolved user thread on reserved CPUs blocks an
accepted timing result. Do not silently relax isolation or retry.

Smoke mode requires no CPU timing isolation and records no performance result.
The launcher allows one process group for at most 540 seconds; on timeout,
terminate the entire group and preserve the last valid benchmark checkpoint.
Provider limits remain one GPU, no retries, 600-second Modal function timeout.
Container exit does not itself stop a RunPod pod's billing.

No new device timings are available when this protocol is committed.
Host-stub tests validate argument wiring, module separation, pairing, and
failure handling only. GPU correctness, image build, and device measurements
remain pending working authorized compute access.
