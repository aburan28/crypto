# Device-readiness follow-up, 2026-09-29

Frozen before the compiler-resource experiment. Repository reference:
`7b7d4ed2880c8619dc82e97630d79941d02dcd74`, including PR #960.

## Question and fixed comparison

Does the existing, unchanged known-scalar CUDA diagnostic assemble for both
Ada (`sm_89`) and Blackwell (`sm_120`), and does its compiler report spilling?
Compile both existing squaring variants at m=83 and m=131 with NVRTC and
ptxas 12.8.93, default optimization, no register cap. The binary and tau
recodings share the same kernel, so resource counts belong to each kernel,
not separate recodings. Preserve the eight compiler logs and code hashes.

Success means all eight kernels assemble. Record registers, stack bytes and
spill-store/load bytes reported by ptxas; leave missing statistics null.
Nonzero spill counts identify a profiling question, not a measured runtime
bottleneck. Zero spilling does not establish high occupancy or good speed.
Stop on compilation/assembly failure and preserve the failing result. Each
ptxas subprocess is limited to 60 seconds. No arithmetic source changes,
register tuning, or performance measurements are part of this experiment.

## Runner correctness work

Add a smoke mode that checks the same 96 archived fixtures against all four
arms (384 outputs), with one copy of each fixture. Before timing, the full
mode must complete this device gate. Smoke observations have no throughput
summary and cannot satisfy the parent timing protocol. Retain the parent's
24,576-evaluation batches, five A/A pairs, seven alternating rounds, and
screening threshold unchanged. Smoke warms modules explicitly; later stage
cost accounting must continue to label reuse rather than claim cold cost.

Write receipts atomically and checkpoint every verified invocation and A/A
pair. Record a device launch before verification separately from successful
verification, so a mismatch cannot be misreported as no device execution.
Keep failed/interrupted partial results readable. Test local failure and
receipt handling without implying device validation.

## Access and claim limits

Unauthenticated HEAD probes of the documented provider endpoints are
diagnostics only. Do not infer credential validity from these results. Use
the normal configured network path; do not bypass access restrictions.
Retain the earlier authenticated failure receipts unchanged. If access is
still blocked, do not launch another paid job or retry for favorable results.

All outputs are stage diagnostics on synthetic known-scalar inputs. GPU
runtime ratios, walk throughput, rho speedup, and end-to-end DLP costs stay
unmeasured. This work does not implement or deploy a discrete-log solver.
