# P-192 / P-224 high-prime construction continuation

Frozen before execution, 2026-10-09. User request: keep trying the new native
isogeny algorithms on P-192 and P-224, extending explicit constructions beyond
degree 1009. The prior 40-degree search and its failures remain unchanged.
Input/source baseline: 7f1b92b1b, with the published source orders and seed 1.

First record a 20-second P-192 degree-509 diagnostic on the original CLI and
sample only that owned child. This diagnostic does not discharge the requested
construction search. Preserve its raw output, process receipt and sampled stack.

Then attempt P-192 degree 509 and P-224 degree 521 with 600 seconds each, and
P-192 degree 1021 and P-224 degree 1031 with 1800 seconds each. The larger
degrees minimize the maximum Frobenius eigenvalue order among split primes in
[1010,1099] for each source (170 and 206 respectively). Their maps still require
construction and independent replay; this selection is a recorded heuristic.

Use one construction process at a time, the existing 8 GiB sampled RSS cap,
25 ms sampling, and at most 100 ms of transient unavailable monitoring during
exit. Persistent monitoring failure, timeout, map-count disagreement or failed
verification is retained as an unresolved outcome. Keep every raw file and
receipt, the algorithm source/executable digests, stage traces and environment.
Elapsed values are operational receipts on a shared host, without a speed claim.

Accept a split-degree construction only with two maps, successful CLI checks
and independent kernel, codomain, exact rational-map and fresh public subgroup
transport replay through the previous study's verifier. Register new models
before citing them in final visuals. Keep candidate and constructed-map counts
separate. Preserve the old reports; this continuation gets its own results,
editable figure, checked PDF and manifest. Existing publication checks remain
mandatory, including the root release-library gate.

Native setup may replace redundant divisor scans with a field-valued divisor
sieve, with exact modular-coefficient tests and the complete standalone suite.
Stage logs are enabled only for this diagnostic/continuation. Do not weaken the
mathematical checks or treat a stage improvement as an end-to-end measurement.
