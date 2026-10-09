# Execution notes

The original 2026-10-08 evidence is preserved. Its large probes had low sampled
resident memory before the timeouts, motivating stage diagnosis. Source review
found that the initial E4 coefficients scan every possible divisor separately
for every series coefficient. It also uses u64 cubes, which overflow at large
indices. Any correction is checked independently before construction results
are accepted.

The ignored algorithm-validation build cache was moved from the previous
study's target directory to /private/tmp/p192-p224-continuation.f9sck4/algorithm-validation.
It is absent from the prior manifest's bound executable list. Prior raw data,
reports, source snapshots and bound executables remain at their original paths.
New build artifacts use this temporary volume. An initial Conductor connection
failure was followed by a successful scope check and T-7 reservation under the
currently configured cryptanalysis project.

The first copied launcher kept its old argument interface; that failed launch
is preserved in validation/supervisor-first-launch-failure.log. The next
20-second diagnostic completed but had no PID file for sampling. Both sampler
attempts failed before collecting a stack. The corrected PID-writing launcher
uses baseline-diagnostic-v2/ and successfully sampled only its own child.
The retained sample places all sampled main-thread stacks inside the Hecke
modular-polynomial constructor; it does not identify a later stage or establish
a speed comparison. New stage traces will distinguish the constructor phases.

The field-valued divisor sieve and its optional stage logging were committed
as 97c80540b. The full release suite passed 114 tests, including the independently
enumerated divisor coefficients, a coefficient beyond u64 cube/sum capacity,
and known normalized j-series coefficients. The new protocol and diagnostics
were committed separately as b8133a7bc after adding the study to sparse checkout.

The extended full-Hecke probes preserve 600-second timeouts for p192/509 and
p224/521 and an 1800-second timeout for p192/1021. Their last completed stage
was q-series construction; later baby/giant powers had not completed within
those budgets. The p224/1031 invocation follows serially. These are operational
stage records, without a matched performance comparison.

An exact additional eigenvalue-order screen through 65537 found a degree-5
eigenline at p224/1471 and a degree-3 eigenline at p192/10453. KERNEL_PROTOCOL.md
freezes a distinct kernel-first route; it does not replace the full-Hecke
two-map obligations. Extension tests passed on both source primes and on an
independent exhaustive small-field point count. The per-map replay imports the
same unchanged existing walker arithmetic and retains every mathematical check.

The first Cargo build attempted to include the manually linked kernel helper
without a crate dependency and failed. That log is preserved. The Cargo
workspace now builds only the independent replay; the helper and supervisor
use their recorded rustc --extern builds. The kernel helper's early argument
arity was corrected before target construction. Kernel source 2fe16c12b,
the frozen executable and the separate single-map supervisor are used for the
queued kernel runs. No new map is promoted without independent replay.

For the degree-10453 certificate, the independent verifier retains the original
walker field, kernel and curve checks. Its polynomial multiplication is adapted
to exact Karatsuba recursion, with the original schoolbook source as the test
oracle on both source fields and balanced/unbalanced lengths. This changes
arithmetic cost, without skipping any torsion, subgroup, normalization or exact
map-identity check. The copied implementation and equivalence test are frozen
as verification_poly.rs and verification_poly_checks.rs.
