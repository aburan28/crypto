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
