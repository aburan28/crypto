# Native prime walk: frozen comparison

The guarded 16-jump relation walk currently uses an affine prime group
addition at every query. The candidate uses `prime_fast` Montgomery
arithmetic with a per-solve cache of the current point and the 16 jumps,
plus lazy restart offsets. The factor base, batched pair table, oracle,
elimination and exact recovery check are unchanged. The new target source
is explicit as `targets=walk-prime-fast`; the old `targets=walk` remains
the baseline.

The one-lane ecbench spec fixes the registered 20-bit prime curve, eight
fresh public targets under seed 20261018, one warmup, three measured rounds,
alternate order, `prime-abscissa:size=48`, `mitm-prime-batch:m=2`,
`incremental-gauss`, max_trials 5000000, and matched negation rho plus A/A.
The primary correctness gate is identical `(recovered, verified, trials,
relations, rank, walk steps, restarts, lookups, guard probes, total GAE)`
on every paired execution. The unit test compares every `(point, a, b)`
on a seeded stream with forced restarts and exceptional additions. Replay
all measured runs. A mismatch is a failed implementation, not a speedup.

Report native field-operation tallies and any conversion overhead. On this
macOS host, all wall readings are L0 diagnostics. A one-lane primitive
timing cannot establish an end-to-end ECDLP advantage. The separate
`add_many_batched` kernel tests independent lane streams and batches their
denominator inversions; its ignored throughput test is a bounded kernel
diagnostic. It excludes lane scheduling, relation probes, matrix work and
memory, so it is not an end-to-end arm.
