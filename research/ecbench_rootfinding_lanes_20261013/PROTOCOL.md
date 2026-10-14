# Rootfinding collector follow-up: frozen protocol

Requested work: integrate synchronized 4/8/16-lane additions; profile probe costs on an isolated host; retune B after collector change; test relation coverage against uniform-point occupancy. The public prime-search subgroup is bits=20, seed=59297, with exact target verification. The reference is `rho.negation`; A/A rho is a noise control. All phases, unsuccessful probes, setup, matrix work, verification, resident memory, and failed runs stay in the record. Prior immutable sessions remain untouched.

## Lane experiment

`specs/lanes.json` freezes 8 new public targets, target/algorithm seed 20261020, 1 warmup and 3 measured rounds, alternate order, 120 s per run, B48, `mitm-prime-batch:m=2`, target pin, and max 5,000,000 trials. Compare 4/8/16 synchronized lanes with either separate native additions or one batch inversion per round. Same-width pairs must replay the same ordered point/coefficient stream, relation rows, rank, target scalar, and verification. Charge prepared but unused tail lanes, restarts, conversions, field operations, setup and memory. This is an engineering test at one curve size; local macOS wall is L0 diagnostic.

## Probe profile and host gate

Profile field inversion, canonical coordinate conversion, jump selection, table lookup, guard hash, and restart-pool work separately. `specs/lanes-l2.json` fixes the same targets and arms at five measured rounds. An admitted wall comparison requires its L2 isolated ecbench session with A/A noise plus exact audit. If no worker is available, retain L0 profiling as diagnostic and keep the L2 gate open. The prior one-lane native/generic wall ratio 1.3898 was L0 and is a symptom to investigate, not a qualified speed difference.

## Base retune

`specs/retune-train.json` and `specs/retune-holdout.json` freeze distinct sets of 8 public targets at seeds 20261021 and 20261022. B is 40, 44, 48, 52, 56, or 64, all with batch-8 lanes and the same collector, oracle, cap, rho, rounds and stopping rule. Select one B on train by cold total counted cost including setup; holdout checks the selection without changing it. Report memory, native field work, phase breakdown, relation yield, independent rank and exact verification. B differences in lane schedules remain target-matched but are not identical probe-stream pairs.

## Coverage test

`specs/coverage-train.json` and `specs/coverage-holdout.json` freeze B48 with abscissa strides 1, 7 and 23 modulo p. Stride 1 is the ordinary base. Both panels use eight distinct public targets, seeds 20261023 and 20261024, and the same batch-8 walk, oracle, cap and measurement schedule. Strides are fixed curve-only choices before target generation; there is no target-dependent optimization. Primary metric is unique signed pair-table keys and independent relations per total cold cost and peak memory. The null is uniform subgroup occupancy, with expected probe hit rate `2K/r` for K unique negation-folded keys. Compare observed hits with binomial intervals and equal probe budgets. A result consistent with the null establishes only that this tested construction adds no detected coverage bias at this size. Any excess must survive independent public holdout targets and cold-cost accounting.
