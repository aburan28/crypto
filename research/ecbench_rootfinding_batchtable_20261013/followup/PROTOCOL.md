# Target-pinned relation collection follow-up: frozen protocol

The question is whether stopping at a verified target pin, retuning the prime
factor-base size, or changing the root-finding oracle lowers cold one-target
work on the registered 20-bit prime subgroup. The reference is `rho.negation`.
Each public target is checked by exact `[d]G = Q`. All setup, failed probes,
relation collection, linear algebra and verification remain in the ecbench
record. Counted GAE and native work are reported separately. macOS L0 timings
are diagnostics, not admitted wall-time comparisons.

The frozen sweep uses prime search bits 20, seed 59297, eight fresh public
targets with target seed 20261016, one warmup and three measured rounds,
alternate arm order, algorithm seed 20261016, and a 120-second per-run limit.
Every IC arm uses max_trials 5000000 and the same guarded 16-jump walk.
The arms are matched rho, its A/A repeat, B64 batched pair table with
checked full rank, B40/B48/B64/B80 batched pair table with target pinning,
and B64 `rootfind-subtract:m=2` with target pinning. No workloads, arms,
units or seeds change after results are inspected. Preserve errors and
timeouts. Replay every measured run and compare paired cold work and phase
counts. Select the minimum *measured cold counted work* among the four base
sizes, then freeze an independent seed for a confirmation panel versus B64
and matched rho. The sweep selects; only the confirmation tests the choice.

The target-pinned change is accepted only if all measured arms verify and
the replay audit is identical. Report relation hit rate, trials to pin,
rank at pin, table occupancy, memory and phase costs. A size that wins only
in an unpriced or noisy metric is not a validated cold-work winner. The
subtraction arm answers a current-binary comparison, including all failures.
No scaling claim follows from this one curve size.

For the native walk experiment, first establish exact point and coefficient
stream equivalence with the generic walk on fixed seeds and restarts. Profile
walk and oracle separately. A one-lane `Fp64` implementation is valid only
if it preserves guard, restart and relation behavior; count field operations
and include conversions and setup. Multiple independent lanes with batched
inversion form a separate arm. They must retain the same generated target
multiset and exact answers, and account for lane fill/drain and memory. A
descriptive L0 timing is not a wall-time advantage. Keep unsuccessful code or
run attempts as separate artifacts; do not overwrite a sealed run.
