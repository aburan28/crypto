# Compact paired whole-solve timing amendment

The V4 order-67 batch independently replayed 198 of 198 pairing scalars.
At distinct-target prefix 64, its phase-summed pairing CPU was
`1,318,707 / 976,037 / 1,105,577 ns` over three rounds, while the
coarse `all_in_pairing_cpu_ns` was
`4,754,960 / 3,890,210 / 4,808,711 ns`. The latter interval subtracts
timed rho calls from process time but includes repeated JSON serialization
and reporting from both arms; it is an upper bound, not a fair solver-only
pairing timer. The phase sum excludes some input and between-phase work.
Neither number alone decides a cold crossover.

V5 keeps all public target lists, H selections, finite-field BSGS,
pairing formulas, same-target rho seeds, restart cap, independent replay and
three rounds fixed. Read and validate the public input once, charge that
cost identically to both methods, then time each complete method from cold
field/group setup through its final scalar verification. Do no per-target
JSON work during either timer; emit one record after both complete. Take
cumulative snapshots at prefixes 1, 4, 16 and each cell endpoint, including
64 and 66 for order 67. Run pairing before rho in rounds 0 and 2, and rho
before pairing in round 1. Preserve phase and primitive counts from V3/V4
as the accounting ledger, and use V5 whole-solve CPU only for a clean L0
paired comparison. Add common public-input validation CPU to both totals.

Stop if either method fails, if any scalar arrays differ, or if independent
replay fails. A V5 CPU ratio remains shared-host L0 evidence with three
rounds; it does not establish a counted-operation or scaling crossover.
