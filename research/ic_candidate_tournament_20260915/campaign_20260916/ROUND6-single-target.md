# Round 0006 pre-registration: single-target rho parity

Written before any measured stage of this round was run.

## Objective

Reach parity with the shipped signed-Frobenius rho on complete cold
**single-target** public-hash ECDLP jobs, on the five small Koblitz cells the
operation measures (n13a0, n17a1, n19a0, n19a1, n23a0). This is the panel the
[continued operation](RESULTS.md) left open at 3.04 times rho's instructions
and 1.52 times its native time (round-0004 selection `folded_lift_batch4`).
The 16-target panel is a separate workload and is not touched.

Parity is the rule already in the contract: candidate/rho upper paired 95%
limits and every curve-cell ratio at most 1.10, in user-space instructions and
in native process wall, on confirmation and on replay. Promotion over the
incumbent keeps the standing gates: at least 20% lower instructions and native
wall, upper paired 95% limits below one, no cell more than 10% worse, on
confirmation and replay, then a fresh-process replay with the independent
Python checker.

## Parent, seed, budget

Incumbent: round-0004's promoted `folded_lift_batch4` — the frozen
`runs/round-0004/source_candidates/folded_lift/source` with
`{batch_trials: 4, linear_algebra: sparse, max_trials: 4096, solver: pair_table,
summands: 3}`. Fresh campaign seed 2026091606; target count 1; pilot profile;
1,800 paired-job budget; one pinned CPU; 8 GiB address-space cap; 60-second
watchdog; Valgrind 3.22.0 `Ir`; `--require-native-progress`. Rho is the
existing per-target signed-Frobenius solver API on the incumbent's executable,
exactly as in every earlier round.

## Where the single-target cost was

Per-phase means from round-0004's confirmation receipts (n19a0, millions of
Ir): rho 15.7 in total, of which 8.7 curve construction and target hashing,
2.1 general-arithmetic final verification and 0.4 startup are shared by both
arms, 4.5 is the rho solve and 0.08 reporting. The IC winner spent 17.4 on
factor-base and table setup, 6.6 on descent, 3.4 on collection, 3.5 on linear
algebra, 1.6 on log certification and 1.6 on reporting. Function-level reads of
the retained Callgrind profiles attribute most of that to general `BigUint`
point arithmetic while sampling the base, Fermat field inversions, thread-pool
start-up for `rayon` iterators with one thread, hash maps keyed by big
integers, and certificate serialisation through `BigUint` formatting and a
fragment-by-fragment JSON formatter. The relation search itself was small.

For parity the IC-specific work must fit in roughly the rho solve plus ten
percent of the job; on n13a0 that is about 2 million instructions for
everything the IC arm does beyond the shared phases.

## Candidates (six arms, 1,584 paired jobs)

1. `incumbent` — as above.
2. `combined_descent` — the 16-target panel winner's source (round-0005),
   unchanged, on this single-target panel: its descent changes were only
   measured with 16 targets per job.
3. `fast_report` — round-0005 source plus a worker-only change: single-word
   coordinates format as `u64` decimal (same text), and the report is
   serialised once into a buffer and written once. Isolates the serialisation
   cost. The algorithm, certificates and phase boundaries are untouched.
4. `tiny` — round-0005 source plus a new library module,
   `koblitz_tiny_ic.rs`, and a worker branch that runs it for
   `pair_table`, three-summand, prime-degree jobs with the `SubgroupOrbits`
   recipe (every job in this operation). The worker's phase dumps, fixture
   construction, rho branch and general-arithmetic final verification are
   unchanged; the report carries the same fields and types, written directly.
   The module's declared mechanisms:
   - the identical factor-base **point set** (same seeded draws, same lift
     test, same cofactor projection and orbit closure), checked by a test
     against `FactorBaseSpec::materialize` on every cell;
   - single-word curve arithmetic with an extended-Euclidean field inverse
     instead of the Fermat inverse;
   - a folded pair-sum table keyed by the least rotation of normal-basis
     coordinates (constant on Frobenius orbits), one stored witness per sum
     orbit, with the rotation count locating the pair;
   - **density-sized table rows**: rows `C_o + P_j` are built for the least
     number of orbits `t` such that a random rest has at least 2.5 expected
     pair-sum witnesses, `λ·(1 − ((K−t)/K)²) ≥ 2.5` with
     `λ = |F|(|F|+1)/(2r)`, or every row. This is a fixed function of the
     public base size and subgroup order; on the sparse cells (n17 and up) it
     builds every row, as before;
   - walked probes `R_{t+1} = R_t + S` with `a_{t+1} = a_t + s` in place of
     one scalar multiplication per probe; decomposition scans the base in
     small blocks and stops at the first witness; each witness is re-added in
     the group before it becomes a relation;
   - incremental row-echelon solve over `Z/rZ` with every column logarithm
     certified by `[ℓ]G = [h]C` inside the solve and again before reporting;
   - descent by walked `[a]G + [b]Q` probes, the answer verified as
     `[d]G = Q` in single-word arithmetic and then, as before, in the general
     arithmetic by the worker.
   Base support, `m = 3`, no direct relations, and every verification
   obligation are unchanged. Configuration keys `linear_algebra` and `sparse`
   have no effect on this path (the system has at most a few columns); the
   registry states the configuration it was run under.
5. `tiny_batch1` — the same source with `batch_trials: 1`.
6. `tiny_fulltable` — the same source with the table rule disabled (every row
   is always built): the control that attributes the small-cell gain to the
   rule rather than to the rest of the rewrite.

Development evidence before freezing (not a claim): on the round-0004
confirmation fixtures, natively, all 60 `tiny` outputs pass the frozen checker
with factor-base fingerprints equal to the incumbent's in every cell; on one
case per cell under Callgrind the total instruction ratio to rho was 0.80 to
0.88.

## Boundary, floor, class

The unit and boundary are the contract's: complete user-space `Ir` from
startup to termination. The weak floor is still `K` instructions for the
`K`-column full-rank collector. Nothing here changes the mathematical
coverage bound `binomial(|F|+2, 3)` of a fixed signed base with three
summands. Class: **engineering** — the algorithm's arithmetic and its
implementation moved; the ratio to the floor did not, and no non-generic
advance, arithmetic-complexity, family-wide or cryptographic-size claim
follows. The rho reference is the shipped implementation whose per-job cost on
these cells is mostly fixed setup; a rho specialised the same way has not been
measured and would be cheaper. Both statements go into the report.

## Honesty conditions

- Every measured stage runs on fresh fixtures from the new seed; no prior
  holdout is reused.
- A failure anywhere is retained; nothing is re-run to obtain a better draw.
- The report states the target count (1), the rho API, the parity rule, the
  worker change, and that this is an instruction/native-time implementation
  result on five small Koblitz cells.
- The single-target and 16-target panels stay separate: no ratio combines
  them.
