# Generic backend qualification — registered before measurement

The three-attempt pair-table improvement budget is exhausted
([bounded goal result](../BOUNDED-GOAL-RESULT.md)). No selected challenger
passed the 0.8 cold gates. This protocol registers the next source-bound
complete-DLP comparison: qualify generic F4/F5/SAT decomposition and larger
sparse relation-LA paths against the optimized incumbent and strong rho on
identical public targets. It is not an improvement round, a held-out claim, a
global optimum, or an ECC2K-130 crossover. Measurement is pending; this file
freezes the comparison before any new solve.

## Hypothesis and stop conditions

**Hypothesis.** At least one of the registered generic F4, F5, `sat_xor`,
`sat_cnf` or `inherited_f4` pipelines, with dense or declared sparse final
relation LA, can complete independently verified one-target solves on every
smoke and development job of the five-cell toy panel under the resource
envelope below, with full scientific admission and exclusive phase closure.

**Success (qualification).** At least one F4/F5 arm and at least one SAT arm
must each complete *every* scheduled smoke and development job. Each verified
job independently binds source/build, usable base, query law, observed engine
dispatch, stored matrix/rank, factor logs, target descent, scalar replay and
native/profile closure. A complete scalar without the checked IC path fails.
Retain separate cold/online IC and rho leaders, cell tradeoffs and all
alternatives. No familywise promotion inference is made on this development
panel. A qualified arm may later enter a new improvement registry only after a
separately reviewed reference-binding update.

**Stop / negative qualification.** If an arm cannot complete verified solves
on the panel within the declared timeout/memory, or fails admission, it remains
a table row with null competitive costs. Do not widen the envelope, drop cells,
or substitute a cheaper stage diagnostic as a complete-DLP win. Incomplete
schedules are incomplete; they are not retries of this registration.

**Inadmissible moves.** Retuning on confirmation/replay from the three sealed
improvement rounds; reusing any exposed public point; changing the operation
unit after seeing results; counting stage diagnostics as end-to-end speedups;
or claiming ECC2K-130 transfer from this toy panel.

## Accepted references and scope

Reuse the version-two reference binding accepted in PR 876 / protocol
[improvement-v2/PROTOCOL.md](../improvement-v2/PROTOCOL.md):

| Role | Selection | Source SHA-256 |
| --- | --- | --- |
| incumbent | prepared `pairinv` | `240b8daa0478aacb1f6fc248de9fd5ed9b0773c59d3ae869ced0bab8be438378` |
| ic_online | prepared `both` | `65df2c4417ec963a2b03e511719df28df51ea266a3c4a9e825307b8d1ac8c053` |
| rho | `rho_incumbent_8` | same prepared `pairinv` source, `rho_parallel_walks=8` |
| rho_online | `rho_generic_dense_16` | generic source `72da739cbc72e351e85c0184664ea60b6d2f23b54aab0e9ff45fdd1a60714d7d`, build `54316bcec91a68a5a08bcc8e92548a46ad75e2a27d8c056fe6520d79f3410962` |

All costs remain those of fully charged instrumented implementations. Missing
totals stay null. GF(2) work inside F4/F5/SAT is distinct from relation linear
algebra modulo the prime subgroup order; both must be charged when used.

## Frozen panel

Seed **2026092901**. Stages: A/A, smoke and development only. Cells: n17a1,
n19a0, n23a0, n23a1, n31a0. The n29a1 holdout is neither generated nor
inspected. One previously supplied public point per job; three process
repetitions. Export every exposed point for later exclusion, including failed
or unexecuted cases.

IC arms (see [panel.json](panel.json)):

1. `incumbent` — prepared pairinv, pair_table, tiny_gauss
2. `prepared_both` — prepared both, pair_table, tiny_gauss
3. `generic_pair_dense` — generic-v1 pair_table + dense LA (continuity control)
4. `generic_pair_sparse` — generic-v1 pair_table + sparse LA
5. `generic_f4_dense` — generic-v1 f4 + dense LA
6. `generic_f4_sparse` — generic-v1 f4 + sparse LA
7. `generic_f5_dense` — generic-v1 f5 + dense LA
8. `generic_sat_xor_dense` — generic-v1 sat_xor + dense LA
9. `generic_sat_cnf_dense` — generic-v1 sat_cnf + dense LA
10. `generic_inherited_f4_dense` — generic-v1 inherited_f4 + dense LA

The registered `prepared_both` alias maps to the accepted `ic_online` reference
role in the executable registry. It is counted once: nine IC candidates
including `incumbent`, one additional IC reference, and two rho references.
Scheduling the same source and configuration again as a candidate would create
a duplicate method identity and artificial replication.

Every generic IC arm uses summands 3, batch_trials 1, max_trials 65536, and
factor-base recipe subgroup_orbits with seed 43 and requested points 6n. Actual
usable support, orbit folding and rank are independently checked. Prepared arms
keep their existing tiny_gauss configuration. Reference arms cannot win
challenger selection.

The follow-on runner pins the new generic worker source *before target
generation*: Git tree `src` object
`6caf5dd2de704de77fd3cae5affab8b0ec3b6911`, example worker object
`745247f3d89d46cdffb6d1f28778971caae7aa21`, `Cargo.toml` object
`4177c1daa3b7f779abcf5fdeecb3b284d93b0f19`, and checked-in CI lockfile
object `a3181ddde6d0a7460f0d9a1c6e87e3bead801980`. The controlled build
then binds all actual compiled source and dependency bytes, compiler, flags and
executable. A later main-branch source change stops this registration before
fixture generation; it cannot silently change the F4/SAT candidate.

Rho arms for continuity: `rho_incumbent_8` and `rho_generic_dense_16` only
(the accepted cold and online rho leaders). Do not expand the width screen in
this registration; a separate rho-sensitivity note may follow if a new IC
leader appears.

## Schedule and resources

With 10 IC + 2 rho = 12 arms:

| Stage | Public points | Arms | Repetitions | Native/profile pairs |
| --- | ---: | ---: | ---: | ---: |
| A/A | 5 | 2 | 3 | 30 |
| Smoke | 5 | 12 | 3 | 180 |
| Development | 15 | 12 | 3 | 540 |
| Total | — | — | — | 750 |

Cap 900 pairs. Linux amd64, Rust 1.94.1, Valgrind 3.22.0, one pinned CPU and
one Rayon thread, 8 GiB and **300 seconds** per child (raised from 180 only for
this qualification; record the change). Failures, timeouts and OOMs remain
rows. Never retry or overwrite a measured execution.

Primary metric: verified one-target native online wall time after reusable
preparation through scalar replay. Supplementary: complete cold native time and
user-space guest instructions with exclusive phase closure. Report
`S = Ir / sqrt(r)`, actual B/K, paired rho ratios and the weak collector floor
where applicable. Use the same supplied public point and frozen resource
envelope for each IC and rho arm. An online gain does not imply a cold-start
or total-operation gain; a stage gain does not establish a verified online
speedup.

Audit every ordinary sampled collection query, including unsuccessful PDP
calls, before reporting natural relation yield. Keep `unresolved` and budget
exhaustion distinct from independently proved UNSAT. Report verified witness
queries, duplicate and accepted rows, novel rank trajectory, final rank,
matrix and final relation-LA costs, target attempts, scalar replay and actual
base size before orbit folding. Report base construction and memory when
measured; absent memory is unknown. Planted correctness vectors are not a
natural-yield sample. Retain the raw reports of incomplete bounded jobs and
independently check any attempted-query history used in rates. A timeout with
no complete report is censored, not a zero-yield query. Process repetitions
on one point are not additional independent points; descriptive rate
uncertainty resamples distinct public points within their fixed cells.
Because a bootstrap can collapse to zero in a zero-yield cell, also report a
conservative 95% Hoeffding interval for the mean of distinct-point yield
rates under the frozen point law. Neither interval turns a censored run into
a measured failure rate.

## Exclusions

Before generating any point, restore and verify the three sealed improvement
archives and extend the historical census with every exposed point from:

- `ic-improvement-round1-20260925` (SHA-256 `f2b0f9b5c226ddc9331449d9e37ace1c2d7a94faa5ef74396d9a0a1231b13593`)
- `ic-improvement-round2-20260928` (SHA-256 `020e0cef5f044e01613625da93d9cedaa48e875a59523f2859477c484f23b6f0`)
- `ic-improvement-round3-20260928` (SHA-256 `a2e88ae540ddfa3430b03873f5cd85df0d5b1ad12b8e59b7dea82a5b33df0a4f`)
- readiness / generic-reference-qualification / adapter-control fixtures
- every other corpus already pinned by the round-three runner

The panel file byte SHA-256 is
`83c640a03b4239b918851f6f1b8450e2fe27dc710fb99f306e3481a99d8875cf`.
An altered seed cannot make an old public point fresh. Keep failed and
interrupted preparation evidence in the exclusions.

## Delivery

Implement a versioned runner that prepares the panel under
`--qualification` with this protocol file, restores the accepted reference
archive, materializes prepared and controlled generic builds, runs
prepare/run/verify once, and retains every receipt. Land the runner,
panel hash and workflow in an implementation PR before dispatch. Preserve
partial attempts. Update the research note and
`docs/index-calculus-scoreboard.html` from stored measurements in the same
evidence PR that archives the result. Classification follows AGENTS.md §3;
stage diagnostics are never reported as end-to-end speedups.
