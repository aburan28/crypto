# PKM tower oracle, round 3: a compact basis, and `m = 4` at `N = 16`

Raw data for §12 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
change, the identity check, the cell, the predictions and the decision rule
(§§12.1–12.6) before any round-3 cell ran, and §12.3 discloses the one
measurement of a new-size system made before that. The engine is
`src/cryptanalysis/f4_fp_tower.rs`, and the systems come from
`examples/pkm_tower_pilot.rs` with `--engine tower`: the pilot's and round 2's
systems wherever the flags coincide (same seed rule).

Everything here is a stage diagnostic on the solver axis (`AGENTS.md` §2). No
row prices a decomposition end to end, so there is no `S` and no scoreboard row.

## Checks

```sh
# the compact-basis build against round 2: every row and every trace step
python3 research/pkm_tower_round3_20260925/compare_builds.py
```

The CI workflow `.github/workflows/pkm-tower-pilot.yml` runs it, with the
pilot's and round 2's checks and the engine's unit tests.

## The change (note §12.1)

The basis keeps the rows the elimination leaves, dense over the step's columns
without a divisor from the lead on, at 4 bytes an entry, instead of sparse
polynomials at 20 bytes a term. Nothing else in the algorithm changes. The
example now prints each step's trace as the step ends (`--trace`), with the
kept elements, their entries, `B'`'s entries and the process's memory.

## The profiling run (the disclosure of note §12.3)

`profile/profile-kummer-m4-p1-N16.log` is the stderr of one run of round 2's
D2 system (Kummer, `m = 4`, `p₁ = 2013265921`, `N = 16`, target 0: flags
`--engine tower --p 2013265921 --kinds kummer --controls tower --planted 0
--random 1 --ladder-t none --m 4 --t-min 4 --max-t-m4 4 --budget 2700
--max-nnz 2000000000 --trace`, environment `F4T_DIAG=1`, 14 GB address-space
cap). It ran on the round-2 engine with the temporary counters of
`profile/profiling.patch` (a diff against the round-2 `f4_fp_tower.rs`, sha256
`04317c7f439992f349996ca7b0f74ea4b57444efe1dcd2551af2c9e3909e4d0e`), from 07:46:45
to 08:12:55 UTC on 2026-09-25, and was stopped by hand inside step 49 (exit 143).
Its `[diag]` lines give, per step, the basis's elements and terms, the retirable
ones, `B'`'s entries, the echelon's entries and the process's memory.

## The identity check (note §12.2)

`replay_round2.sh` re-ran, with the new build and round 2's flags, every system
that round 2 and its cross-check finished. Its rows and logs are in `replay/`,
named after the round-2 files they repeat, and `progress.txt` records when each
run started and ended. `compare_builds.py` compares them with round 2's.

It ran from 08:20 to 09:04 UTC and, after a restart of the machine, from 14:11
to 14:13 UTC on 2026-09-25. The machine restarted during the confirmation run
`C-K1-kummer-m2-p1-t11`, which had written no row. Its partial log is kept as
`replay/C-K1-kummer-m2-p1-t11.interrupted.log`. The script was started again,
skipping the runs that had ended, and re-ran that one from the start.

| part | round-2 files | rows | identical |
|:--|:--|--:|--:|
| cross-check | XV1–XV5 | 308 | 308 |
| round 2's cells, to the sizes they finished | K1, K1n, K0, I0, M4 (`N ≤ 12`), M3 (`N ≤ 15`) | 36 | 36 |
| confirmations | the five `C-*` cells | 10 | 10 |
| diagnostic | D1 | 1 | 1 |
| **total** | | **355** | **355** |

Every field of every row is identical except the wall clock (`ms`, `wall_s`),
and no round-2 row of a replayed file is missing. That includes the
multiply-add counts, the nonzeros and the residue sizes. **106 trace steps are
identical** on every field round 2 printed, the times aside:
- D1 (29 steps);
- the replayed M4 system at `N = 12`, target 0, against D1 (29);
- K0's target 0 against both of round 2's profiling traces of it (P0, 2 × 24).

**Paired runs.** `pair_memory.py` ran four systems that round 2 finished, each
with round 2's build (the baseline) and then this round's (the candidate), one
after the other, from 14:15 to 14:48 UTC. It recorded each process's peak
resident memory (`ru_maxrss` through `wait4`) and its wall time. The rows,
logs and summary are in `pairs/`, and the two rows of every pair agree on every
field but the wall clock.

| system (target 0) | peak, baseline | peak, candidate | candidate / baseline | wall, baseline | wall, candidate |
|:--|--:|--:|--:|--:|--:|
| K0: Kummer, `m = 2`, `p₀`, `N = 20` | 668 MB | 463 MB | 0.69 | 62 s | 59 s |
| D1: Kummer, `m = 4`, `p₁`, `N = 12` | 338 MB | 171 MB | 0.51 | 11 s | 10 s |
| M3: Kummer, `m = 3`, `p₁`, `N = 15` | 1,738 MB | 1,084 MB | 0.62 | 380 s | 374 s |
| K1: Kummer, `m = 2`, `p₁`, `N = 22` | 2,213 MB | 1,523 MB | 0.69 | 585 s | 572 s |

One run per build is not a timing claim, and wall time is a practicality note
(`AGENTS.md` §2). Every multiply-add is the same in both builds.

## Builds

The replay and the paired runs use the example built (`cargo build --release
--example pkm_tower_pilot`, stable Rust 1.94.1, x86-64 with AVX-512) from the
sources committed with the pre-registration:

| file | sha256 |
|:--|:--|
| `src/cryptanalysis/f4_fp_tower.rs` | `c5253bcb4296246eb443f0be6851511eeece1e731d80fcc8afcb7630460a747a` |
| `examples/pkm_tower_pilot.rs` | `51eb7fb3d09d068b015211d34f7d839f5fb61f45bf3ce2730b737318c30f6934` |
| `src/cryptanalysis/mod.rs` | `30e44b6f9facd323b75e8cb4943274f4bc63ec4983e40e5b4727fd0b5d36d1b1` |
| the binary | `ae387b056c74f35bcc18bb7b6ae4b2d2d1d1fdf9ba79f4b81f6aeb7acbc430bb` |

The baseline of the paired runs is the example built the same way from the
merge of round 2 (`8f7589f2`), whose engine, example, module list and lock file
are those of `main` when this round began (`d9243908`): binary sha256
`623655a36224096fa1f30d221e57cea99cb3b960eee2710236b5468ee57125a1`.
