# PKM tower pilot, 2026-09-24

The raw data behind §10 of
[`RESEARCH_PKM_TOWER_ORACLE.md`](../notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md):
a design-validation pilot of that note's Stage A. It measures how the F4
solving degree of Petit–Kosters–Messeng tower systems grows with the number of
tower unknowns `N = m·t`. It is a stage diagnostic of the solver axis
(`AGENTS.md` §2). It prices no pipeline and reports no `S`.

## Reproduce

```bash
cargo build --release --example pkm_tower_pilot
B=target/release/examples/pkm_tower_pilot
$B <flags from the table below> --out runs/<name>.jsonl 2> runs/<name>.log

python3 research/pkm_tower_pilot_20260924/analyze.py research/pkm_tower_pilot_20260924/runs/*.jsonl
python3 research/pkm_tower_pilot_20260924/verify.py  research/pkm_tower_pilot_20260924/runs/*.jsonl
```

- `analyze.py` prints one table per (kind, `m`, control), the fit that §5.2 of
  the note pre-registered, the slopes that replace its degenerate interval, and
  the generator-count ladder.
- `verify.py` recounts every finished tower row by exhaustive search over
  `V^m`, without F4. It exits non-zero on any disagreement.

Every cell draws its tower, curve and targets from the seed `0x504B4D54`
(`1347112276`) combined with (kind, `m`, `t`, `g`). A cell is therefore the same
system in every run that contains it. Both scripts count such a system once.
`analyze.py` also checks that every repeat agrees on every deterministic field
and prints any that do not. Rows from the sparse tower engine (round 2,
`research/pkm_tower_round2_20260925/`) carry `"engine": "f4_fp_tower"`. The
scripts count one system measured by both engines once per engine, and
`analyze.py` compares the engines system by system. It fails if they disagree on
whether a system has a solution.

## Runs

`p = 786433 = 3·2^18 + 1` unless stated. The runs queued in the background had a
13 GB address-space cap (`ulimit -v 13000000`). Budgets are per system.

| file | flags | build | what happened |
|:--|:--|:--|:--|
| `smoke.jsonl` | `--kinds kummer --max-t 4 --max-t-m3 3 --planted 2 --random 2 --budget 20 --ladder-t 4 --ladder-g 1,4` | b0 | First smoke run. Written before `solving_degree_max` existed; `analyze.py` skips its rows. Kept as the record of the accounting correction (note §10.1) |
| `smoke2.jsonl` | `--kinds kummer --max-t 6 --max-t-m3 3 --planted 2 --random 2 --budget 60 --ladder-t 4 --ladder-g 1,4` | b1 | Complete |
| `A1-kummer-m2-tower` | `--kinds kummer --m 2 --controls tower --max-t 9 --planted 2 --random 2 --budget 300 --ladder-t none` | b1 | Stopped by hand (SIGTERM) in its first `N = 16` planted system. That run was still certifying its basis at about 10 GB. Rows through `N = 14` are complete; `N ≥ 16` was re-run without planted targets as B1 |
| `A2-kummer-m3-tower` | `--kinds kummer --m 3 --controls tower --max-t-m3 5 --planted 2 --random 2 --budget 300 --ladder-t none` | b1 | Ended without a message in its first `N = 15` planted system (no panic in the log; killed for memory). Rows through `N = 12` are complete; two `N = 12` planted systems hit the 300 s budget in their certification tails. `N = 15` re-run without planted targets as C1 |
| `B1-kummer-m2-tower-large` | `--kinds kummer --m 2 --controls tower --t-min 8 --max-t 9 --planted 0 --random 3 --budget 900 --ladder-t none` | b2 | Complete |
| `A3-dickson-m2-tower` | `--kinds dickson --m 2 --controls tower --max-t 9 --planted 2 --random 2 --budget 300 --ladder-t none --planted-max-t 7` | b2 | Complete |
| `A4-isogeny-m2-tower` | `--kinds isogeny --m 2 --controls tower --max-t 8` + the A3 settings | b2 | Complete |
| `A5-kummer-m2-controls` | `--kinds kummer --m 2 --controls null,naive --max-t 9` + the A3 settings | b2 | Aborted (allocation failure) in its first null `N = 18` system. Null rows through `N = 16` are complete; the naive cells never started. Re-run as F1 and F2 |
| `A6-ladder` | `--kinds kummer,isogeny --m none --planted 2 --random 0 --budget 300 --ladder-t 6 --ladder-g 1,2,4,8,12` | b2 | Complete |
| `A7-m3-dickson-isogeny` | `--kinds dickson,isogeny --m 3 --controls tower --max-t-m3 4 --planted 2 --random 2 --budget 300 --ladder-t none --planted-max-t 3` | b2 | Complete |
| `C1-kummer-m3-tower-large` | `--kinds kummer --m 3 --controls tower --t-min 5 --max-t-m3 5 --planted 0 --random 2 --budget 1800 --ladder-t none` | b2 | Aborted (allocation failure) in its first `N = 15` system after about 3.5 minutes. No row: `m = 3` at `N = 15` is beyond the dense `f4_fp` in 13 GB. Its `.jsonl` is empty. Re-run with the step trace as C1b |
| `C2-isogeny-m2-tower-large` | `--kinds isogeny --m 2 --controls tower --t-min 9 --max-t 9 --planted 0 --random 2 --budget 1800 --ladder-t none` | b2 | The first system finished (`D = 5`). The second was aborted (allocation failure) after 19 minutes, with no row. Re-run with the step trace and the staircase stop as C2b |
| `D1-kummer-m2-tower-p32` | `--p 3221225473 --kinds kummer --m 2 --controls tower --max-t 8 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none` | b2 | Complete |
| `E1-kummer-m2-tower-p32-N20` | `F4_DEBUG=1`, `--p 3221225473 --kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --planted 0 --random 2 --budget 3600 --ladder-t none` | b2 | Aborted (allocation failure) in its first `N = 20` system after 4.5 minutes. No row; the log holds the step trace, 23 steps, none above degree 5 |
| `R1-repro-final-build` | `--kinds kummer,dickson,isogeny --m 2,3 --controls tower,null,naive --max-t 5 --max-t-m3 3 --planted 2 --random 2 --budget 60 --ladder-t 4 --ladder-g 1,4` | b3 | Complete. It re-measures cells of smoke2 and A1–A7 with the committed source, and adds Dickson and isogeny nulls up to `N = 10` |
| `R3-repro-m3-b4` | `--kinds kummer,dickson,isogeny --m 3 --controls tower --max-t-m3 3 --planted 2 --random 2 --budget 120 --ladder-t none` | b4 | Complete. It re-measures `m = 3` cells of A2, A7 and R1 with b4, after the chain was generalised to any `m` |
| `R2-null-N16-stop` | `--kinds kummer --m 2 --controls null --t-min 8 --max-t 8 --planted 0 --random 2 --stop-below 8 --budget 600 --ladder-t none` | b3 | Complete. It is the A5 null cell that timed out at 300 s, re-run with the staircase stop: F4 pinned that system to at most one point at degree 5 |
| `F1-kummer-m2-null-N18-stop` | `--kinds kummer --m 2 --controls null --t-min 9 --max-t 9 --planted 0 --random 2 --stop-below 8 --budget 900 --ladder-t none` | b3 | Complete; both systems refuted, so the stop never fired |
| `F2-kummer-m2-naive-N14` | `--kinds kummer --m 2 --controls naive --t-min 7 --max-t 7 --planted 2 --random 2 --budget 300 --ladder-t none` | b3 | Complete |
| `C1b-kummer-m3-N15-trace` | `F4_DEBUG=1`, `--kinds kummer --m 3 --controls tower --t-min 5 --max-t-m3 5 --planted 0 --random 1 --budget 1800 --ladder-t none` | b2 | Aborted (allocation failure) after 3 minutes. No row; the log holds the step trace, 19 steps, none above degree 6 |
| `m4-smoke` | `--kinds kummer --m 4 --controls tower --max-t-m4 2 --planted 2 --random 2 --budget 120 --ladder-t none` | b4 | Complete |
| `G1-m4-tower-random` | `--kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 0 --random 2 --budget 900 --ladder-t none` | b4 | Complete |
| `G2-m4-tower-planted-stop` | `--kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 2 --random 0 --stop-below 8 --budget 900 --ladder-t none` | b4 | Stopped by hand after its Kummer `N = 8` cell. At `m = 4` every ordering of a decomposition solves the chain, so a planted system has at least 24 solutions and a stop at 8 cannot fire. Re-run as G2b |
| `G2b-m4-tower-planted-stop64` | `--kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 2 --random 0 --stop-below 64 --budget 900 --ladder-t none` | b4 | Complete. The stop fired at 24–43 standard monomials, and every planted solution was kept |
| `C2b-isogeny-m2-N18-stop-trace` | `F4_DEBUG=1`, `--kinds isogeny --m 2 --controls tower --t-min 9 --max-t 9 --planted 0 --random 2 --stop-below 8 --budget 1800 --ladder-t none` | b4 | Complete. Both systems reach `D = 5` in about 71 s. The second, the one C2 lost, is decomposable: the stop fired at a two-point staircase, where C2 had spent 19 minutes certifying the basis |

A run's `.log`, where there is one, holds the example's stderr: one line per
system, and any panic. The smoke runs and R2 were run with stderr discarded or
on the terminal, and have none.

**Builds.** All five builds construct the systems the same way and use the same
seed rule. They differ in the flags and the output fields that were added over
the day.
- b0: before `solving_degree_max`.
- b1: adds `solving_degree_max`, `max_cols_to_solution` and `steps_to_solution`.
- b2: adds `--planted-max-t` and records `V` for isogeny towers.
- b3: adds `--stop-below`, `--cap` and `--dump`.
- b4 (the source committed with this directory): adds the `m ≥ 4` chain and
  `--max-t-m4`.

R1 re-measured the cells it shares with the b1 and b2 runs using b3, and every
one agreed on every deterministic field. b4 changes nothing for `m ≤ 3`.

**Timings.** `ms` and `wall_s` are wall-clock times on a 4-core container.
F4's elimination uses every core. Some cells overlapped other runs, builds and
diagnostics, so the timings are practicality notes, never the metric.

## Row fields

| field | meaning |
|:--|:--|
| `solving_degree_max` | highest degree of a step at which F4 learned something (a new basis element, or `1`): the solving degree of `docs/ic/FRAMEWORK.md` §3 |
| `last_productive_degree` | the degree of the last such step (`f4_fp`'s older `solving_degree`); lower than the above when F4 climbs and descends |
| `degree_reached` | the highest degree F4 processed, including steps that only certified the basis |
| `max_cols_to_solution`, `steps_to_solution` | the widest matrix and the step count up to the last productive step |
| `inconsistent` | F4 derived `1`: no solution on the grid |
| `staircase_at_stop` | set only with `--stop-below`: the number of standard monomials F4 stopped at, which bounds the number of solutions; the basis is then not certified complete |
| `planted_ok` | on planted targets, whether the planted solution satisfies F4's basis |
| `tower` | the tower's public data: `g`, `zeta` (Kummer), `u`, `zeta` (Dickson), or `V` itself with the auxiliary curve (isogeny) |
