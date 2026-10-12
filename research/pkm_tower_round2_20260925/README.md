# PKM tower oracle, round 2: the sparse tower engine and `N` past the pilot

Raw data for §11 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
cells, predictions and decision rule (§§11.4–11.6) before any round-2 cell ran,
and §11.3 discloses the one new-size measurement made before that. The engine
is `src/cryptanalysis/f4_fp_tower.rs`; the systems come from
`examples/pkm_tower_pilot.rs` with `--engine tower`, and are the pilot's systems
wherever the flags coincide (same seed rule).

Everything here is a stage diagnostic on the solver axis (`AGENTS.md` §2): no
row prices a decomposition end to end, so there is no `S` and no scoreboard
row.

## Checks

The pilot's scripts read these files as they read the pilot's:

```sh
# every finished tower row against exhaustive search over V^m
python3 research/pkm_tower_pilot_20260924/verify.py research/pkm_tower_round2_20260925/runs/*.jsonl
# tables, repeats, planted solutions, and f4_fp against f4_fp_tower system by system
python3 research/pkm_tower_pilot_20260924/analyze.py \
    research/pkm_tower_pilot_20260924/runs/*.jsonl research/pkm_tower_round2_20260925/runs/*.jsonl
```

`analyze.py` exits non-zero if a repeated system disagrees with itself, a
planted solution was lost, or the two engines disagree on whether a system has a
solution. The CI workflow `.github/workflows/pkm-tower-pilot.yml` runs both, and
the engine's unit tests.

## Runs

`p₀ = 786433 = 3·2^18 + 1` (the pilot's prime) and
`p₁ = 2013265921 = 15·2^27 + 1`. Every run except the `P0` profiling runs used a
14 GB address-space cap (`ulimit -v 14000000`). A run's `.log` holds the example's stderr: one line per
system, and the step trace where `--trace` was given.

**Cross-check (engine validation, note §11.2).** The pilot's cells are re-run
with `--engine tower` and the pilot's flags, some widened to every kind and
control, with no size above the pilot's largest for its `m`. 38 systems have no
`f4_fp` counterpart, from the widened cells and from one cell whose targets were
drawn differently (note §11.2).

| file | flags (all with `--engine tower`) |
|:--|:--|
| `XV1-m2-tower-null` | `--kinds kummer,dickson,isogeny --m 2 --controls tower,null --max-t 9 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none` |
| `XV2-m3-tower` | `--kinds kummer,dickson,isogeny --m 3 --controls tower --max-t-m3 4 --planted 2 --random 2 --planted-max-t 3 --budget 900 --ladder-t none` |
| `XV3-m4-tower-random` | `--kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 0 --random 2 --budget 900 --ladder-t none` |
| `XV3b-m4-tower-planted-stop64` | `--kinds kummer,isogeny --m 4 --controls tower --max-t-m4 3 --planted 2 --random 0 --stop-below 64 --budget 900 --ladder-t none` |
| `XV4-ladder` | `--kinds kummer,isogeny --m none --planted 2 --random 0 --budget 300 --ladder-t 4,6 --ladder-g 1,2,4,8,12` |
| `XV5-kummer-m2-p32` | `--p 3221225473 --kinds kummer --m 2 --controls tower --max-t 8 --planted 2 --random 2 --planted-max-t 7 --budget 900 --ladder-t none` |

**The disclosure of note §11.3.** `P0-profiling-kummer-m2-p0-N20` holds the
two rows and step traces of the one new-size measurement made while tuning the
engine: the Kummer `m = 2` cell at `N = 20`, `p₀`, random target 0
(`--kinds kummer --m 2 --controls tower --t-min 10 --max-t 10 --planted 0
--random 1 --budget 3000 --ladder-t none --trace`). They come from two
intermediate builds of the engine that were not committed. The first formed
every row in memory; the second forms rows on demand, as the committed engine
does. Their traces agree line for line. Cell K0 re-measures the same system with
the committed build, and `analyze.py` compares the copies.

**Round 2 (note §11.4).** `run_round2.sh` ran the cells in the note's order,
one system at a time, from 01:54 to 04:32 UTC on 2026-09-25, with
the example built from the commit that pre-registered them.

| file | cell | what happened |
|:--|:--|:--|
| `K1-kummer-m2-p1` | Kummer, `m = 2`, `p₁`, `N = 12–22` | Complete: 12 systems, all refuted. `D = 5` up to `N = 18`, 6 at 20 and 22 |
| `K1n-kummer-m2-p1-null` | its null control, `N = 12–20` | Complete: 10 systems, all refuted, with the tower's `D` at every `N` |
| `K0-kummer-m2-p0-N20` | Kummer, `m = 2`, `p₀`, `N = 20`, staircase stop at 8 | Complete: both systems refuted at `D = 6`, so the stop never fired. Target 0 repeats `P0` exactly |
| `I0-isogeny-m2-p0-N20` | isogeny, `m = 2`, `p₀`, `N = 20`, staircase stop at 8 | Complete: both refuted at `D = 6` |
| `M4-kummer-m4-p1` | Kummer, `m = 4`, `p₁`, `N = 8, 12, 16` | `N = 8, 12` complete (`D = 6, 7`). The first `N = 16` system ran out of memory after 29 minutes ("memory allocation of 131072 bytes failed" in the log). The process ended there: no row, and the second `N = 16` target never ran |
| `M3-kummer-m3-p1` | Kummer, `m = 3`, `p₁`, `N = 9–18` | `N = 9, 12, 15` complete (`D = 6, 6, 7`). Both `N = 18` systems stopped for size, the `--max-nnz` cap of `2·10⁹`, after about 44 minutes each. They had reached degree 7 at 104,600 columns (`D ≥ 7`; rows with `"oversize": true`), and with every system at that `t` stopped, the example ended the cell |

The exit codes in `progress.txt` are all 0 and carry no information:
`run_round2.sh` wrote them with `$(date …)` in the same line, which resets `$?`
first. The logs show how each run ended. `confirm_round2.sh` reads the status
before anything else runs, and so does `run_round2.sh` since the round (a
logging fix, made after its cells ran).

**Confirmations (note §11.4).** `confirm_round2.sh` re-runs, once and whole,
every `m = 2` cell with a finished system at `D ≥ 6`, with the degree bound at 5
(`--cap 5`). It lists them in `confirm_cells.txt`. `analyze.py` keeps these
rows out of the tables (their degree bound is below the default) and reports
them in their own section.

They ran after the D2 diagnostic, from 04:57 to 05:01 UTC. The chain meant to
start them straight after the cells got "Permission denied", because the script
had been committed without its executable bit.

| file | cell re-run with `--cap 5` | what happened |
|:--|:--|:--|
| `C-K1-kummer-m2-p1-t10` | K1, `N = 20` | Both systems ended with 30,335 pairs above the bound, not refuted: confirmed |
| `C-K1-kummer-m2-p1-t11` | K1, `N = 22` | The same: 30,335 pairs above the bound, confirmed |
| `C-K1n-kummer-m2-p1-null-t10` | K1n, `N = 20` | The same: confirmed |
| `C-K0-kummer-m2-p0-N20-t10` | K0 (`--stop-below 8` as in the cell) | The same, with no stop: confirmed |
| `C-I0-isogeny-m2-p0-N20-t10` | I0 (`--stop-below 8`) | 14,828 pairs above the bound, not refuted, no stop: confirmed |

**Diagnostics (after the cells, labelled as such).**

| file | flags (all with `--engine tower --p 2013265921 --kinds kummer --controls tower --planted 0 --random 1 --ladder-t none --trace`) | what it is |
|:--|:--|:--|
| `D1-kummer-m4-p1-N12-trace` | `--m 4 --t-min 3 --max-t-m4 3 --budget 600` | The step trace of M4's `N = 12` target 0, which it repeats exactly. It shows how `m = 4` refutes, through a long tail of degree-6 steps (note §11.7) |
| `D2-kummer-m4-p1-N16-trace` | `--m 4 --t-min 4 --max-t-m4 4 --budget 1500 --max-nnz 2000000000` | M4's first `N = 16` system, which ran out of memory after 29 minutes, re-run with a 25-minute budget so that the engine stops itself and prints its trace. It stopped at the budget after 48 steps, all at degree 7 or below and every one from step 12 on at 7, without refuting: `D ≥ 7`. Its row says `timed_out` |

## Builds

The cross-check runs used the example built (`cargo build --release --example
pkm_tower_pilot`, stable Rust 1.94.1, x86-64 with AVX-512) from the sources
committed with the pre-registration:

| file | sha256 |
|:--|:--|
| `src/cryptanalysis/f4_fp_tower.rs` | `c54ee3025b708b006e946f467fa15a146ff578d7881cbd7a31f8452ccf586f56` |
| `examples/pkm_tower_pilot.rs` | `dbd4b72f780b2665d85bab158f120d17ed9aa0a0ac47225aba40e33fb91b30d3` |
| `src/cryptanalysis/mod.rs` | `c059b21b4bc4ecfc8a4cf650894641daa4969789bc98b34d9259987b3cb2b5fa` |
| the binary | `02656958303f6d9e7e383619ff2d0a325380d57f1e0c0738379be76da99973b1` |

Before the round-2 cells, `origin/main` was merged into the branch (merge
commit `d8b1df8f`). It touched none of these sources except the module list in
`mod.rs` (sha256 `30e44b6f9facd323b75e8cb4943274f4bc63ec4983e40e5b4727fd0b5d36d1b1`
after the merge). The round-2 cells, the confirmations and the diagnostics ran
on that build, binary sha256
`0ef59debbfc54d89373bd6ddb8296d3263042d5b32aae67e762a03980b289751`. It
reproduces `XV2` and `XV5` exactly (60 systems, every deterministic field). The
commit with the results changes only comments in `f4_fp_tower.rs` and
`pkm_tower_pilot.rs`.

`xv-progress.txt` and `progress.txt` record when each run started and ended.
