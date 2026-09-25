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
`p₁ = 2013265921 = 15·2^27 + 1`. Every run used a 14 GB address-space cap
(`ulimit -v 14000000`). A run's `.log` holds the example's stderr: one line per
system, and the step trace where `--trace` was given.

**Cross-check (engine validation, note §11.2).** The pilot's cells, re-run with
`--engine tower` and the pilot's flags, so every system here is one the pilot
measured with `f4_fp`.

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

**Round 2 (note §11.4).** `run_round2.sh` runs the cells in the note's order,
one system at a time, with the example built from the commit that
pre-registered them.

| file | cell |
|:--|:--|
| `K1-kummer-m2-p1` | Kummer, `m = 2`, `p₁`, `N = 12–22` |
| `K1n-kummer-m2-p1-null` | its null control, `N = 12–20` |
| `K0-kummer-m2-p0-N20` | Kummer, `m = 2`, `p₀`, `N = 20`, staircase stop at 8 |
| `I0-isogeny-m2-p0-N20` | isogeny, `m = 2`, `p₀`, `N = 20`, staircase stop at 8 |
| `M4-kummer-m4-p1` | Kummer, `m = 4`, `p₁`, `N = 8, 12, 16` |
| `M3-kummer-m3-p1` | Kummer, `m = 3`, `p₁`, `N = 9–18` |

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

The round-2 cells use the same build. `xv-progress.txt` and `progress.txt`
record when each run started and ended.
