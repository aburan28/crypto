# PKM tower oracle, round 7: `m = 3` at `N = 18`

Raw data for §16 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
stages, the rule for the size cap, the predictions and the decision rule
(§§16.1–16.6) before any of these runs, and §16.4 discloses what was run before
them.

The cell is round 2's M3 at `t = 6`: Kummer, `m = 3`, `p₁ = 2013265921`,
`N = 18`, its two random targets. Round 2 stopped both systems for size in step
24 (`D ≥ 7`).

Everything here is a stage diagnostic on the solver axis (`AGENTS.md` §2). No
row prices a decomposition end to end, so there is no `S` and no scoreboard row.

## Runs

```sh
cargo build --release --example pkm_tower_pilot
research/pkm_tower_round7_20261006/run.sh stage1
research/pkm_tower_round7_20261006/run.sh stage2 C      # C by note section 16.3's rule
research/pkm_tower_round7_20261006/run.sh confirm C D T # target T finished with D
```

`run.sh` holds every flag. Every run is under `ulimit -v 14000000` (KiB), at
the host's four threads.
- **Stage 1, `M3b-sizing-t0`.** Target 0 runs with round 2's flags
  (`--max-nnz 2000000000`, `--budget 7200`).
  - It must stop for size in step 24 and reproduce round 2's row.
  - Its `tower stop in step 24` line gives `B'`'s entries, the columns without
    a divisor and the live heap, from which `C` follows.
- **Stage 2, `M3b-t0`, `M3b-t1`.** Each target runs as its own process
  (`--targets`), with `--max-dense C` and `--budget 43200`.
- **Confirmation, `C-M3b-t<T>-cap<D−1>`.** A system that finished with `D` is
  re-run with the degree bound at `D − 1`.

`runs/<name>.jsonl` holds the rows. `runs/<name>.log` holds the step trace, one
line per step with the process's memory, and the stop line if a stop cut a step
off. `progress.txt` records each run's start, end and exit status.

## Build

The engine and the example are as committed with the pre-registration.

| file | sha256 |
|:--|:--|
| `src/cryptanalysis/f4_fp_tower.rs` | `f1452c496e2b198886d5382d56b2fbe17d78c07b3e4da226a27324e45245d517` |
| `examples/pkm_tower_pilot.rs` | `5511895760ffccff97cb4390d8f0e96eec1127751f9edf6bcc56b777009faa6f` |
| `run.sh` | `18713435e43c7f85d0ef1c653f8a340b2b8e6697e99d99678f9753fd3777346e` |
| the binary | `297f63dd3dcc5abec8131623bed8cf6c886fc00aefbb1cac6f060425ba064d66` |

## Host

Intel Xeon at 2.10 GHz, 4 cores without SMT, AVX-512 (F, BW, CD, DQ, VL, IFMA,
VBMI, VBMI2, VNNI, BF16, FP16, BITALG, VPOPCNTDQ) and AVX2, 16 GB without swap.
Linux 6.18.44, x86-64, glibc 2.39, rustc 1.94.1. It is the host class of round
6, not that of rounds 2–5 (2.80 GHz), so no wall time here is compared with
theirs.

## Checks

The checks run on a native checker, `examples/pkm_tower_check.rs`, which
replaces the Python scripts `verify.py` and `analyze.py` (`AGENTS.md`, "no
Python"). Note §16.3 says what it must reproduce before it reads this round.
