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

## Stage 1: the sizing run (note §16.8)

`M3b-sizing-t0`, 02:02–02:42 UTC on 2026-10-06, exit 0. Its row equals round
2's target-0 row on every compared field, multiply-adds included
(`pkm_tower_check compare-rows`, note §16.8). Its trace, from
`pkm_tower_check trace-table runs/M3b-sizing-t0.log` (memory in MB of 2^20
bytes):

| step | degree | S-rows | pivot rows | columns | without a divisor | residues | new elements (lowest degree) | `B'` entries | kept entries | multiply-adds | s | memory MB |
|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|
| 1 | 4 | 1 | 2 | 15 | 13 | 1 | 1 (3) | 25 | 35 | 25 | 0.0 | 3 |
| 2 | 4 | 3 | 4 | 24 | 20 | 3 | 2 (4) | 72 | 74 | 178 | 0.0 | 3 |
| 3 | 5 | 9 | 17 | 99 | 82 | 5 | 5 (4) | 1203 | 452 | 2996 | 0.0 | 3 |
| 4 | 5 | 24 | 48 | 262 | 214 | 17 | 12 (4) | 8536 | 2816 | 24925 | 0.0 | 3 |
| 5 | 5 | 29 | 84 | 451 | 367 | 26 | 16 (4) | 25325 | 8117 | 86122 | 0.0 | 3 |
| 6 | 5 | 33 | 97 | 556 | 459 | 25 | 19 (5) | 37002 | 16455 | 200476 | 0.0 | 4 |
| 7 | 6 | 288 | 517 | 1940 | 1423 | 231 | 90 (4) | 626555 | 131130 | 12162295 | 0.0 | 8 |
| 8 | 5 | 11 | 227 | 1282 | 1055 | 11 | 10 (4) | 192326 | 140869 | 1243714 | 0.0 | 8 |
| 9 | 5 | 11 | 232 | 1258 | 1026 | 11 | 11 (4) | 193502 | 151212 | 1336568 | 0.0 | 8 |
| 10 | 5 | 11 | 229 | 1130 | 901 | 10 | 8 (5) | 174842 | 158059 | 1266094 | 0.0 | 8 |
| 11 | 6 | 841 | 1451 | 4706 | 3255 | 797 | 212 (5) | 3905562 | 720167 | 312931510 | 0.1 | 26 |
| 12 | 6 | 1367 | 2286 | 7064 | 4778 | 1334 | 334 (5) | 9246184 | 2105680 | 2049605052 | 0.3 | 53 |
| 13 | 6 | 1331 | 2856 | 7852 | 4996 | 1331 | 343 (5) | 12416872 | 3646942 | 3049171370 | 0.4 | 72 |
| 14 | 6 | 1077 | 4750 | 15432 | 10682 | 1075 | 584 (5) | 40351764 | 8845242 | 8636336218 | 1.2 | 204 |
| 15 | 6 | 918 | 5486 | 15916 | 10430 | 893 | 406 (5) | 46903432 | 12603954 | 12017448139 | 1.5 | 254 |
| 16 | 6 | 417 | 5719 | 16312 | 10593 | 393 | 243 (5) | 49909297 | 14933233 | 8907618119 | 1.2 | 266 |
| 17 | 6 | 62 | 5951 | 16419 | 10468 | 62 | 39 (5) | 52039204 | 15310171 | 6857994705 | 1.0 | 279 |
| 18 | 6 | 16 | 5907 | 16401 | 10494 | 15 | 15 (6) | 51860617 | 15458618 | 6320858989 | 0.9 | 283 |
| 19 | 7 | 22318 | 22000 | 45792 | 23792 | 20738 | 2354 (5) | 452256413 | 62880699 | 1378038867950 | 150.8 | 2004 |
| 20 | 6 | 23 | 8291 | 26660 | 18369 | 17 | 17 (5) | 126322407 | 63163033 | 46236043808 | 4.7 | 2016 |
| 21 | 6 | 7 | 98 | 227 | 129 | 1 | 1 (6) | 11402 | 63163162 | 47004 | 0.0 | 2016 |
| 22 | 7 | 29963 | 44892 | 88958 | 44066 | 29951 | 4343 (6) | 1624882277 | 212726463 | 8249356789831 | 912.6 | 7181 |
| 23 | 7 | 33337 | 54274 | 98183 | 43909 | 33328 | 5453 (6) | 1978877004 | 407681007 | 11469056026077 | 1272.0 | 9223 |
| 24 (stopped) | 7 | 28321 | 64289 | 104600 | 40311 | — | — | 2227368915 | 407681007 | — | 20.3 | 9300 |

The rule of note §16.3 (`pkm_tower_check size-cap runs/M3b-sizing-t0.log`):

```
stopped step 24: E = 2227368915 entries of B', q = 40311 columns without a divisor, live heap H = 1583 MB = 1659895808 bytes
A = 14336000000 bytes, reserve 2.5e9 bytes: C = 2500000000 (C >= E: the step fits by the rule)
```
