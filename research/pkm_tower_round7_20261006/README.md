# PKM tower oracle, round 7: `m = 3` at `N = 18`

Raw data for §16 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
stages, the rule for the size cap, the predictions and the decision rule
(§§16.1–16.6) before any of these runs, and §16.4 discloses what was run before
them.

The cell is round 2's M3 at `t = 6`: Kummer, `m = 3`, `p₁ = 2013265921`,
`N = 18`, its two random targets. Round 2 stopped both systems for size in step
24 (`D ≥ 7`).

**Result: `D = 7` on both targets**, each refuted in step 28 and confirmed by a
re-run with the degree bound at 6 (note §§16.9–16.11).

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

## Stage 2: the measurement (note §16.9)

Each target runs as its own process: `run.sh stage2 2500000000`, so `--max-dense
2500000000` and `--budget 43200`.

**A restart.** The session's container was reclaimed during the first attempt
at target 0, in step 23, between 03:00 and 03:59 UTC. The run ended without an
exit status. Its files are `runs/M3b-t0.restart1.{log,jsonl}`, and
`progress.txt` records the restart. The run was repeated whole with the same
flags and binary. Its first 22 steps equal the killed run's
(`pkm_tower_check compare-trace`).

### Target 0: `D = 7`, refuted

`M3b-t0`, 03:59–06:56 UTC, exit 0. Its first 23 steps equal stage 1's. Its
trace (`pkm_tower_check trace-table runs/M3b-t0.log`):

| step | degree | S-rows | pivot rows | columns | without a divisor | residues | new elements (lowest degree) | `B'` entries | kept entries | multiply-adds | s | memory MB |
|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|
| 1 | 4 | 1 | 2 | 15 | 13 | 1 | 1 (3) | 25 | 35 | 25 | 0.0 | 3 |
| 2 | 4 | 3 | 4 | 24 | 20 | 3 | 2 (4) | 72 | 74 | 178 | 0.0 | 3 |
| 3 | 5 | 9 | 17 | 99 | 82 | 5 | 5 (4) | 1203 | 452 | 2996 | 0.0 | 3 |
| 4 | 5 | 24 | 48 | 262 | 214 | 17 | 12 (4) | 8536 | 2816 | 24925 | 0.0 | 3 |
| 5 | 5 | 29 | 84 | 451 | 367 | 26 | 16 (4) | 25325 | 8117 | 86122 | 0.0 | 4 |
| 6 | 5 | 33 | 97 | 556 | 459 | 25 | 19 (5) | 37002 | 16455 | 200476 | 0.0 | 4 |
| 7 | 6 | 288 | 517 | 1940 | 1423 | 231 | 90 (4) | 626555 | 131130 | 12162295 | 0.0 | 8 |
| 8 | 5 | 11 | 227 | 1282 | 1055 | 11 | 10 (4) | 192326 | 140869 | 1243714 | 0.0 | 8 |
| 9 | 5 | 11 | 232 | 1258 | 1026 | 11 | 11 (4) | 193502 | 151212 | 1336568 | 0.0 | 8 |
| 10 | 5 | 11 | 229 | 1130 | 901 | 10 | 8 (5) | 174842 | 158059 | 1266094 | 0.0 | 8 |
| 11 | 6 | 841 | 1451 | 4706 | 3255 | 797 | 212 (5) | 3905562 | 720167 | 312931510 | 0.1 | 26 |
| 12 | 6 | 1367 | 2286 | 7064 | 4778 | 1334 | 334 (5) | 9246184 | 2105680 | 2049605052 | 0.4 | 54 |
| 13 | 6 | 1331 | 2856 | 7852 | 4996 | 1331 | 343 (5) | 12416872 | 3646942 | 3049171370 | 0.5 | 77 |
| 14 | 6 | 1077 | 4750 | 15432 | 10682 | 1075 | 584 (5) | 40351764 | 8845242 | 8636336218 | 1.6 | 207 |
| 15 | 6 | 918 | 5486 | 15916 | 10430 | 893 | 406 (5) | 46903432 | 12603954 | 12017448139 | 1.7 | 250 |
| 16 | 6 | 417 | 5719 | 16312 | 10593 | 393 | 243 (5) | 49909297 | 14933233 | 8907618119 | 1.4 | 278 |
| 17 | 6 | 62 | 5951 | 16419 | 10468 | 62 | 39 (5) | 52039204 | 15310171 | 6857994705 | 1.0 | 284 |
| 18 | 6 | 16 | 5907 | 16401 | 10494 | 15 | 15 (6) | 51860617 | 15458618 | 6320858989 | 0.9 | 295 |
| 19 | 7 | 22318 | 22000 | 45792 | 23792 | 20738 | 2354 (5) | 452256413 | 62880699 | 1378038867950 | 202.3 | 1998 |
| 20 | 6 | 23 | 8291 | 26660 | 18369 | 17 | 17 (5) | 126322407 | 63163033 | 46236043808 | 8.4 | 2006 |
| 21 | 6 | 7 | 98 | 227 | 129 | 1 | 1 (6) | 11402 | 63163162 | 47004 | 0.0 | 2006 |
| 22 | 7 | 29963 | 44892 | 88958 | 44066 | 29951 | 4343 (6) | 1624882277 | 212726463 | 8249356789831 | 989.9 | 7225 |
| 23 | 7 | 33337 | 54274 | 98183 | 43909 | 33328 | 5453 (6) | 1978877004 | 407681007 | 11469056026077 | 1399.1 | 9343 |
| 24 | 7 | 28321 | 64289 | 104600 | 40311 | 28248 | 5118 (6) | 2227368915 | 586673662 | 11126680454368 | 1387.8 | 10998 |
| 25 | 7 | 26804 | 69642 | 105985 | 36343 | 26804 | 7048 (6) | 2253999269 | 808782146 | 13010822579882 | 1572.4 | 11965 |
| 26 | 7 | 64767 | 76826 | 106288 | 29462 | 64767 | 14093 (5) | 2145052750 | 1124685222 | 36338267877533 | 4052.4 | 12898 |
| 27 | 6 | 25042 | 33262 | 48633 | 15371 | 25042 | 14356 (3) | 505282469 | 1242307743 | 8238985427103 | 966.2 | 12896 |
| 28 | 4 | 3097 | 4213 | 5228 | 1015 | 1015 | 1015 (0) | 4271270 | 1242307743 | 7073729193 | 0.8 | 12794 |

The row: `D = 7`, refuted in step 28, 28 steps, widest step 106,288 columns and
141,593 rows, at most 2,344,243,162 nonzeros and 2,253,999,269 `B'` entries,
89,912,686,130,244 multiply-adds, 10,587 s at four threads. The address space
reached 13,411,820 KiB of the 14,000,000 KiB cap, a figure read from `/proc`
during step 27.

**Confirmation.** `C-M3b-t0-cap6` (`run.sh confirm 2500000000 7 0`), 06:56–06:57
UTC, exit 0. With the degree bound at 6 it ends after 18 steps with 22,382
pairs above the bound, without refuting and without a staircase stop. Its 18
steps equal the full run's.

### Target 1: `D = 7`, refuted, with target 0's trace

**Two more restarts.** The container was reclaimed during each of the first
two attempts at target 1. Neither ended with an exit status, and
`progress.txt` records both restarts.

| attempt | ran (UTC) | how it ended | log ends at | files |
|:--|:--|:--|:--|:--|
| 1 | 10-06 06:56; last seen alive 08:21, in step 25 | reclaimed; container back 10-07 00:09 | step 24 | `runs/M3b-t1.restart1.{log,jsonl}` |
| 2 | 10-07 00:12; last seen alive 01:06, in step 24 | reclaimed; container back 14:04 | step 23 | `runs/M3b-t1.restart2.{log,jsonl}` |
| 3 | 10-07 14:06–17:15 | exit 0 | step 28 (`1`) | `runs/M3b-t1.{log,jsonl}` |

Each attempt repeats the run whole with the same flags and binary. Attempt 3's
first 24 steps equal attempt 1's, and its first 23 attempt 2's, multiply-adds
included (`pkm_tower_check compare-trace`).

**Host.** Attempts 2 and 3 ran on another host class: an Intel Xeon at
2.80 GHz, 4 cores without SMT, AVX-512 (F, BW, CD, DQ, VL, VNNI) and AVX2,
16 GB without swap, Linux 6.18.44. The binary is the one above (same sha256).
It is built for baseline x86-64, and the engine's runtime dispatch (`avx512f`,
then `avx2`) takes the AVX-512F path on both hosts. No time of target 1's is
compared with target 0's.

**The trace is target 0's.** Attempt 3's 28 steps equal target 0's in every
count but two (`pkm_tower_check compare-trace runs/M3b-t0.log runs/M3b-t1.log
--except muladds` leaves the second):
- the multiply-adds of steps 22 and 26, 2,097 and 25,511 more than target 0's;
- step 26's nonzeros, 2,344,243,163 against 2,344,243,162.

Both counts skip zero entries, so an entry that vanishes by chance in one
target and not the other moves them. That is no surprise: step 26 has about as
many entries as `F_p` has elements. The degrees, pairs, columns, residues, new
elements, `B'` sizes and kept entries are equal in every step.

The row: `D = 7`, refuted in step 28, 28 steps, widest step 106,288 columns and
141,593 rows, at most 2,344,243,163 nonzeros and 2,253,999,269 `B'` entries,
89,912,686,157,852 multiply-adds, 11,349 s at four threads. The trace's peak
resident memory is 12,933 MB. The address space reached 13,426,460 KiB of the
14,000,000 KiB cap, read from `/proc` during step 27.

**Confirmation.** `C-M3b-t1-cap6` (`run.sh confirm 2500000000 7 1`), 17:15–17:16
UTC on 2026-10-07, exit 0. With the degree bound at 6 it ends after 18 steps
with 22,382 pairs above the bound, without refuting and without a staircase
stop. Its 18 steps equal the full run's.

## Checks

From the repository root, with `C=target/release/examples/pkm_tower_check` and
`R=research/pkm_tower_round7_20261006/runs`:

| command | result |
|:--|:--|
| `$C verify $R/*.jsonl` | 2 rows checked against exhaustive search over `V³` (262,144 triples each), 2 agree; the two capped rows have no verdict to check |
| `$C compare-trace $R/M3b-sizing-t0.log $R/M3b-t0.log` | 23 steps, 23 identical |
| `$C compare-trace $R/M3b-t0.restart1.log $R/M3b-t0.log` | 22 steps, 22 identical |
| `$C compare-trace $R/M3b-t0.log $R/C-M3b-t0-cap6.log` | 18 steps, 18 identical |
| `$C compare-trace $R/M3b-t1.restart1.log $R/M3b-t1.log` | 24 steps, 24 identical |
| `$C compare-trace $R/M3b-t1.restart2.log $R/M3b-t1.log` | 23 steps, 23 identical |
| `$C compare-trace $R/M3b-t1.log $R/C-M3b-t1-cap6.log` | 18 steps, 18 identical |
| `$C compare-trace $R/M3b-t0.log $R/M3b-t1.log --except muladds` | 28 steps, 27 identical: step 26's nonzeros differ by one |
| `$C analyze` over rounds 1–5 and 7 (CI's command) | no repeated measurement disagrees on a deterministic field; the `m = 3` line at `p₁` reads `D = 6, 6, 7, 7` and A1 inconclusive; both confirmation runs confirm |
