# PKM tower oracle, round 6: stop a step's elimination at a full echelon

Raw data for §15 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
change, the identity check, the timing protocol, the predictions and the
decision rule (§§15.1–15.6) before any of these runs. §15.4 discloses what was
read and run before them.
- The change is `full_rank_exit` in `src/cryptanalysis/f4_fp_tower.rs`, on by
  default, with `--no-full-rank-exit` in `examples/pkm_tower_pilot.rs`.
- Once a step's residue echelon has a pivot on every column, the step's
  remaining S-rows are skipped.
- The basis, the trace and `D` stay the same. Only the multiply-adds, the
  residual-row counts and the time may fall.

Everything here is a stage diagnostic on the solver axis (`AGENTS.md` §2). No
row prices a decomposition end to end, so there is no `S` and no scoreboard row.

## Checks

```sh
# the identity check of note section 15.2, with the multiply-adds saved
python3 research/pkm_tower_round6_20260930/compare_exit.py
# the timing protocol of section 15.3, summarized
python3 research/pkm_tower_round6_20260930/bench.py summary
# the engine's unit tests, the exit's two among them
cargo test --release --lib f4_fp_tower
```

## The identity check (note §15.2)

```sh
cargo build --release --example pkm_tower_pilot
research/pkm_tower_round6_20260930/replay.sh
```

`replay.sh` runs every system with the flags it was committed with, at the
host's default thread count. `progress.txt` records when each run started and
ended.
- **`off/`, the exit off:** round 2's cross-check (XV1–XV5), D1, and M4 at
  `N ≤ 12`. These must reproduce the committed rows exactly.
- **`replay/`, the exit on:** everything round 3 replayed (round 2's
  cross-check, cells, confirmations and D1), then round 3's `m = 4`, `N = 16`
  cells (M4b, C-M4b). Only the multiply-adds and the residual-row counts may
  differ, and only downwards.

## The timing protocol (note §15.3)

```sh
python3 research/pkm_tower_round6_20260930/bench.py run BASELINE_BIN CANDIDATE_BIN
```

`bench.py` times three systems, target 0 each, at one thread, every run
through `tools/isolated_bench.py` on CPU 3:
- first the baseline against a byte-identical copy of itself (A/A), five
  interleaved rounds;
- then the baseline against the candidate (A/B), five interleaved rounds.

The isolated-bench records are in `bench/records.jsonl`, the rows in
`bench/rows/` and the builds' hashes in `bench/builds.json`.

## Host

Intel Xeon at 2.10 GHz, 4 cores without SMT, AVX-512 (F, BW, CD, DQ, VL, IFMA,
VBMI, VNNI, BF16, FP16) and AVX2, 15 GB. Linux 6.18.44, x86-64, rustc 1.94.1.

The baseline is the example built from `main` at `9f4df528`, whose engine is
the one rounds 3–5 measured (`f4_fp_tower.rs` sha256
`6fc2408987bc049c3390a9c61b24e7c3a1b8244ca284cb9dd6b19ee245b5526d`). Its binary
sha256 is
`830cfad997e249499674f9e1b1ae923f0b83b7d41375ea2a998df944f3f40764`.

The candidate is the example built from the pre-registration commit
`00983dfc`:

| file | sha256 |
|:--|:--|
| `src/cryptanalysis/f4_fp_tower.rs` | `2080e3db300387add1c8faa3c9bf36b491dae295a72f40fd14d0086f6b6bd0f8` |
| `examples/pkm_tower_pilot.rs` | `e5e0efb7671ec223736a24dc94f952af3463237f979b3ad22ba5698c3c56cbf9` |
| `replay.sh` | `c08cc1f04df25beb2c263a992989cb33560aa25377763a8fed9af499da2191f7` |
| `bench.py` (as first run) | `2badaa01a114116b75db7256878b8d1afa9062dd1e8ef4a3a601a50e372a46bb` |
| the binary | `644ac8664e7369dccb11ede3fe158defbd295a05dbc27b8d87b18b1c31be140f` |

## Timing results (note §§15.7–15.9)

Every run was at one thread, pinned to CPU 3 by `tools/isolated_bench.py`, and
none was contended. The multiply-adds did not vary between repeats or between
the two copies of the baseline.

| system | multiply-adds, baseline → candidate | first measurement (5 pairs): B/A [95%] | second measurement (15 pairs): B/A [95%] | A/A (15 pairs) [95%] |
|:--|:--|:--|:--|:--|
| `m = 4`, `N = 12` | 85,477,602,933 → 51,426,567,808 (−39.8%) | 0.664 [0.649, 0.679] | **0.669 [0.663, 0.675]** | 1.001 [0.990, 1.012] |
| `m = 3`, `N = 12` | 10,595,733,281 → 10,595,732,141 | 1.023 [0.979, 1.067] | 1.019 [0.995, 1.043] | 0.999 [0.979, 1.018] |
| `m = 2`, `N = 16` | 6,178,203,100 → 5,655,592,357 (−8.5%) | 0.908 [0.850, 0.966] | **0.916 [0.900, 0.931]** | 1.006 [0.988, 1.025] |

- `bench/` holds the first measurement: 60 runs, 5 refused starts retried.
  - The literal timing clause of note §15.6 failed on it at `m = 3`: its "A/A
    spread", the difference of two five-run medians, came out at 0.19%.
  - The note's §15.7 records that. §15.8 then fixed a second measurement and a
    paired-interval test before running it.
- `bench2/` holds the second measurement: 180 runs, 18 refused starts retried.
  - The session's background limit stopped the driver once, inside the run
    `m4-N12-ab13-A`.
  - `bench.py` then gained a resume, and the run was repeated (note §15.9).

## Identity results (note §15.10)

`replay.sh` ran from 05:40 to 08:18 UTC on 2026-09-30, with the candidate built
from `00983dfc`, at four threads. `compare_exit.out` holds `compare_exit.py`'s
output, including the per-system table of multiply-adds saved.
- **`off/`:** 313 rows, exactly as committed, multiply-adds included; 129 trace
  steps identical.
- **`replay/`:** 359 rows as required; 1,111 trace steps as required, 37 of
  which skipped rows.
- **`verify.py`:** agrees with all 200 distinct rows it checks.
- **Multiply-adds:** they fall on exactly the 179 rows whose runs skipped
  S-rows. The largest savings are 39.2% at `m = 4`, `N = 12` and 24.6% at
  `m = 4`, `N = 16`.

`rows_skipped_full_rank` and the trace's "skipped" in these files run up to one
chunk short. The round's build set the count to the later chunks' rows instead
of adding them, losing the rest of the chunk in which the echelon filled. The
exact count is the step's S-rows minus the fill point, both printed in the
trace. The counter was fixed after the runs; nothing else was affected.
