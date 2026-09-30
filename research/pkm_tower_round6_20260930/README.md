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
