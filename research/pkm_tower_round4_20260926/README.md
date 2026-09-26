# PKM tower oracle, round 4: signature criteria (F5/GVW), checked against F4

Raw data for §13 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
engine, the cells, the predictions and the decision rule (§§13.1–13.6) before
any round-4 cell ran, and §13.3 discloses the runs made while the engine was
built. The signature engine is `src/cryptanalysis/sig_fp_tower.rs`, the F4
engine `src/cryptanalysis/f4_fp_tower.rs`, and the systems come from
`examples/pkm_tower_pilot.rs` with `--engine sig` and `--engine tower`: round
2's systems (same seed rule).

Everything here is a stage diagnostic on the solver axis (`AGENTS.md` §2). No
row prices a decomposition end to end, so there is no `S` and no scoreboard row.

## Checks

```sh
# the signature rows beside F4's, this round's F4 rows against round 2's,
# and the predictions and adoption rule of note sections 13.5-13.6
python3 research/pkm_tower_round4_20260926/compare_engines.py
# every finished round-4 row, both engines, against exhaustive search over V^m
python3 research/pkm_tower_pilot_20260924/verify.py research/pkm_tower_round4_20260926/runs/*.jsonl
# tables and repeats, with the earlier rounds
python3 research/pkm_tower_pilot_20260924/analyze.py \
    research/pkm_tower_pilot_20260924/runs/*.jsonl \
    research/pkm_tower_round2_20260925/runs/*.jsonl \
    research/pkm_tower_round3_20260925/runs/*.jsonl \
    research/pkm_tower_round4_20260926/runs/*.jsonl
# the engine's unit tests
cargo test --release --lib sig_fp_tower
```

The CI workflow `.github/workflows/pkm-tower-pilot.yml` runs all of them, with
the earlier rounds' checks.

## Round 4 (note §13.4)

`run_round4.py` runs each size as one process per engine, F4 (cells T2–T4)
first and then the signature engine (cells G2–G4), each measuring the size's
two targets. Every process runs under the 14 GB address-space cap, and its
peak resident memory, wall time and exit status go to `memory.jsonl`.
`progress.txt` records when each process started and ended.

```sh
cargo build --release --example pkm_tower_pilot
python3 research/pkm_tower_round4_20260926/run_round4.py
```

It ran from 05:22 to 06:34 UTC on 2026-09-26, with the example built from the
pre-registration commit (hashes under "Builds").

| files | what happened |
|:--|:--|
| `runs/T2-*`, `runs/T3-*`, `runs/T4-*` | F4 on all 10 sizes, both targets: 20 rows, each identical to round 2's except the wall clock |
| `runs/G2-*` (`N` = 12–20) | the signature engine refutes all 10 systems, with `D_sig` = 6, 6, 6, 7, 7 against F4's 5, 5, 5, 5, 6 |
| `runs/G3-*` (`N` = 9, 12) | refutes all 4, with `D_sig` = 7, 8 against 6, 6 |
| `runs/G3-kummer-m3-p1-N15` | ran out of the 14 GB cap after 1,381 s, on target 0, in the split of 25.4 million pending pairs after step 47 (exit −6, no row). The `.log` holds the 47 steps' traces and the allocation failure. The `m = 3` sizes stop there, as fixed in advance |
| `runs/G4-*` (`N` = 8, 12) | refutes all 4, with `D_sig` = 7, 9 against 6, 7 |

`compare_engines.py` prints the full table: widths, multiply-adds, rows built,
zero rows and peak memory. It also evaluates the note's predictions (all four
confirmed) and the adoption rule (not met, so F4 stays the measuring engine).
`verify.py` agrees with all 38 finished rows.

Each `.log` holds one trace line per step, printed as the step ended. For the
signature engine that line gives:
- the step's signature degree and position, and its pairs;
- the rows the criteria skipped;
- the rows built, by kind;
- the columns and nonzeros;
- the zero, singular and new rows;
- the late pairs;
- the times and multiply-adds;
- the process's resident memory and its peak so far.

## Builds

The round's example was built (`cargo build --release --example
pkm_tower_pilot`, stable Rust 1.94.1, x86-64 with AVX-512, 4 cores, 15 GB) from
the pre-registration commit `bf41cee2`:

| file | sha256 |
|:--|:--|
| `src/cryptanalysis/sig_fp_tower.rs` | `1061c62e896e2ed52f1c1c6b08f7fe78f98bc5d8539cbddcb1cca2b349c7d335` |
| `src/cryptanalysis/f4_fp_tower.rs` | `6fc2408987bc049c3390a9c61b24e7c3a1b8244ca284cb9dd6b19ee245b5526d` |
| `examples/pkm_tower_pilot.rs` | `2689e1519cc22785cba91032d253a46b27a0616fc621a5654c67bd583b16b170` |
| `src/cryptanalysis/mod.rs` | `7c0f7dd57646a48d7ca506107381dfbb444f4f2c804f264f96d78f51ce830156` |
| `run_round4.py` | `a0799d2df81e49d58fa7c66c0d3ffa5cb4688c0aa8322cb1b8f171e5178be552` |
| the binary | `2e1c91f3a6b72149254cd181e1b6c66d8782ba7e2e121e75839ee9f3ff45fd8a` |

The next commit, `6117ec4e`, made two changes for CI's newer toolchain, and
neither changes a row:
- it sorts a step's pairs with `sort_by_key` on the same key, the same stable
  sort;
- it moves one line of `mod.rs` into rustfmt's order.

The example built from that commit reproduces the round's G4 `N = 8` and G2
`N = 14` rows, on both targets, on every field except the wall clock.

After the runs, `compare_engines.py` gained the listing of processes that
ended without writing every row.
