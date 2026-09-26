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
