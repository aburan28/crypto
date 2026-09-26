# PKM tower oracle, the gate of §13.10: signature criteria with F4's steps

Raw data for §14 of
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md`. The note fixes the
variant, the gate, the context runs and the prediction (§§14.1–14.4) before
any of these runs, and §14.4 discloses the one plumbing check made before them.
The engine is `src/cryptanalysis/sig_fp_tower.rs` with
`Steps::PolynomialDegree` (`--sig-steps degree` in
`examples/pkm_tower_pilot.rs`). F4's rows are round 4's T cells
(`research/pkm_tower_round4_20260926/`), which repeat round 2's.

Everything here is a stage diagnostic on the solver axis (`AGENTS.md` §2). No
row prices a decomposition end to end, so there is no `S` and no scoreboard row.

## Checks

```sh
# the gate of note section 14.2, and the context of 14.3
python3 research/pkm_tower_round5_20260926/compare_gate.py
# every finished row against exhaustive search over V^m
python3 research/pkm_tower_pilot_20260924/verify.py research/pkm_tower_round5_20260926/runs/*.jsonl
# the engine's unit tests
cargo test --release --lib sig_fp_tower
```

## The gate (note §14.2)

`run_gate.py` runs three signature-engine variants, one process per variant and
size, each measuring the size's two targets:

| cell | steps | module order | role |
|:--|:--|:--|:--|
| GP | by polynomial degree | position over term | the gate |
| GD | by polynomial degree | signature degree first | the gate |
| S4 | by signature degree (round 4's engine) | position over term | context: its `D_lm` |

The gate's systems (`m = 3`, `N = 9` and 12) run first, then the context sizes
(`m = 2`, `N = 12`, 14, 16; `m = 4`, `N = 8`). The budget is 600 s per system,
with the 14 GB address-space cap. Each process's peak resident memory and exit
status go to `memory.jsonl`, and `progress.txt` records when it started and
ended.

```sh
cargo build --release --example pkm_tower_pilot
python3 research/pkm_tower_round5_20260926/run_gate.py
```
