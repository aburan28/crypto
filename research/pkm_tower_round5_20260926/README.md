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

It ran from 07:16 to 07:42 UTC on 2026-09-26, with the example built from the
pre-registration commit (hashes under "Builds").

| files | what happened |
|:--|:--|
| `runs/GD-*` | all 12 systems refuted. At `m = 3`, `D_sig` = 6 at `N = 9` (F4: 6) and 7 at `N = 12` (F4: 6): **the gate fails** |
| `runs/GP-*` | 10 systems refuted. At `m = 3`, `N = 9`, `D_sig` = 6 (F4: 6). Both `N = 12` systems reached the 600 s budget at `D_sig ≥ 7`: **the gate fails** |
| `runs/S4-*` | round 4's engine, rebuilt: all 12 rows repeat round 4's G rows on every deterministic field and signature counter. The new field `lm_degree_max` (`D_lm`) shows round 4's extra degree was, at every size but `m = 2`, `N = 12`, a larger leading-monomial ideal, not bookkeeping |

`compare_gate.py` prints the full table (verdicts, `D`, `D_sig`, `D_lm`,
widths, multiply-adds, zero rows, peak memory) and applies the gate.
`verify.py` agrees with all 34 finished rows.

Two shared scripts changed with this round, so that two signature variants of
one system count as two measurements:
- `research/pkm_tower_pilot_20260924/verify.py`, which checked one row per
  system and engine;
- `research/pkm_tower_pilot_20260924/analyze.py`, which compared the variants
  as repeats.

A signature row's engine now includes its module order, rewrite order and
steps. Rows written before the steps field existed count as round 4's variant.
The earlier rounds' output is unchanged, apart from the heading of the
signature engine's tables.

## Builds

The example was built (`cargo build --release --example pkm_tower_pilot`,
stable Rust 1.94.1, x86-64 with AVX-512, 4 cores, 15 GB) from the
pre-registration commit `718d91d2`:

| file | sha256 |
|:--|:--|
| `src/cryptanalysis/sig_fp_tower.rs` | `cd2ea69a7b6faffd5ee2500a71e249e2378363e52e524a639880c8430489c691` |
| `src/cryptanalysis/f4_fp_tower.rs` | `6fc2408987bc049c3390a9c61b24e7c3a1b8244ca284cb9dd6b19ee245b5526d` |
| `examples/pkm_tower_pilot.rs` | `d790573864bb266288bcc7dd35c3de77c41381e997b48aa0590a9e2cec446f9b` |
| `src/cryptanalysis/mod.rs` | `a2efafa0bc5658e9c377d98f853f4b028d335c06535c8eede384ae93512d897a` |
| `run_gate.py` | `8d32431d67b956d7faa083acf46ee6da508bb0bf15ffee35f61f4e9a04261bcf` |
| the binary | `70b3d80575c5bd29c90f42577911fe56bdd44de7ad55fe7bff28806059491184` |
