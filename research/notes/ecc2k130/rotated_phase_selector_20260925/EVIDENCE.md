# Frozen n19 phase-selector evidence

`receipt.json` records the exact frozen source/input hashes, command lines, UTC intervals, exit codes, child wall/CPU/high-water RSS, stream hashes, and hashes of every evidence file. `producer/outcome.json` contains both fixed phase orders, all cap counts and probe totals, the frozen 64-case row certificates, and the preregistered priority-screen verdict. Each `producer/first_hit_*.u8` is exactly 130,873 bytes in target-index order `Q=[j]H`, `0<=j<q`. A byte 0–18 is the first **position** in that order, whose phase exponent is in `outcome.json`; 255 means all 19 phases missed. The independent `verification.json` replays both complete arrays and every fixed-case row. The four stdout/stderr files are retained even though they are empty.

The input archive is the committed `../rotated_pdp_corpus_20260925/evidence/raw.tar.gz`, SHA-256 `39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`. The fixed point-only target file is `../rotated_four_base_joint_rank_20260925/inputs/point_only.json`, SHA-256 `3b9cadf3513eb1e504d7d517acda57db0ea6d35698439a79a842444165cae10e`. The complete source/input hash map is in `FROZEN.json`, SHA-256 `351486905c13a3a422a7adbaf4ab099032d68556bf62c94a59f23e55610eff14`.

From the repository root, with Python 3.13:

```sh
python3 research/notes/ecc2k130/rotated_phase_selector_20260925/ci_replay.py
python3 research/notes/ecc2k130/rotated_phase_selector_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_phase_selector_20260925/evidence
```

The first command is hash-only and reads no phase outcome. The second checks the committed receipt and independently rebuilds all q phase decisions and fixed-case group rows. `run.py --out NEW_EMPTY_PATH` can reproduce the bounded producer/verifier measurement; the archived result is the single preregistered local run and no input/order/cap was changed after it.
