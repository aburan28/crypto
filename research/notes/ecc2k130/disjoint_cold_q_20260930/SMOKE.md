# Local integration smoke before hosted timing

The input freeze is at commit `63a695d7`. It contains 45 Q/label pairs, 15,390 unique new Q at their respective degrees, and [INPUT_RECEIPT.json](INPUT_RECEIPT.json), which independently verifies all 15,390 equations `[d]G=Q` and disjointness from the 31-file prior inventory. The exact frozen source materialized with all twenty SHA-256 pins; an offline Cargo release build produced the two examples. This is a build and correctness smoke on macOS, **not** an eligible CPU comparison.

The same frozen binaries then ran one complete three-arm block at n37/L1 with the off prefilter and one at n41/L1024 with the blocked prefilter. Both processes in each compact pair reached full rank, recovered every public Q, and passed independent base, rank-row, witness and scalar replay. The rho arm independently recovered the same Q. Both verifier receipts are `SMOKE_PASS` and `timing_eligible=false`.

| Smoke cell | Public Q | Compact ranks replayed | Compact target logs checked | Rho logs checked | Receipt SHA-256 |
|:--|--:|--:|--:|--:|:--|
| n37/L1 block 0 | 1 | 2 | 2 | 1 | `3e762608de4d469090b57b24ed13a5e48432cd53da8c25d7c969b12a60558e80` |
| n41/L1024 block 0 | 1,024 | 2 | 2,048 | 1,024 | `f6e566f874dca7d22553a39cf9e669e48e701dc37aa0e5f1249834099e20cf78` |

[LOCAL_SMOKE.tar.gz](LOCAL_SMOKE.tar.gz) holds the 30 raw stdout, stderr, base, rank, target, materialization, run and receipt files from these two blocks. [LOCAL_SMOKE_FILES.json](LOCAL_SMOKE_FILES.json) pins each member's byte count and SHA-256, plus archive SHA-256 `e6327d4a0959d5bb59bd862f39f7b202f7780a9f607319942692499b6e5f1bea`. The archive used sorted two-component member names, fixed metadata and gzip mtime zero. CI rehashes every member, extracts it under the runner temporary directory, and runs the same independent verifier in relocated mode. The original local command used `run_cold.py --mode smoke` with `--source-root`, `--compact`, `--rho` and `--materialization` pointing to the materialized frozen tree, followed by `verify_cold.py --cell <cell> --run-dir <raw-directory> --out <receipt>`; both commands and every environment value are also in each `cold_run.json`.

No hosted timing arm had run at the time of this smoke. The six-cell paired CPU decision remains pending the isolated Linux workflow, an all-raw second-host replay, and outcome PR.
