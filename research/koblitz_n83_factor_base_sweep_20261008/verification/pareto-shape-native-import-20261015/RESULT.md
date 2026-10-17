# Native import panel for the retained N83 source Pareto shapes

Status: `PASS_NATIVE_IMPORT_PANEL_ONLY`. The frozen [plan](plan.json)
(SHA-256 `ced6a67a75de8ee785d34cbd9aa9e79d9767683e84c6390d4f115e39896d20ba`)
selects one representative for each of the 24 distinct source-metric shapes
in the retained N83 Pareto receipt. Those representatives use 18 unique
signed-Frobenius factor-base archives, totaling 27,822,670 compressed bytes.
All 18 S3 HEAD and GET pairs and all 24 native worker imports exited 0 within
the 1,200-second panel cap. Every download matched its manifest size and
compressed BLAKE3 identity. The 60 measured stage calls sum to 232.496 seconds
of L0 wall time; the import calls account for 107.455 seconds of that sum.
The maximum recorded native import RSS was 168,902,656 bytes.

The panel reused the executable with SHA-256
`c5b9af9930f293533de819c0a673a0b7149ad2418800fa3fee7509d772d207e7`
from the earlier clean Rust 1.98.1 build. It pinned the worker source, full
Rust `src/` tree, pilot manifest, source Pareto receipt, and exact S3 URIs
before launch. The source commit was `97b3cde9ef7ca9d5fb8b7329cb47cd5d3ba1c4c0`.
`result.json` records every S3 HEAD response, compressed-object SHA-256,
command, cap, exit status, wall time, RSS and raw-stream hash. The byte-preserved
120 stdout/stderr streams are in `raw-logs.tar.gz`; their SHA-256 digests are
also in `result.json` and the copied archive was checked stream by stream.
`artifact_manifest.json` hashes every retained output.

A separate [row-level audit](audit.py) parsed every downloaded archive,
checked each signed-orbit label and point uniqueness, and recomputed all 24
full-orbit overlaps and residual x-domain partitions. It also rehashed the
source commit, binary, result, S3 objects and raw worker streams. Its
`audit.json` status is `PASS_S3_BYTES_RAW_STAGES_AND_24_ORBIT_PARTITIONS`;
the bounded audit call exited 0 in 9.813 L0 seconds, with 357,810,176-byte
maximum RSS. The K256 hash-base/K600 sequential-pool representative that
source-dominates the earlier K600 Gray control was among the passing imports:
the native receipt and independent replay both found zero shared and 600
residual signed orbits.

The retained package was replayed again by extracting `raw-logs.tar.gz` next
to a copy of `result.json` and linking the same verified local S3 downloads.
The fresh audit output was byte-identical to `audit.json`; the bounded check
passed in 5.249 L0 seconds. `retained_replay.json` records that packaging
check and the source, tar and audit digests.
An isolated copy with `arm-00.stdout` changed from `PASS` to `FAIL` was then
rejected before an audit receipt could be written; `tamper_control.json` and
`tamper.stderr` retain that fail-closed control.

This panel measures import compatibility and cost on one macOS host. Its
representatives are 24 metric shapes drawn from 311 tied or nondominated
source pairs, not an exhaustive runtime sweep. All import calls ran in a
fixed order after common S3 downloads, so their timings include cache and
order effects. No SAT relation search, useful rank, linear algebra, fresh
target-log stage or matched rho ran here. `selected_best_total_runtime` is
null in the plan, worker result and independent audit; no runtime-winning
factor base follows from these receipts. The local Docker daemon still
returns HTTP 503, which prevents the frozen Linux paired search from running.

The archives remain in S3 under
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/`.
No new S3 object was written. The 18 local copies in
`/private/tmp/n83-pareto-shape-native-import-20261015/archives/` are scratch;
the manifest identities and S3 response metadata in `result.json` let a
reviewer fetch and replay them again. Extract `raw-logs.tar.gz` as `logs/`
beside `result.json` before rerunning the audit against those archive bytes.
