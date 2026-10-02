# GPU-wide exact-v3 hint queue: first panel

The frozen 128-thread resolver experiment completed from clean source
`7e5281e6b276b89428dbb74a1c05a50ca5b61a02` on one NVIDIA RTX PRO6000 Blackwell
Server Edition, CUDA 13.3.73. Its terminal producer decision is
**DO_NOT_PROMOTE**. The independent native post-run audit passes; see
`INDEPENDENT-AUDIT.md` and `results/independent-audit.json` (SHA-256
`948903619feb3c92f99aac4afd20265d6a52a903ec8b37278c912eceb491d5dc`).

| Complete-update benchmark arm | Median B updates/s | Ratio to matched control | Correctness |
|:--|--:|--:|:--|
| Exact-v3 block queue | 7.021655 | 1.000000 | PASS |
| GPU-wide queue,128-thread resolver | 5.610974 | 0.799018 geometric mean | PASS |

Five paired ratios are 0.798980,0.798811,0.798973,0.799153 and0.799175.
A/A maximum symmetric drift is 0.056254%; every candidate pair loses about20%.
Each ranked row completes 201,863,462,912 scalar updates, including counter
reset, raw selection, full cold resolution, hot arithmetic, host launch and
synchronization costs. The goal of 26 B/s and the 1.10 engineering gate both fail.
This table compares the same exact-v3 point map. The separate sigma-fused
selected preset remains faster at about15.4–15.6 B/s in its matched panels;
these are whole-walk kernel measurements, not an end-to-end challenge solution.

All preflight gates passed before any timing:

- 149,940 native ownership/phase model cases and native field/reference controls.
- 49 production selector/resolver device cases with forced empty/full/sparse
  histories, canaries, unique owners, partial blocks and repeated resets.
- 300/300 initial reference reports per arm, zero drops, odd 95-step launches.
-Independent spread 300 replay per arm including at least 299 nonzero trails.
-Byte-identical sorted full corpora containing 1,710,384 records per arm,
  plus matching 511/513-worker partial corpora.
-Identical prefix checkpoints/corpora, four matching continuation checkpoints
  and sorted continuation corpora across both source and destination modes.

The resource split is measured, not inferred:

| Kernel | Registers/thread | Local bytes/thread | Static shared bytes/block |
|:--|--:|--:|--:|
| Control hot/block kernel |128|400|1040|
| Candidate hot kernel |126|0|0|
| Candidate selector |56|0|0|
| Candidate resolver |124|400|0|

Each table stage also uses 57,052 dynamic shared bytes. The resolver's initial
188 blocks of128 threads expose only four warps per SM. This is a plausible
scheduling limitation to test separately; it is not established as the cause
of the full regression. A different resolver geometry remains a follow-up,
without changing this negative result or the frozen admission thresholds.

Raw provenance and compact evidence are in `results/`. The immutable archive
has SHA-256 `f677fa729e662a7718ed0beb1d7880473697bb25588d54f19f93c6f39e5caa6f`
and 81,726,971 bytes, retained in the Modal `ecc2k130-jobs` volume under token
`afd2fd1cde744544b68e96328f276dec`. The authenticated retrieval command is:

```sh
ECC_GPU=RTX-PRO-6000 modal run modal_job.py \
  --fetch-token afd2fd1cde744544b68e96328f276dec --out /tmp/global-hints-replay
```

The first local GPU-name override failed before any function spawn or GPU
kernel and remains in `attempt-0-dispatch-failure.json`. The measured run
completed exit 0. Its original collector encountered transient ConnectionError;
an explicit volume fetch recovered the completed archive without another GPU
run. The scientific manifest verifies 79 immutable files; a post-run manifest
also binds the manifest itself and the wrapper's final job.log and exit-code
(82 files). All 24 recorded source hashes match the measured checkout.

The frozen `results/launch.json` records the earlier explicit volume fetch.
The original collector later received the same completed
archive and wrote a supplementary `collector-completion-receipt.json` with a
different completion timestamp. Both identify the same function call, token,
clean source and exit code. The independent audit binds the later collector
receipt; the root rerun also passes using the earlier fetch receipt. Their
scientific audit fields are identical; both immutable receipts are retained.

The later source-manifest glob correction only prevents future attempts from
trying to hash the newly added results directory. It changes no measured
source, binary, sample or decision input. The default mode remains off.

The added Linux checker compilation exposed a GCC range-loop-construction
warning promoted to an error. Making that local string construction explicit
repairs portability without changing any audit predicate. The repaired checker
again passes the full sealed-archive audit; all scientific fields match the
frozen audit. The original source/audit hashes remain evidence, and the additive
`results/gcc-portability-repair.json` binds the repaired source and failed CI
job. This repair requires no additional GPU measurement.
