# Retained dependency CI source-closure failure

PR #1362 at `30c84566f782e6503901d187348709dfaad05cc2` has four failed
legacy compact-orbit jobs. At 2026-10-05 03:06 UTC, its exact-head rollup had
six successful checks, four failures and seventeen queued checks. These failures
are not a passing dependency gate. No manual retry or scientific panel was
launched locally to replace them.

| Original job | Retained complete log | SHA-256 |
| --- | --- | --- |
| [S3 batch, n41](https://github.com/aburan28/crypto/actions/runs/37254746170/job/111589313455) | `ci-1362-paired-s3-n41-original.log` | `3ada8dc69732b5e6cad0303a6a786153badf96a5d0d89f4bafaecb6291f19528` |
| [S3 batch, n53](https://github.com/aburan28/crypto/actions/runs/37254746170/job/111589313679) | `ci-1362-paired-s3-n53-original.log` | `0ac6881982d95260c5d31faad55c846b01b56e7717d479be23b2faeb2dede897` |
| [Swap grid, n41](https://github.com/aburan28/crypto/actions/runs/37254746161/job/111589313330) | `ci-1362-paired-swap-n41-original.log` | `b8ff9c4486855490b2d480afab4a66eb7c2397468dd60fbb1e07ddc430dd010a` |
| [Swap grid, n53](https://github.com/aburan28/crypto/actions/runs/37254746161/job/111589313463) | `ci-1362-paired-swap-n53-original.log` | `2fd636bc37f8c899f79193e315fa2f68c09ec006e37a57bf9f7528310570ae89` |

The workflows `.github/workflows/compact-s3-batch.yml` and
`.github/workflows/compact-swap-quotient.yml` copy the checked-out crate into
`runner.temp/frozen-src`, then overlay a selected old pinned source subset. Their
result is a hybrid of current and historical modules. The original n41 S3 log
checked out PR merge `509268435509486297773e885598d7c2e90d4137`.

All four original logs fail at compilation, before any measurement panel:

- Current embedded-source dependencies are absent from the sparse/copied tree:
  `docs/ic/boundary_targets.json`, `docs/ecbench/schema.sql` and
  `docs/curves/registry.json`.
- The current `prepared_n17_target` imports `IndividualLogObserver` and
  `IndividualLogInterruption`, introduced by PR #1362. The overlaid historical
  `koblitz_index_calculus` lacks those types. This new import exposes the hybrid
  dependency defect; the failures cannot be described solely as old missing files.

These are build failures, not zero-yield results, successful solver comparisons,
or evidence that the new runtime was admitted. No measurement artifact was
produced. Original manifests, workflow outputs and source snapshots stay intact.
Earlier PR #1342 also retains four legacy missing-source failures; this note
does not declare those resolved.

The repository requires native research execution. Repair requires a new
versioned path using a complete coherent native source snapshot and all embedded
data, with portable data-only replay where the question is historical custody.
Do not update a historical frozen manifest to include new code, silently compile
current-code fallbacks, extend/execute the legacy Python research harness, or
rerun a closed comparison panel just to turn CI green. That migration remains
pending and must preserve the original failure and its scientific boundaries.
