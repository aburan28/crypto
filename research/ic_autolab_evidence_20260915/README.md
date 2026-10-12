# AutoLab pull-request evidence

This PR adds fast Koblitz scalar multiplication, fast odd-degree factor-base point lifting, the bounded AutoLab/Harbor harness, and the recorded benchmark results.

## Result and scope

The final frozen campaign measured **47.72% fewer complete cold-worker instructions** against the prior fast-lifting batch1/window8 winner. Its 60 fresh targets, 1,356 profile/native pairs, independent audit and fresh-process replay passed. The earlier campaign measured 26.87% against its original batch16 baseline; the unsuccessful configuration-only campaign and all development variants are retained. These percentages use different target sets and are not combined into a cumulative measurement.

The numerical claims bind the archived source snapshots, compiler, fixtures and worker binaries. The arithmetic changes are also ported onto main at `90bb7a02` in this PR. Integration tests exercise that port separately; they do not assign the historical percentage to a different source tree. The worker emits the original field-certificate JSON using explicit `degree` and `low_terms` fields, so it builds without unrelated serialization changes.

The metric is Valgrind 3.22 amd64 Ir for the whole cold worker, including setup, failed attempts, verification, matrix work, descent, final recovery and cleanup. This is implementation engineering on small Koblitz cases. The final candidate uses 4.98 times the matched rho instructions. Native times remain diagnostics.

## Restore and independently audit

Run from the repository root with Python 3.12 and zlib 1.3, the recorded environment:

```bash
python3 research/ic_autolab_evidence_20260915/evidence.py
ROUND=research/ic_autolab_next30_20260915/runs/autolab-20260915T224645Z/jobs/autolab-20260915T224645Z/task__4brp5nH/verifier/round-20260915t224645z
python3 "$ROUND/evaluator/tournament.py" verify --round "$ROUND"
```

Use `--only round-20260915t224645z` to restore just the final campaign, or `--destination /tmp/ic-review` for an isolated extraction. Each bundle restores paths under `research/`. Existing identical files are accepted; differing files are preserved and cause an error. Frozen reports and original machine paths remain unchanged as provenance.

The raw Callgrind intervals were already individually gzip-compressed. Simply archiving them produced hundreds of megabytes. The evidence format stores their decoded payloads together, plus each original gzip header/trailer and SHA-256. Restoration recreates the exact original gzip bytes with zlib before accepting each file. Files that could not be recreated exactly are stored verbatim. Every archived file has a size, mode and original SHA-256; each outer archive has another SHA-256 in [bundles.json](bundles.json). Use the restoration tool, not a plain tar extraction. The development archive is stored in two parts to fit GitHub's upload limit; the tool joins them automatically, checks each part, then verifies the original archive hash before extraction. A different zlib output is rejected, never silently accepted.

All raw profile/native outputs, job inputs, receipts, source manifests, complete frozen sources, measured binaries, evaluator copies, stage summaries and decisions are included. Build caches are omitted. Full original working directories remain preserved on the experiment host. The bundles also retain unsuccessful development probes and the historical pricing catalog.

## Validation

- [PR arithmetic tests](validation/pr-tests.json): nine focused Rust tests, including the newly landed degree-31 divisor-base regression; 16 harness tests and 15 tournament tests.
- [PR worker certificates](validation/pr-worker-validation.json): complete IC and rho recovery for all 60 recorded confirmation targets, independently checked; [raw outputs](validation/pr-worker-certificates.tar.xz).
- Ten evidence-restoration tests cover byte identity, archive corruption, reconstructed-byte mismatch, existing-file preservation, path traversal, missing payloads and multipart integrity.
- The final campaign is re-audited after restoration to verify that the portable bundle contains the full evidence needed by the original evaluator.

## Recorded harness and new runs

The frozen harness scripts are preserved as used. Their controller records the experiment host's AutoLab installation at `/home/ubuntu/crypto-work/autolab-upstream-20260915`; update `UPSTREAM` for another installation. AutoLab commit `4127da3dde8449be61a1cf9859473b9fbbd51751` used Harbor 0.3.0. The final run used the nop agent with a locked candidate, 2 CPUs, 8 GiB, a 2-hour supervisor cap and a 55-minute verifier cap.

The image recipe under `task/environment/Dockerfile` expects the recorded Rust distribution, a locked Cargo vendor directory, the optional OpenCode executable, baseline source, evaluator and harness in its generated build context. The local image IDs in launch records identify the actual experiment images; they are not public registry tags. Restoring and auditing the saved evidence needs no Docker image, network access or paid service. Future experiments must build their own image, lock a new candidate and use fresh excluded-target sets before making another measured claim.

The final instruction comparison used the previous winner as baseline, with unchanged batch1/window8, sparse linear algebra, pair-table solver, m=3 and a 4,096-trial cap. Only the promotion threshold changed to 0.70. [Proposal](../ic_autolab_next30_20260915/proposal.json), [result](../ic_autolab_next30_20260915/README.md), [costs](../ic_autolab_next30_20260915/COSTS.md).

The reported 0.4662 container CPU-hours cover the additional-30% experiment, excluding subsequent PR preparation and integration tests. Local compute and this assistant session have no supplied dollar price.
