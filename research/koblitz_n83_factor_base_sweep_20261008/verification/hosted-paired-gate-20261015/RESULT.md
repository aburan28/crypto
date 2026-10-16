# Hosted Linux fallback for the N83 K600 Gray/hash first trial

Status: `FIRST_HOSTED_RUN_PRODUCER_FAILURE_AWS_CREDENTIALS`; a credential-free
retry is prepared. [Run 38961528522](https://github.com/aburan28/crypto/actions/runs/38961528522)
checked out the frozen source and pulled the pinned ARM image, then failed in
`aws-actions/configure-aws-credentials@v4` with a request-signature mismatch.
S3 fetch, import, and both searches were skipped. The provider retained the
[setup artifact](https://github.com/aburan28/crypto/actions/runs/38961528522/artifacts/12072513691)
with the image ID `sha256:c1999f81da8e1f787006d58cbaece8ac4837d7ae04ad16784d829470f4bebc8d`.
The failed attempt provides no relation-yield or factor-base timing evidence.

The [original paired plan](../k600-gray-hash-paired-protocol-20261014/plan.json)
remains frozen at SHA-256
`fbd3fb943bc0f181797d106c41ee9961b8d70768b5c65a5584002e3c1b5e488f`.
This fallback retains its two S3 factor-base pairs, public target, first trial,
3+2 split, 5,000-conflict and 32-model caps, arm order and 4 GiB zero-swap
one-CPU search envelope. It uses the official Linux ARM `rust:1.98.1-bookworm`
image with immutable multi-architecture manifest digest
`sha256:93ce27a88655056a51dbdd8f5f2d7ddc071c7b0070fb288a37b5a285fc83971e`.
The image ID and cgroup readings are recorded at run time. This is a new host
protocol; its wall times must be compared within the hosted pair, not spliced
into the earlier local Linux or macOS timing series.

The revised [manual workflow](../../../../.github/workflows/n83-hosted-paired-gate.yml)
uses `workflow_dispatch` on the guarded research branch or `main`. It checks
out a source commit carrying this revision while retaining the frozen Rust
source tree and mathematical plan. It no longer consumes AWS secrets. Instead,
it builds the checked-in exporter in the pinned image, regenerates the complete
54-base public pilot in an isolated networkless container, and selects the
three exact K600 objects. [The independent selector](prepare_archives.py)
checks policy, seed, curve, columns, compressed BLAKE3, compressed SHA-256,
plain BLAKE3 and copy integrity against the frozen manifest. The preparation
receipt binds the producer commit and workflow-dispatch commit, and the runner
rejects a mismatched receipt before any worker stage. This is a byte-equivalent
input route to the three archived S3 objects, whose storage remains separate.
The job then prepares an offline Cargo registry, separately builds each
arm, verifies the imported signed-orbit domains, and runs one bounded natural
trial per arm. Each search runs inside a networkless 4 GiB, zero-swap Docker
cgroup. The independent receipt auditor checks candidate point membership and
the group equation if a candidate appears. An upload step retains setup logs,
stage stdout/stderr, Docker states, source and binary hashes, worker receipts,
and independent audit outcomes even when a stage fails.

The per-arm conservative budget is 600 seconds build, 360 import, 1,200 search
and 240 cleanup/audit, totaling 2,400 seconds; two arms use the previously
approved 4,800-second stage envelope. Fresh regeneration has its own 1,800
second cap and observed wall receipt; the revised 110-minute hosted timeout
allows both the preparation and stage envelopes. Any total-runtime model must
charge input construction or explain how reusable setup is amortized.
The worker's search timer is 1,100 seconds, and its outer Docker wall is 1,150
seconds. An incomplete job, cap, missing receipt or failed replay remains
censored. A first-trial partial or solver refutation alone cannot rank total
index-calculus runtime: useful relation yield, graph/rank growth, sparse linear
algebra, individual logarithm and matched rho are still required.
