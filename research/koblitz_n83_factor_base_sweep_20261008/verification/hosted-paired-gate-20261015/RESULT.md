# Hosted Linux fallback for the N83 K600 Gray/hash first trial

Status: `SOURCE_PREPARED_NO_HOSTED_SOLVER_RUN`. The local Docker daemon is
unavailable, so this directory prepares a separately identified GitHub-hosted
Linux ARM comparison. The [original paired plan](../k600-gray-hash-paired-protocol-20261014/plan.json)
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

The [manual workflow](../../../../.github/workflows/n83-hosted-paired-gate.yml) reacts
only when label `n83-hosted-pair` is added to the existing same-repository PR
on branch `codex/n83-mixed-orbit-supervisor-20261012`. It checks out the exact
PR head SHA. Repository AWS secrets are used only inside that guarded job to
read the three private S3 archives. The job checks their compressed SHA-256
values and the prior source and orbit audit before any worker stage. It pulls
the pinned image, prepares an offline Cargo registry, separately builds each
arm, verifies the imported signed-orbit domains, and runs one bounded natural
trial per arm. Each search runs inside a networkless 4 GiB, zero-swap Docker
cgroup. The independent receipt auditor checks candidate point membership and
the group equation if a candidate appears. An upload step retains setup logs,
stage stdout/stderr, Docker states, source and binary hashes, worker receipts,
and independent audit outcomes even when a stage fails.

The per-arm conservative budget is 600 seconds build, 360 import, 1,200 search
and 240 cleanup/audit, totaling 2,400 seconds; two arms use the previously
approved 4,800-second stage envelope. Setup and S3-transfer wall are recorded
separately by the hosted job and must be included in any total-runtime model.
The worker's search timer is 1,100 seconds, and its outer Docker wall is 1,150
seconds. An incomplete job, cap, missing receipt or failed replay remains
censored. A first-trial partial or solver refutation alone cannot rank total
index-calculus runtime: useful relation yield, graph/rank growth, sparse linear
algebra, individual logarithm and matched rho are still required.
