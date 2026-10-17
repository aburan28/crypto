# Structured-orbit import and first N83 search gate

The frozen 332,166-point factor base imported successfully. The one authorized
five-summand trial returned **`UNKNOWN_solver_cap`** after exactly 10,000 SAT
conflicts. It produced zero models and zero relations, so the retained matrix
has zero rows and rank zero. This is a bounded unknown: it gives neither a
natural relation-yield estimate at a larger search limit nor a complete cold
index-calculus runtime. There is no selected factor-base winner.

The user authorized 1,800 local seconds beyond the previously exhausted
7,200-second pilot for this exact gate. The offline build was charged 600
seconds and the import/search gate 1,200 seconds, regardless of their shorter
observed duration. The cumulative conservative charge is **9,000 of 9,000
seconds**. The local authorization file stays outside the checkout; its SHA256
is bound in the committed receipts and [evidence manifest](evidence-01/evidence-manifest.json).

## Frozen source, object and environment

| Item | Frozen value |
| --- | --- |
| Source commit | `5e798094406e11a137f6ce6c906bc95570256c51`, clean before build and search |
| Compressed object | `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v3-structured-orbit/a0/objects/aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389.jsonl.gz` |
| Downloaded-object SHA256 | `c5de0b26e84ad80bfb7be6a214b89c9682e356effb6079c7abf28cd1151ffb03` |
| Worker ELF SHA256 | `afe0a8df6d9d761bb6a6e111512ea57a02bebbc3839fc15c000538d619ed2956` |
| Linux arm64 image ID | `sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31` |
| Trial | Public `a=0, fixture=0`, `m=5`, `ChainedS3FrobeniusOrbit`, one trial, 10,000 conflicts, at most 32 models |
| Search envelope | One CPU, 4 GiB memory, zero swap, no network, 710-second worker and 720-second container wall caps |

The [protocol](PROTOCOL.md) defines the launch and failure classifications.
[Build](evidence-01/build/build.json), [gate configuration](evidence-01/gate/config.json),
[outer](evidence-01/gate/outer.json), [worker configuration](evidence-01/gate/run/config.json)
and [worker result](evidence-01/gate/run/result.json) preserve the exact source,
input and command bindings. The source tree changed only after all frozen-input
hashes were checked, to add this report and update the study summary.

## Observations and charge

| Phase | Receipt outcome | Observed wall or peak | Conservative charge |
| --- | --- | ---: | ---: |
| Offline release Linux build | `PASS_build_only`; exit 0; cleanup verified | 305.088 s | 600 s |
| Archive import in isolated container | `PASS`; 332,166 points, 2,001 signed orbits, 166,083 x coordinates | 6.027 s | Included below |
| One five-summand search container | Exit 0; `UNKNOWN_solver_cap`; cleanup verified | 569.745 s; 590,553,088-byte cgroup peak | Included below |
| Import/search supervisor and cleanup | `UNKNOWN_solver_cap`; no input change or OOM | 577.661 s | 1,200 s |
| Entire additional gate | Bounded unknown | 882.748 s observed build plus supervisor | 1,800 s |

The [import receipt](evidence-01/gate/verify.stdout) reports all archived rows
and generic orbit-map checks passing. The [search result](evidence-01/gate/run/result.json)
reports one SAT call, 10,000 conflicts, one solver unknown, zero models, zero
invalid models, zero verification failures, zero relations, zero independent
relations, zero dependent relations and rank 0 of 2,002 columns (2,001 orbit
columns plus the target column). Its one attempt is explicitly `Unknown` with
no decomposition indices. The worker reported 568,491.285 ms wall time and
recorded relation collection and rank activity in
[events.jsonl](evidence-01/gate/run/events.jsonl). These phase measurements
describe this censored gate and are not a complete index-calculus runtime.

An [independent receipt audit](evidence-01/independent-audit.json) rehashed
the committed raw evidence, all frozen input paths before report edits, the
build receipt and worker binary. It invoked the supervisor's import and search
receipt validators, checked matching statuses and the 9,000-second cumulative
charge, and confirmed that both experiment containers and the build container
were removed. All checks passed. The [budget receipt](evidence-01/gate/budget.json)
records zero additional seconds remaining. No further N83 search was launched
under this authorization.

The result leaves relation yield at a useful cap, independent matrix-rank
growth, sparse linear algebra, the individual-log phase and matched cold
comparison across bases unresolved. A rank-zero result from a single capped
trial cannot rank the structured base against the other stored candidates.

## Post-gate validation

The [validation manifest](validation-01/validation-manifest.json) hashes the
retained test logs. The study and boundary Python suites pass 29/29 and 39/39;
the touched worker and exporter examples pass 2/2 and 1/1. The required
`cargo test --release --lib` run was started with a larger Rust thread stack.
It completed 3,824 passing, 71 ignored and 12 failing tests before being
interrupted during unrelated long-running tests. Nine failures were in the
loopback collaboration group, which then passed 37/37 with local network
access. Three JV-cover failures in the broad concurrent run each passed when
rerun alone. The broad library suite therefore remains **unverified as a full
pass**; neither its failures nor the isolated reruns change the frozen N83
gate outcome.
