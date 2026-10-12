# Bounded ordered first-pair N83 comparison

The source-matched five-summand ordered and unrestricted Frobenius-orbit arms
both reached `UNKNOWN_solver_cap` after exactly 1,000 SAT conflicts on the same
public target and K=2,048 base. Each produced zero SAT models, zero relations
and rank zero. These are censored one-trial observations; neither arm completed
relation collection, column logs, an individual log or a total cold
index-calculus run. No minimum-runtime factor base is selected.

The selected base contains 339,968 points and 2,048 signed-Frobenius orbit
columns. Its generic multi-limb replay checked every point and representative,
and its S3 upload/download receipt reports a matching compressed hash. The
stored object is
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851.jsonl.gz`.
An S3 HEAD check on 2026-10-11 confirmed its 10,313,647-byte object and
compressed BLAKE3 metadata; the [HEAD receipt](evidence-01/panel/s3-head.json)
is included with the original upload/download receipt.

The clean source commit was `b7f55cf1528f13419ba079711d3c817d764a7dae`.
The offline Linux arm64 build passed with worker SHA-256
`b84cc4a6036ea8041362bf29fcf2388a9c0f9c73d8c4f92e376968d97d893b3b`.
Both arms used that binary, the same object SHA-256
`aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`,
the same public target and fixture hashes, one CPU, 4 GiB memory, zero swap,
no network, at most 32 models and a 280-second external wall cap. The only
SAT-domain difference was the optional `x0 <= x1` constraint on the first S3
summand pair. The ordered arm uses 82 more source-count SAT variables and
emits up to 247 extra CNF clauses; installed counts depend on root
simplification and were not recorded by this search worker.

| Stage | Receipt | Descriptive process wall | Cgroup peak | Conservative charge |
| --- | --- | ---: | ---: | ---: |
| Linux worker build | `PASS_build_only` | 147.139 s | Not recorded (8 GiB cap) | 280 s |
| Ordered orbit search | `UNKNOWN_solver_cap` | 80.704 s | 518,672,384 bytes | 280 s |
| Unrestricted orbit search | `UNKNOWN_solver_cap` | 53.394 s | 506,798,080 bytes | 280 s |
| Cleanup and receipt reserve | Audited | Included in stage walls | — | 60 s |

The prior conservative charge was 9,900 seconds. This protocol charges the
full 900-second extension, so the cumulative charge is **10,800 of 10,800
seconds**. Both search containers and the build container were independently
confirmed absent after completion. The [audit](evidence-01/audit.json)
rechecks frozen hashes, exact receipt contracts, source/binary/object identity,
the same query point, original-curve query replay, S3 and generic replay
receipts, rank counts and cleanup. The [raw evidence](evidence-01/) retains
the build log and both supervisors' configs, source attestation, events, worker
summaries and outer receipts. Frozen source bytes are rehashed from the build
commit during audit. The worker ELF and 9.9 MiB compressed base stay
outside Git; the latter is stored at the S3 URI above.

The 27.310-second difference between these two outer wall observations is
descriptive only. There was one run per arm on a shared host, both searches
stopped at the same conflict cap, and neither yielded an equation. It is not
evidence of a complete-runtime speedup or of one factor base outperforming
another. The next meaningful gate requires a separate compute authorization:
measure natural relation yield and rank novelty at a useful solver limit,
then column-log and individual-log completion with a matched rho control.

Validation: the exhaustive first-pair comparator release tests passed 2/2;
the degree-7 ordered/unrestricted lifting test passed 1/1; the exporter
release tests passed 19/19; and the study Python suite passed 36/36. The
required broad `cargo test --release --lib` run reached its 900-second wall
cap without a complete summary. Its eight-thread run reported three failures
in `jv_cover`, outside the touched N83 paths; all three passed in the serial
`jv_cover` partition, which finished with 12 passed and one ignored. The
broad library suite remains unverified as a complete pass. The
[validation receipt](evidence-01/validation/validation.json) and
[serial log](evidence-01/validation/serial-jv-cover.log) retain the exact
commands and focused result.
