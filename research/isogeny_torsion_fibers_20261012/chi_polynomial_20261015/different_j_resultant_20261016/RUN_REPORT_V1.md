# Run receipt: direct χ chord resultant exactness

Run ID: `different_j_resultant_20261016/exhaustive-v1`. Input:
`different_j_order_batch_20261015/PUBLIC_V1.json`, SHA-256
`f6c3126b4d95e3c457c9e552471bc2cdff322de0e0a9de9fbba701d7c05897ba`.
Source and checker versions are pinned in `RUN_SPEC_V1.json` and
`CHECKER_SPEC_V1.json`. Host: Darwin arm64, shared L0.

The requested run covered 229 subgroup points and two torsion points.
All 231 completed. Independent replay matched every returned pair set,
including 132 empty sets. The run counted 3,640 chord resultants,
7,410 polynomial gcds, 4,422 rational roots, 6,598 candidate additions,
and 6,358 rejected candidates. The metered replay used 7,094,272 bytes
peak RSS. Its real time was 0.39 seconds on the shared host.

The first meter command failed with `sysctl kern.clockrate: Operation not
permitted` in the sandbox; its exit, stderr and partial output are preserved
as `EXHAUSTIVE_V1_METERED.*`. The escalated repeat completed and passed
independent replay. The 512 MiB ceiling was satisfied. No comparative
performance ratio is computed from this exactness run.

Raw records: `EXHAUSTIVE_V1.jsonl`,
`EXHAUSTIVE_V1_METERED_ESCALATED.jsonl` and their corresponding exit,
stderr and `.check.json` files. The graph is
[`RESULTANT_GRAPH_V1.svg`](RESULTANT_GRAPH_V1.svg); the full result
analysis is [`RESULTS_V1.md`](RESULTS_V1.md).
