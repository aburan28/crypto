# N83 balanced-order factor-base records

This directory seals the four-arm comparison in
[`PROTOCOL.md`](../n83_k_abba_20261015/PROTOCOL.md). Each of the two
candidate manifests names an exact factor base and kernel. Four measurement
rows retain both fresh processes per candidate. The workload record fixes
the common 4,000-point public panel and the experiment's process order.

The executed source commit is `d5e6fc717bac432efba63564ba12081d464eebb6`.
Its public-point SHA-256 is
`35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491`;
the known-scalar sidecar SHA-256 is
`6e1181ba958708803316903a645463425bc79130074041d39c8094778e2f4379`.
The sidecar enters independent replay after production, outside measured
target intervals.

The sealer requires the exact raw ZIP and independently generated analysis
digests. It checks source, fixture, base, resource, and analysis identities
before writing canonical candidate, workload, measurement, selection, and
index records. Repeating the command verifies byte equality of existing
records. The source-pinned CI workflow and result report name the artifact
and exact replay command after the experiment completes.
