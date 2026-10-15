# N83 root Bloom profile records

This directory seals the six-arm comparison in
[`PROTOCOL.md`](../n83_bloom_profile_20261015/PROTOCOL.md). Three candidate
manifests distinguish the exact Bloom size and hash count on the same K3300
factor base and root11 index. Six measurement rows retain two fresh processes
per profile and the common 4,000-point public panel.

The executed source commit is `42fc71945cd7087f78464abc65c4bae7c52a93e8`.
The public-point SHA-256 is
`35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491`;
the known-scalar sidecar SHA-256 is
`6e1181ba958708803316903a645463425bc79130074041d39c8094778e2f4379`.
The sidecar is used for independent replay outside the measured intervals.

The sealer checks the exact raw ZIP and independent analysis digests, source,
fixture, base, Bloom profile, resource envelope, and analysis identities.
Repeating the command verifies byte equality of each canonical record. The
artifact and replay commands will be recorded after the run completes.
