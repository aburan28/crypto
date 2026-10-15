# N83 root4 high-K balanced-order records

This directory seals the four-arm K3550/K3596 root4 comparison in
[`PROTOCOL.md`](../n83_root4_highk_20261015/PROTOCOL.md). Candidate manifests
name the exact base and root index. Four measurement rows preserve both fresh
processes per candidate and the common 4,000-point public panel.

The executed source commit is `8a37349f3d3cc2bc2960b03167adcf57da83ba77`.
The public-point SHA-256 is
`48a243619e270abc096e172fc669bb03494fd8c09d4e39c55935121408d44eea`;
the known-scalar sidecar SHA-256 is
`fa660a434e2618a6fb3a9c251843c6bc4c9cc596c612c244ca2bb2ff8d54262d`.
The sidecar is used for independent replay outside the measured intervals.

The sealer checks the exact raw ZIP and independent analysis digests, source,
fixture, base, resource envelope, and analysis identities. Repeating the
command verifies byte equality of each canonical record. Artifact and replay
commands will be recorded after the run completes.
