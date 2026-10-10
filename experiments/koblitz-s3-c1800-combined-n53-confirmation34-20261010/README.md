# Fresh N53 targets for C1800 combined confirmation

The checked Sage launcher generated and independently verified **34 new N53
public subgroup targets**, T247–T280. Each has a distinct canonical one-target
workload ID and public-point file. The target-only [freeze](TARGET_FREEZE.json)
fixes **136 planned single-target process rows**: two C1200 packed and two
C1800 prefix-prefetch-projective runs per point, in alternating ABBA/BAAB
order. Its target and run order SHA-256 is
`c3183b0beff2b6bd815cc5b45c04ca8e0d523b8c48081a7ed0d12b98e3e17c87`.

The comparison uses the existing N53 candidate manifests and public-point
input law from the [eight-point combined study](../koblitz-s3-c1800-combined-20261010/README.md).
The C1200 baseline candidate is
`IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hf8b850f65ece`
with 127,200 actual usable base points; the C1800 combined candidate is
`IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h1ab23390e700`
with 190,800. Both use one target worker and the same relation seed. Their
source, manifests, base digests, and search schedules are bound in the freeze.
The primary comparison is each target's verified online interval after
reusable preparation and through scalar recovery check. The predeclared gate
requires both a geometric paired baseline/combined ratio of at least 2× and
a target-bootstrap 95% lower endpoint of at least 2×.

The target indices continue the parent study's fixed SHA-256 scalar law.
The [generator](make_targets_sage.py) replayed existing T200 and T246 fixtures
before writing T247–T280 and checked that their public points and workload
IDs were distinct from earlier N53 targets. The
[independent checked Sage receipt](local/verification.json) verifies every
scalar multiplication, subgroup membership, point file hash, and canonical
workload identity; all **34/34** passed. The
[generation record](GENERATED_TARGETS.json) contains local fixture scalars
for later verification. A remote input bundle should include only public
points and workload metadata. The [Sage runtime receipt](sage-runtime-info.json)
was saved before generation and reports a verified accepted build.

The 34-point count comes from the [repeat-aware variance estimate](../koblitz-s3-c1800-combined-20261010/robustness.json)
of the earlier eight-point panel: under its observed target and repeat
variance, 34 distinct points with two fresh runs per arm give approximately
80% power for a two-sided 95% lower bound above 2×. That is a planning
estimate; the actual confirmation host's fidelity and variance must be
checked by its pilot. The target-only freeze deliberately leaves process
timeout, memory cap, host-isolation receipt, and measured costs null until
that preflight is available. Every launched row, including a timeout or
failure, must be retained in the later execution receipt.

`GENERATED_TARGETS.json` has SHA-256
`6464e55889bd902eaa43508d4b069a5e65d53a409e365af6cb309a52d93ad225`;
the independent verification receipt has SHA-256
`66c3551212df0bf580d58eb79b64bbc513f4fcda384d99312ec83580e6e8d576`.
To regenerate fixture verification, use the checked repository launcher:

```sh
/Volumes/SSD990/cryptanalysis/sage -python /absolute/path/to/verify_targets_sage.py --output /absolute/path/to/new-verification.json
```
