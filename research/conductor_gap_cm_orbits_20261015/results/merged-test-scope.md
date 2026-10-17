# Combined-tree validation scope

Pinned Rust 1.98 checks on the merged tree:

- Strict Clippy for the library, ecbench binary, research examples and all test
  targets passed. The final helper Clippy also includes the posthoc rho module.
- Release `cm_orbit`: five passed.
- Release `public_fixture`: three passed.
- Release `fixture_rejects`: three passed; these filters overlap, so their
  counts are not summed as independent tests.
- The attempted `ecbench::fp4::tests` filter matched zero tests. That receipt
  is retained and is not counted as FP4 runtime validation. Main's FP4 code
  compiles in the combined tree; this workstream does not change its algorithms.

The first optimized rebuild was deliberately interrupted as shared memory
approached its limit. A subsequent test-check compile overlapped the large
constructor and caused a cgroup OOM kill of that constructor. Its receipt is
retained separately. After the constructor exited, the release rebuild and
the three relevant test filters completed successfully; see the `*-solo.log`,
`merged-fixture-tests.log` and `merged-fixture-rejection.log` files.

The earlier full suite had 4,116 passes and 117 ignored tests. That is a
historical full run, before the base merge and final fixture guards; it is
not represented as a new full-suite execution on the combined tree.
