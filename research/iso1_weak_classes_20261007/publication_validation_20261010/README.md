# Publication validation after merging current main

The two-branch class checks and all 34 focused native tests pass against the
merged source. The seven identity checks and desktop/mobile dashboard checks
also pass. All 5,509 files in the original experiment seal remain unchanged.
The census denominator is the ordinary trace classes with `t = 2 mod 4`;
the independent source populations retain their model-weighted input laws.

The publication merge combines research commit
`2ed59097cbb251c29cbc85f0d3ae3d331217850b` with main commit
`845fd323015e28d32b73fa25e551f830d5826e6c`. Its two conflicts concern generated
dashboard data and the curve-alias mirror. Both outputs were regenerated from
their authoritative sources. The original experiment directory and its PDF
were preserved byte for byte.

The owner-required Copilot request is recorded at
[the merge request](https://github.com/aburan28/crypto/pull/1588#issuecomment-6094165548).
After no response during local validation, the documented fallback is recorded
at [the resolution notice](https://github.com/aburan28/crypto/pull/1588#issuecomment-6094216972).
The base branch is merged into the PR head without rebase or force push.

The prescribed full-library command, `cargo test --release --lib`, was run
after materializing its tracked sources and eight compile-time fixtures. It
compiled and began 3,923 tests. Seven failures were reported before the process
aborted on stack overflow in
`isogeny_walk::million_store::tests::staged_parts_reconstruct_the_exact_verified_certificate`.
The complete stdout/stderr are retained as `root_final.stdout` and
`root_final.stderr`; this run is not a passing full-library check.

The five registry-trait failures repeat independently in the filtered replay.
They reject the existing Montgomery-form record
`icv1-fp221-tm3509210517603025598879416729141978-c324b837` as an unsupported
model. That record is present in both parent registries; the two-branch study
adds short-Weierstrass extension models. The other reported failures are the
staged j=0 test and the legacy-record deserializer test. The failure list is
preserved in the root output, rather than replacing it with a passing subset.

Earlier attempts that stopped on sparse source/fixture omissions are retained
separately. Materialization restored tracked files without changing their
contents. Root-library failures and current-head CI remain publication merge
gates; the mathematical and bounded-construction results are supported by their
separate native and PARI receipts in the sealed report directory.
