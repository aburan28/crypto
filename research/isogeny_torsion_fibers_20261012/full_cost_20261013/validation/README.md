# Native run and verification commands

Host: Darwin 25.6.0 arm64, shared/L0, rustc 1.93.1, Typst 0.15.1.
The protocols were committed before the corresponding changed runs:
V1 `6d645eeb30714e03d62e3e5915c4a619f0ad0180`, V2 `345fe6b42f`,
V3 `ef5a5f3cd4`. All target points and source/destination identities are
emitted into the raw TSV and `identities.tsv`.

The run command for each frozen version was:

```sh
cargo build --release --example isogeny_fiber_full_cost_20261013
target/release/examples/isogeny_fiber_full_cost_20261013 > research/isogeny_torsion_fibers_20261012/full_cost_20261013/raw-vN.tsv 2> research/isogeny_torsion_fibers_20261012/full_cost_20261013/run-vN.log
```

The source snapshots `runner-v1.rs`, `runner-v2.rs`, `runner-v3.rs` are the
exact three run versions. `runner-v1.rs` has SHA-256
`8038f0e7cd857c4336c7924804243091432d170b44223b1ddb337af5cc6412db`,
which matches the recorded pre-V2 source hash. V2's source was
reconstructed by reversing V3's only arm-order patch; the reconstructed
binary matched all 768 deterministic V2 records. V3's source is byte-identical
to the committed example path and has SHA-256
`0edc04eea9109b7cd5400cc5c2c1ee958b03e23ee3b21e66f94aa340ce197220`.
The original V3 release binary SHA-256 is
`364366018ac1d34e680f1c3d12258bbe125f43b8e92bac6b4f4abbf4bbbf39d3`.
The original V1 release binary hash recorded before V2 was
`0efdc98b93e7c512db91ad9cb741eb16fee21a0c1f5f55309ced6f26e6fd5c5c`.
V2's original binary hash was not captured; its source reconstruction and
full deterministic replay are preserved. Original raw hashes are in
`../ARTIFACTS.sha256`.

The independent checker imports no curve or relation code from the runner.
It implements its own affine group law and modular inverse, reconstructs
the deterministic public targets, factor bases, pair tables, relation
trials, row rank/recovery, target descent and scalar checks. It also checks
phase partitions and A/A count identity. Commands:

```sh
rustc --edition=2021 -O research/isogeny_torsion_fibers_20261012/full_cost_20261013/verify.rs -o /private/tmp/isogeny-fiber-full-cost-verify-v3
/private/tmp/isogeny-fiber-full-cost-verify-v3 research/isogeny_torsion_fibers_20261012/full_cost_20261013/raw-v3.tsv
rustc --edition=2021 -O research/isogeny_torsion_fibers_20261012/full_cost_20261013/summarize.rs -o /private/tmp/isogeny-fiber-full-cost-summarize-v3
/private/tmp/isogeny-fiber-full-cost-summarize-v3 research/isogeny_torsion_fibers_20261012/full_cost_20261013/raw-v3.tsv
```

The checker output is `independent-v3.log`; V1/V2 were rechecked by this
stronger checker too. To reproduce the source snapshots without modifying
the tree's Cargo example, build each with `rustc --edition=2021 -O`,
`--extern crypto_lib=target/release/deps/libcrypto_lib-2f648d1af79d9c46.rlib`,
`--extern libc=target/release/deps/liblibc-34af9988a9ad958c.rlib` and
`-L dependency=target/release/deps`; run it to `replay-vN.tsv`. Build
`replay.rs` with `rustc --edition=2021 -O` and compare the original and
repeat. Each `replay-check-vN.log` reports 768 identical deterministic rows.
Clocks and process RSS are excluded from byte-identity comparison.

Visuals were generated from the two `.typ` graph sources with
`typst compile graph-*.typ graph-*.svg` and
`typst compile --ppi 150 graph-*.typ graph-*.png`; `report.typ` and the
three `run-vN.typ` files yielded PDFs. The same PDF/SVG/PNG bytes were
copied to `benchmarks/20261013-isogeny-fiber-full-cost/` and compared with
`cmp`. The native ICV1 and EC1 outputs came from `identity.rs` linked to
the repository's release `crypto_lib` implementation. The phase summaries
are produced by the checked Rust `summarize.rs`, not by Python.

The pinned Lean 4.34.1 module uses only `Init` and the standard kernel
foundations. In `lean/`, `lake build > build-final.log 2>&1` succeeded and
`#print axioms` reported no axioms for `ratio617`, `ratio2633`, or the four
`coprime*` theorems. `lean4checker --fresh` was unavailable on this host.
The first full eight-theorem attempt is preserved as
`NonisomorphicFour-v1.lean` and `build-v1.log`: it exhausted the default
recursion depth on the larger finite exclusions. The raised-limit source
`NonisomorphicFour-v2.lean` was interrupted after a bounded wait, with
`build-v2.log` (empty because Lean had not emitted buffered output). The
larger exclusions and elliptic-curve classification/map bridge remain
unformalized; neither interrupted attempt is a proof.

The first `cargo test --release --lib` pass is retained in `release_lib.log`:
3,995 passed, 14 failed, 113 ignored, exit 101 after 997.16 s. Nine
`pollard_collab` tests failed at local socket creation with OS error 1,
and the Z3 test could not open its generated solver input under the
sandbox. Four `jv_cover` checks exceeded their internal wall budgets under
heavy shared-host contention (including measured F4 times 75–86 s).
These are validation failures, not a clean suite result. A second pass,
`RUST_TEST_THREADS=2 cargo test --release --lib` with local access, reached
the same Nagao wall-sensitive failure under shared-host contention. It
was stopped after that failure (exit 130); its partial output is retained
as `release_lib_controlled.log` and is not a suite pass.

Targeted serial reruns with local access passed: the single Nagao check
(1/1 in `targeted_nagao.log`), the entire `jv_cover` module (14 passed,
one ignored in `targeted_jv_cover.log`), all 37 `pollard_collab` tests
(`targeted_pollard_collab.log`), and the exact Z3 round trip (1/1 in
`targeted_z3.log`). A complete serial `cargo test --release --lib` attempt
was stopped during its first costly ElGamal case after more than three
wall minutes and about ten seconds of test-process CPU under host load;
`release_lib_serial.log` and `release_lib_serial.status` preserve its
incomplete status. To cover all 4,122 registered library tests without
running the wall-sensitive module concurrently, the passing serial
`jv_cover` module was paired with a two-thread, local-access run of every
other library test (`release_lib_non_jv.log`). That run exited zero:
**3,995 passed, zero failed, 112 ignored, 15 filtered out** in 3,568.65 s.
The serial module passed **14, zero failed, one ignored** in 175.04 s.
The two partitions therefore cover every registered case with **4,009
passed, zero failed, 113 ignored**. This is complete split coverage, not
an unfiltered full-suite pass. `cargo test --release --example
isogeny_fiber_full_cost_20261013` compiled and exited zero; the example
defines no unit tests, so its harness reported zero run (`example_test.log`).
The native runner, independent checker, source-snapshot replay, Rust format
check, and Lean partial-module build have their separate receipts above.
No Python execution path was added or changed, so no legacy Python suite
is relevant to this Rust-only comparison.

From the repository root, `shasum -a 256 -c
research/isogeny_torsion_fibers_20261012/full_cost_20261013/ARTIFACTS.sha256`
verified all 75 listed study files after the final edits. The manifest
excludes itself and the ignored Lean build directory.
