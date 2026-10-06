# n=83 F6 wide algebraic block: exact admission, search still open

The [registered protocol](PROTOCOL.md) and [bounded-search amendment](AMENDMENT_1.md)
tested a two-word field table in the 128-variable Boolean wide backend.
The existing one-word table remains available for smaller fields and is
now explicitly `unsupported` above degree 64; `exhausted` stays set for
that incomplete result so existing callers cannot call it a proof of
absence. This is an algebraic block for the pinned K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, not a complete F6-IC method.

All eight [final release wide-backend tests](wide_tests_final.log.gz) passed.
The new table agreed with the legacy table at degrees 9, 17, and 31 and
with reference field products and squares at degree 83, including high
bits. A three-summand block from the 12-dimensional source subspace had
exactly 119 variables. For planted source points `[0,2,4]`, every
Boolean equation vanished under the exact 83-bit intermediate point,
and the three points replayed in the curve group. The source subspace
had 4,057 curve points; cofactor-four projection produced 4,054
distinct nonidentity subgroup points. The test independently verified
the four rational 4-torsion offsets and the identity `S=U+T` for exactly
one offset, where `R=4S` and `4U=R`. The algebraic source points must be
projected before a subgroup relation is counted.

The [native root probe](../../examples/f6_n83_wide_block_probe.rs)
used that planted source sum with `node_budget=1`. Its [raw output](root_probe.jsonl),
[exit status](root_probe.status), and empty [stderr](root_probe.stderr.txt)
give the following single exploratory construction diagnostic:

| Actual source points | Variables | Maximum matrix rows × columns | Setup | Root search | Peak RSS | Outcome |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 4,057 | 119 | 166 × 15,259 | 75.0 ms | 467.9 ms | 768,655,360 B | `exhausted` |

The solver counted two visited nodes and two reductions before the
one-node budget stopped the search. It returned no witness and made no
refutation claim. The probe finished before the registered 120-second
limit, so the **block construction and root reduction** clear that
admission gate. No conclusion follows about the cost or success of a
complete three-summand search, much less eight summands. The macOS host
rejected `ulimit -v 8388608`; the exact [limit attempt](limit_check.txt)
is preserved, and peak RSS was observed after the run rather than
enforced as a hard memory limit.

The final commands were `cargo test --offline --release --lib
cryptanalysis::wide_groebner::tests::`, `cargo build --offline --release
--example f6_n83_wide_block_probe`, and `gtimeout -k 5s 120s
target/bench_bins/f6_wide_block_final`, with the worktree's SSD-backed
`TMPDIR` and shared `CARGO_TARGET_DIR`. The [final build log](probe_build_final.log.gz)
and the [first compile failure](compile_failure_1.log.gz),
[test-fixture failure](test_failure_1.log.gz), and successful pre-status
logs are retained with lossless `gzip -n` compression. The compile
failure was a missing `Sync` bound in the new field-table trait; the
test failure came from using a Koblitz-curve constructor where only an
irreducible polynomial was needed. Both were fixed before the final
tests and probe. Physical arm64 macOS 26.6 used Rust 1.93.1 release
builds. No x86 hardware run is claimed.

Final source commit `31b8c9ba3246ba16989cb883aad012e677583fdd`;
wide backend SHA-256
`5b3e421bada26a451c645dceb9a0ac3ee38b9c30c4bdf7d3f22ba13f9a94603d`,
probe SHA-256
`fde4a969957de325ab0803adba429382165b119dd97b70ecc0c0d0763166659c`,
factor-base/curve source SHA-256
`1c98b492fafe80ae3a9dbea5400e7fceb114dec7b9a548e4d98f33a691f05b5d`,
and frozen probe binary SHA-256
`17678d80e0bb221c43b3008e28c73557e79b0b25fec682f3fac6067abb52f16e`.
The [protocol](PROTOCOL.md) and [amendment](AMENDMENT_1.md) hashes are
`e3676f1a2d6bebb39a75a7b541d145991e4e045333a4b0b188ec607f0a5438bc`
and `3c44b23e989ecf3fe99b15a4753828547ad20f2003d7d8e9f22a5e0eecfd3d66`.

This host has no CPU-isolation receipt, so the wall times are
exploratory. No ordinary relation, natural yield, recovered DLP,
complete F6/F4/F5 comparison, one-target IC online interval, or paired
rho run was measured. The full F6 and end-to-end speedups remain
unknown. The next algorithmic gate is a compositional high-arity
solver that can turn this exact source block and cofactor bridge into
one independently verified **ordinary** subgroup relation while
charging every failed attempt.
