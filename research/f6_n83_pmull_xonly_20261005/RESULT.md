# n=83 F6 packed PMULL residuals: retained as an exact alternate method

The [preregistered protocol](PROTOCOL.md) compared the retained x-only
four-summand query with an ARM64 PMULL implementation of its degree-83
field arithmetic. The candidate packs field elements into `u128`, uses
four carry-less limb products per multiplication and two per square,
and reduces modulo `z^83 + z^45 + z^2 + z + 1` with fixed sparse folds.
One reference inversion per batch, all exceptional group operations,
the signed x-key lookup, sign check, and four-point group replay remain.
The new `solve4_xonly_pmull83` method uses the packed path only after
checking ARM64 AES/PMULL support and the exact field polynomial; it
otherwise calls `solve4_xonly`. The index build path is unchanged.

The [field test log](field_tests.log.gz) passed 10,000 deterministic
pseudorandom multiplication and squaring comparisons plus edge values.
The [geometry log](geometry_tests.log.gz) passed all eight module tests,
including exceptional residual x keys and exhaustive small-base witness
replay. All six native benchmark processes exited 0 with the same exact
`no_witness` T001 outcome. The [status table](status.tsv), [runner](run.sh),
six JSONL outputs, and six empty stderr files preserve each run. The
baseline and candidate [build logs](baseline_build.log.gz) and
[candidate log](candidate_build.log.gz) record offline release builds;
logs are losslessly compressed with `gzip -n`.

| Actual base points | Metric | X-only baseline | Packed candidate | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 33.392 ms | 23.823 ms | 1.402× |
| 258 | Exact query, three-run median | 14.016 ms | 2.294 ms | 6.110× |
| 1,048 | Index build, three-run median | 1,234.932 ms | 404.923 ms | 3.050× |
| 1,048 | Exact query, three-run median | 401.228 ms | 58.416 ms | 6.868× |
| 4,054 | Index build, two-run median | 6.641 s | 7.284 s | 0.912× |
| 4,054 | Exact query, two-run median | 3.644 s | 1.216 s | 2.998× |
| 4,054 | Maximum process RSS | 1,276,051,456 B | 1,276,297,216 B | ~1.000× |

The full runs were interleaved baseline, candidate, candidate, baseline.
Their exact query times were 3.485, 1.303, 1.128, and 3.803 s;
adjacent paired baseline/candidate ratios were 2.674× and 3.372×.
These two pairs describe observed spread, not a statistical confidence
interval. Full build times were 7.198, 7.084, 7.484, and 6.085 s in
the same order. Every full process used 4,054 distinct usable points,
8,219,485 unordered pairs, 4,108,723 signed-sum representatives, and
2,027 signed columns. The candidate full-query median was **66.6%
lower**, clearing the registered 20% exploratory gate. Both small-base
query medians improved. Full-base build median rose **9.7%**, narrowly
within its registered 10% limit; maximum RSS rose by 245,760 bytes,
within its 10% limit. The retained method is an alternate API path,
pending controlled timing and complete solver integration.

A separate [planted control](planted_control.jsonl) on the full K0 base
constructed a target from points `[0,2,4,6]`. Both the reference and
packed methods returned those indices, and both four-point sums replayed
to the target. The corrected fixture [passed once](planted_control_pass_2.jsonl)
and then [passed again](planted_control.jsonl) with the formatted probe
source; its [final build log](planted_build_final.log.gz) is preserved.
This is a correctness control, **not**
an ordinary-yield observation. The first planted fixture used `[0,1,2,3]`, whose sum was
the identity; its nonidentity fixture assertion exited 101 before either
solver ran. The [planted status](planted_status.tsv), empty first JSONL,
and [failure stderr](planted_control_failed_1.stderr.txt) preserve that
fixture error. The corrected probe and both build logs are preserved.

The baseline branch head was
`f30cc085ed0c0a947b9c501c438f644634034239`, with pair-index source
SHA-256 `3f2ef8bd0f2bb07dde0bdb7a9ad7f3da76e0ff1924a9ec04dbd200102ae4125f`.
The candidate source SHA-256 was
`9ce12438b77d0a77fbaa414d15a2d5f59e8555a1e3dd383368eb6e22454a887c`.
The candidate small, full, and corrected planted probe source hashes
were `bc335edb405ae87d3c3c40e0dd078a09efb3ed843f1687226f4fe8d26363a9cc`,
`cf6cc7e63fd26283b44f2edf11f7b7a3b6ba00bfc2460a0030189797c50b7916`,
and `935a1deac69ac7fc1affc1c13bd40c86b6fd4bb30eaa1f6d9f02ba514655aeb2`.
The protocol and runner hashes were
`b91bd76caf025498c9e6ae27866d1d097f9f72fcd041a1d6c361fef8fc384816`
and `83dab264abf82a7e961c7aaedfb40d66fed2a1a71ae8d13d0d8f87e8ba8b3515`.
The baseline small and full probe binaries were SHA-256
`7e7c14747539c4437a6274e7cbc8e3e38bb4bbbe99b0721c094adee1b2d738ca`
and `09c1d01c2ed9433b587f0986ebd00fdd818bdf8c2b82feb7c03e42e1399fcce1`.
The candidate small, full, and corrected planted binaries were SHA-256
`83b8f59d66422cfcbbe47dd840c6b476fac64cec6baea74d0c7debb6805a521c`,
`8923553ad79d2647383fd8a929b12fd4ddbfe0f880d95b9d641de295fc6e541b`,
and `9cc0c6bda7c4e6dd2f97852f37a6fc36b7276e48fbd410bb07d9812cb07809d4`.
Both arms used Rust 1.93.1 release builds on physical arm64 macOS.

The build and test commands ran from the repository worktree with
`TMPDIR` on the SSD and a shared `CARGO_TARGET_DIR`, using `--offline`:

```sh
cargo test --offline --release --lib n83_pmull -- --nocapture
cargo test --offline --release --lib cryptanalysis::f6_wide_geometry::tests::
cargo build --offline --release --example f6_wide_n83_pmull_index_probe --example f6_wide_n83_pmull_full_probe
cargo build --offline --release --example f6_wide_n83_xonly_index_probe --example f6_wide_n83_xonly_full_probe
sh research/f6_n83_pmull_xonly_20261005/run.sh
cargo build --offline --release --example f6_wide_n83_pmull_planted_probe
```

The local `cargo clippy --offline --lib -- -D warnings` attempt failed
on pre-existing unrelated macOS warnings and two `ecbench` Boolean
expressions; its [failure log](clippy.log.gz) contains no finding in the
F6 geometry module. Applicable CI checks remain the merge gate.

This shared host has no isolation receipt, so the timings and ratios
are exploratory component diagnostics, not controlled CPU speedups.
The full-base build gate is particularly close to its limit, and two
pairs cannot bound host noise. The exact four-sum query missed T001 in
every run; this base's uniform-target coverage ceiling is
`4.662e-12`. No natural relation yield, recovered logarithm, complete
F6/F4/F5 comparison, one-target IC online interval, or paired rho run
was measured here. The end-to-end IC effect remains unknown.
