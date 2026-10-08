# Exit-status regression receipt

The fix is based on parent PR #1338 head
`1b1949e31c6938bd3374a5358978e4a29e119623`. The reproducer is the
published n41/K20 pilot point from PR #1353. Before this patch, each of its
three retained producer rows had `targets_solved: 0`, `targets_failed: 1`,
target-row `exit_code: 0` and process status 0; each independent replay
rejected the row. Those historical artifacts remain unchanged in #1353.

On `arm64` macOS 26.6 with Homebrew `cargo 1.93.1` and `rustc 1.93.1`, the
following completed successfully:

```sh
cargo build --offline --locked --release --example koblitz_orbit_dlp_fast_online
experiments/koblitz-failed-target-exit-20261004/verify.sh target/release/examples
rustfmt --check examples/koblitz_orbit_dlp_fast_online.rs
sh -n experiments/koblitz-failed-target-exit-20261004/verify.sh
shellcheck experiments/koblitz-failed-target-exit-20261004/verify.sh
```

The end-to-end script exited 0 with `PASS: failed target exits 1 after
preserving evidence; solved target exits 0`. It asserted that the frozen
point at K20 now yields process status 1, target-row `exit_code: 1`, full
rank 20, zero solved and one failed target, while summary, target row,
base dump and rank trace still exist. The same point at K85 yields process
status 0, target-row `exit_code: 0`, full rank 85 and one solved target.
Both arms completed under the 30-second cap. This validates the reporting
contract, not performance.

`VALIDATION_SHA256SUMS` pins the patched source, regression script, exact
point and locally built release binary. It passed `shasum -a 256 -c` before
publication. The binary is reproducible from the committed source and lockfile
but is not checked into Git; the receipt records its local hash.
