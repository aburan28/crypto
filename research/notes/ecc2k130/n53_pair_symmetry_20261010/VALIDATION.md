# Pre-input validation

The pair-swap index preserves every n53 root-table key and value for the
23,320-point K220 base. Its focused release test passed both n13 and n53
cases with source SHA-256
`b7b290961b0137ce728b80b9ff92934013328c930dbf32fe0be00c01ae788471`:

```sh
cargo test --offline --release --example koblitz_orbit_dlp_fast_online pair_symmetry_tests -- --nocapture
```

`validation_logs.tar.gz` preserves the raw `validation_example.log`, which
records `2 passed; 0 failed`; `validation_logs_manifest.json` verifies its
SHA-256 and byte count. `rustfmt --check
--edition 2021 examples/koblitz_orbit_dlp_fast_online.rs`, Python syntax
checks for the freeze, generation, runner, and archive scripts, `git diff
--check`, and `python3 archive_controls.py` also passed.

The required broad command `cargo test --offline --release --lib` compiled
then exited 101 after the test
`cryptanalysis::isogeny_walk::million_store::tests::staged_parts_reconstruct_the_exact_verified_certificate`
overflowed its stack and the process aborted with SIGABRT. A focused rerun
of that test reproduced exit 101; its raw
`validation_lib_failure.log` is in the same verified archive. This branch
changes only the standalone n53
example and its research evidence; `src/cryptanalysis/isogeny_walk` has no
diff. The broad suite therefore has an explicit failing status rather than
a pass. The focused n53 example test and archived independent rank/target
replays establish the checks relevant to this index change.
