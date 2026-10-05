# Linux replay build-size failure and focused correction

On original PR #1374 head `f1267f8c9f0422094afe9ab28a506e49e169c75f`,
the Ubuntu `disclosed-sat-replay` check failed while running the file-only
watchdog control. The [original job](https://github.com/aburan28/crypto/actions/runs/37263655866/job/111615852069)
shows 22/23 worker tests passed; `journal::tests::watchdog_killed_writer_retains_start_and_original_drain_receipt`
panicked at `journal_tests.rs:211` when `measured_child_request` returned
`input is not a bounded regular file`. The helper was not launched. The
macOS counterpart and the other 15 checks passed on that head.

The reader rejects a symlink, a nonregular file or a file larger than its
128 MiB bound before launching the child. The test uses its own Cargo test
executable as the file-only helper. On Ubuntu, the previous workflow built it
with Cargo's default debug information; excess binary size is a likely cause,
but the failed job did not record its file metadata. That cause is an
inference, not a measured size from the failed runner.

The versioned replay workflow now compiles its development and test targets
with debug information disabled and incremental compilation off. This keeps
the existing 128 MiB reader guard, exact executable hashes, test assertions,
watchdog and process-group drain intact. The test prints the helper's actual
size and regular/symlink status; `--nocapture` retains the values in the next
Ubuntu/macOS replay logs. Source tests still use the exact same code and
fixtures. No scientific registration or closed confirmation set is rerun.

The [next Ubuntu replay job](https://github.com/aburan28/crypto/actions/runs/37353937752/job/111911442372)
passed both worker and controller watchdog controls. Its worker helper was
measured at 7,254,928 bytes, regular and nonsymlink; the worker suite passed
23/23 and the controller suite 24/24. The macOS replay job passed too. This
verifies the correction on that Linux runner. The original failed job did not
retain a helper size, so excess debug-binary size remains the likely cause,
not a proven measurement of that original file.

On the local physical macOS ARM64 host, both exact watchdog controls passed
with the revised Cargo profile environment (`CARGO_PROFILE_DEV_DEBUG=0`,
`CARGO_PROFILE_TEST_DEBUG=0`, `CARGO_INCREMENTAL=0`), using the shared native
busy lock and offline locked dependencies. The worker test measured a
4,593,712-byte regular, nonsymlink helper; the `icprog` copy measured
12,730,208 bytes with the same file properties. Each control exercised the
original two-second watchdog and process-group drain. This local result does
not establish the Linux cause or a controlled performance measurement.
