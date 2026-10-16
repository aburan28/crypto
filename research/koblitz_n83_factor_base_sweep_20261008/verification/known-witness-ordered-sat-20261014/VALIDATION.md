# Validation and host-contended rerun

The [bounded control audit](audit.json) rechecks with `audit.py --check`.
Independent re-extraction of the two retained S3 archives reproduced
`domains.json` at SHA-256
`63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb`.
The current example passes `rustfmt --check`; its release test harness compiled
and exited zero ([log](example-release.log.gz), zero example unit tests).
The focused two-residual library test added after the last documented complete
suite passes on the current source ([log](focused-library.log.gz), one test).
The N83 study passes 63 Python tests ([log](study-python.log.gz)); boundary
autolab passes 43 ([log](boundary-python.log.gz)).

The earlier [complete release-library result](../cross-stream-partial-worker-source-20261013/RESULT.md)
reports 4,001 passed and 113 ignored at parent commit
`bbbf03e94782aca600de5345ce65e7739c9cf4fd`. A diff from that commit to
the current branch changes no production-library code, `Cargo.toml`, or
`Cargo.lock`: the sole `src/` change is the added SAT-candidate assertion in
`cryptanalysis::koblitz_large_prime::tests::mixed_orbit_domains_lift_disjoint_partial_and_fail_closed_on_caps`.
That exact test passes in the focused current run and appears as `ok` in the
retained partial full rerun. The earlier complete-pass claim is documented in
the linked result, but its raw full-suite log was not retained there.

The current `cargo test --release --lib --offline -- --test-threads=1` rerun
did not complete. At host load averages of 228.04, 191.53, and 153.60, it had
completed 1,836 test lines in 44 wall minutes while receiving about eight
minutes of CPU. I sent SIGTERM to that test process only; Cargo recorded exit
101 and signal 15 in the [partial log](release-lib.log.gz), with an exact
[exit receipt](release-lib-exit.json). This is an interrupted infrastructure
attempt, not a passing current full suite or a failed code assertion. The
[attempt ledger](validation-attempts.json) also retains an earlier lost-handle
partial log and a blocked detached launch. A complete current full-suite
result is therefore unverified; the prior complete run plus current focused
test supply the validation coverage for the unchanged production library.

The [validation receipt](validation.json) binds every compressed log, the
partial exit record, and the exact source diff from the earlier complete pass.
Recheck it with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-ordered-sat-20261014/check_validation.py --check
```
