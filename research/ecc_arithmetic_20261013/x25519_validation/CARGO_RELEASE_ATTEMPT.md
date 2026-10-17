# Release library test attempt, 2026-10-13

The focused release test was started from the repository root with:

```sh
CARGO_BUILD_JOBS=1 CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/ecc-arith-x25519 python3 tools/isolated_bench.py busy -- cargo test --offline --locked --release --lib ecc::f25519::tests::arithmetic_matches_biguint_at_boundaries_and_seeded_inputs -- --nocapture
```

It waited behind the shared build lock, then began compiling dependencies and
the large library test binary. The last visible compiler output consisted of
unrelated unused-variable/dead-code warnings; no test summary was emitted.
After more than twenty minutes without a test result, `uptime` reported load
averages `256.06 228.91 189.67`. The interactive command was interrupted and
returned exit code 130 with `KeyboardInterrupt` from `isolated_bench.py` while
waiting for its Cargo subprocess. The terminal stream was not captured as a
complete file; no full compiler log or test result is claimed.

This is an infrastructure-limited attempt, not a failed arithmetic test. The
standalone source-bound field and ladder checkers in this directory completed
separately. The focused Cargo test and complete protocol regression remain
open on a machine with sufficient CPU headroom.
