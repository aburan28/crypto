# Implementation validation

This is a mathematical data-audit implementation, not a new preparation result.
The immutable historical 216-query F5 certificate is a retained control only.

* Rust formatting and `git diff --check`: passed for the implementation draft.
* Original local Cargo controls: the old handle and temporary output no longer
  exist after the interruption. Their terminal outcome is unverified and is
  not counted as a pass. No scientific run is restarted from that observation.
* Portable controls passed on Linux and macOS at source head
  `0cf07c63b5bd618f7057afb88d8e61bd6c0f53d5`: seven ordinary-preparation tests,
  zero failures on each host. The full PR passed all 49 applicable checks.
  The raw CI logs are retained losslessly as `linux-ci-20261004.txt.gz` and
  `macos-ci-20261004.txt.gz`; decompress with `gzip -cd`.
* Integration with the current-main F5 result dependency is now present.
  Its new exact-head CI remains pending; the earlier head's pass is not
  acceptance of this integration.
* Native target-free producers, their source-bound execution and per-attempt
  costs, natural-panel registration, and fresh paired comparison: not performed.

The first full worktree checkout failed because the volume ran out of space.
Git rolled that empty checkout back; the owned worktree was recreated sparsely
from the accepted dashboard merge. No scientific inputs or results were deleted
and no frozen registration was invoked. Local build commands reuse an owned
inactive debug Cargo target cache, separate from immutable scientific capsules.
That temporary cache is now absent; the historical command below is retained
as the original attempted local validation, not a verified terminal receipt.

Validation command (native, no solver dispatch):

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  env CARGO_TARGET_DIR=/private/tmp/ic-native-f5-build-cache-20261003 \
  cargo test --locked --offline --bin icprog ordinary_preparation::tests -- --nocapture
```

The [Linux job](https://github.com/aburan28/crypto/actions/runs/37163794798/job/111322478302)
and [macOS job](https://github.com/aburan28/crypto/actions/runs/37163794798/job/111322478259)
both list all seven controls and their terminal pass. These are portable
correctness receipts, not hardware speedups or natural-yield measurements.
