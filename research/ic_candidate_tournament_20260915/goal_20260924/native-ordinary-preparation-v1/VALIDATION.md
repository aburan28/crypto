# Implementation validation

This is a mathematical data-audit implementation, not a new preparation result.
The immutable historical 216-query F5 certificate is a retained control only.

* Rust formatting and `git diff --check`: passed for the implementation draft.
* Focused local Cargo controls: queued through the shared repository busy lock;
  another ecbench process owns that lock. No competing benchmark is interrupted.
  This is not a test pass and establishes no measurement isolation.
* Portable Linux/macOS controls and current-head CI: pending.
* Native target-free producers, their source-bound execution and per-attempt
  costs, natural-panel registration, and fresh paired comparison: not performed.

The first full worktree checkout failed because the volume ran out of space.
Git rolled that empty checkout back; the owned worktree was recreated sparsely
from the accepted dashboard merge. No scientific inputs or results were deleted
and no frozen registration was invoked. Local build commands reuse an owned
inactive debug Cargo target cache, separate from immutable scientific capsules.

Validation command (native, no solver dispatch):

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  env CARGO_TARGET_DIR=/private/tmp/ic-native-f5-build-cache-20261003 \
  cargo test --locked --offline --bin icprog ordinary_preparation::tests -- --nocapture
```

Replace the pending entries with exact terminal outcomes before accepting this
implementation. Do not report the retained fixture as new natural yield.
