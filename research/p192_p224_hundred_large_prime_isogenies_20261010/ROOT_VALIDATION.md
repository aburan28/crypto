# Required library validation and preserved failures

The initial required release library compilation failed with 643 errors. Three
already merged upstream type/API repairs were applied unchanged as commits
62a96b226, e76b73238 and 7e0bc7c11. The first runtime attempt then exposed these
independent integration issues:

- The trait reader rejected catalogue Montgomery/Edwards forms and extension
  coefficient arrays. It now preserves their exact characteristic and degree,
  checks the recorded order/trace identity, and explicitly marks coordinate/order
  certification `not_evaluated`. Existing binary and prime short-Weierstrass
  adapters retain their actual certification checks.
- A legacy-float regression used serde_json's feature-dependent direct float
  decoder as its reference. The reference now uses Rust's standard correctly
  rounded decimal parser. The production compatibility decoder and frozen
  records are preserved.
- Streaming stored-certificate hashing placed a one MiB buffer on the stack.
  Nested replay exceeded the default test-thread stack. The same buffer and byte
  sequence now live on the heap. The original abort and a larger-stack diagnostic
  that isolated the cause remain retained.
- Sparse checkout omitted the committed SOURCE42 and CryptoPro-B chain fixtures.
  Those exact tracked paths were restored; their contents were not replaced.
- The Hasse curve-generation regression compared generators from independent
  square-root implementations as identical representations. They can choose
  opposite roots. The corrected curve-level check compares p, a, b and the exact
  group order, requires the generated point to equal the frozen point or its
  negative, and verifies subgroup annihilation. Frozen benchmark generators and
  their EC1 identities stay unchanged. The generated representation is not
  claimed to have the frozen generator's EC1 identity.
- Local socket tests fail under the filesystem/network sandbox with
  `Operation not permitted`. The isolated socket control passed with localhost
  access. The final full suite runs with that same access.

The earlier default-stack/full-suite run was interrupted after known failures
had been corrected and while slow tests remained active. Its stdout/stderr and
SIGTERM outcome remain preserved. It is a failed/interrupted diagnostic, not a
completed test pass. The final complete suite is recorded separately.

Every scientific map acceptance check remains intact. These library and test
integration repairs do not alter the frozen producer/verifier executable used
for the 100-map experiment, source commit 49a280dc60972c9d16b2205b4652cc6c5d87284e.
