# Public input and single-target intervals

This follow-up to [PR 733](https://github.com/aburan28/crypto/pull/733) adds
supplied-point admission and native online intervals to the restored optimized
producers. Linux admission is pending. This is accounting and correctness work;
no candidate has been promoted and no improvement round has been consumed.

Measured jobs now require exactly one public point. Fixture construction runs
separately. IC starts its online clock after reusable logs are certified and
stops after final scalar replay. Rho includes its target-dependent setup, solve
and replay on the same point. Target queries, PDP, witness checks, descent and
recovery have exclusive native intervals. Cold process time and complete
instruction cost remain separate supplementary measurements. The pre-existing
bounded goal's cold improvement gates also remain in force.

## Local controls and retained failures

macOS arm64 controls are not the calibrated Linux performance panel. The native
source `e91457938cc3e6742ab7d1b3cc06203be0d2849c82851412cb2dc21bcd595844`
completed IC and rho on `13a0,23a1,37a0,43a1,61a1`: ten independent certificate
checks, all online/cold native sums closed. Their public points were generated
in separate processes. These five fixed vectors are neither new confirmation
targets nor estimates of natural relation yield.

The new concurrent worker tests exposed corruption in the archived arena's
relaxed publication of reused memory. The initial worker tests failed 9/20
repetitions, including an abort. A paired control compiled the same eight tests
and same library with only allocator ordering changed: relaxed ordering failed
18/20 runs (including a segmentation fault); acquire/release passed 20/20.
The normal Cargo build with acquire/release also passed 20/20 repeated runs.
These are bounded regression controls, not proof of arbitrary concurrent safety.
The final patch also removes a test's invalid second deallocation of a block
already released by `realloc`; all derived sources remain distinct from the
historical incumbent.

Raw reports, failed outputs, successful outputs, exact control executables,
source files, manifests and build commands are in the registered archive
`evidence/ic-public-native-controls-20260925.tar.zst`:

- SHA-256: `af5154ad4443a157c9a36ea552cb52def611c6b27669e3cc0740c83ecee78bd8`.
- 1,172,054 compressed bytes; 178 files; 4,948,881 uncompressed file bytes.
- Extract with `zstd -dc ARCHIVE | tar -xf - -C NEW_DIRECTORY`.
- `ic-public-native-controls/files.json` hashes every retained file except itself.
- `v2` and `v3` retain manifests and changed sources; unchanged dependencies are
  reconstructed from the sealed round-0023 source and PR 733's phase patch.
- `native-v3-worker` is the exact measured native executable. The two
  `worker-*-control` binaries are test executables, not benchmark candidates.

The [producer protocol](../../producer/PROTOCOL.md) requires 39 IC native/profile
pairs and 13 rho native/profile pairs across the three restored source policies.
`producer/audit.py ARTIFACT_DIRECTORY` reconstructs their canonical records,
phase ledgers and native intervals from transported artifacts. The old schema-2
auditor and archive remain unchanged.

## Remaining admission gates

Linux integration and fresh artifact replay must pass before these producers
are admitted. The development/promotion drivers still need canonical admission
migration; the current `icx` source, instrumentation overhead and a strong rho
reference still need qualification. Only then can a frozen comparison panel
launch. This checkpoint does not satisfy the three-round improvement objective.
