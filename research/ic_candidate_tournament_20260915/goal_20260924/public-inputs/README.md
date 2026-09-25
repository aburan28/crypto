# Public input and single-target intervals

Merged [PR 748](https://github.com/aburan28/crypto/pull/748), following
[PR 733](https://github.com/aburan28/crypto/pull/733), adds supplied-point admission
and native online intervals to the restored optimized producers. The final Linux
admission and independent transport replay passed. This is accounting and correctness work;
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

Subsequent source review found that this first rho control included reusable
field/Frobenius preparation inside its online clock. The archive is preserved,
but the corrected interval validator explicitly rejects those rho intervals.
`rho-preparation.patch` adds a hook after reusable preparation and before the
target is read, reports the actual arithmetic kernel, and preserves the old
entry for walk/counter equivalence checks. Linux validation must use that final
boundary; an earlier green run is insufficient.

The corrected local source
`9c3d4a51ec4e52bfd4f9dbb83675cadb95889aea7370d6ffc91abe380ed48ba0`
passes the same ten IC/rho public-point controls, eight worker release tests,
and 33 rho release tests (three pre-existing large tests ignored), including
the new prepared-entry counter/walk equivalence test. Its complete control
reports, executable, manifests and build logs are in the separately registered
`evidence/ic-rho-prepared-controls-20260925.tar.zst`. Its `files.json` checks all
retained file contents; use the archive manifest for the archive hash and size.
These local results remain correctness controls only.

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

The final implementation `159ec5ff086702a7fd5603d32eeacb48b6ef7361` merged as
`4a898fcd3464b71439bd9451bd466a6ec217ddc9`; all 20 reviewed files match the merge.
Every applicable check passed at that head. The
[final Linux run](https://github.com/aburan28/crypto/actions/runs/36113511435)
passed the following fixed controls, all reconstructed from downloaded artifacts:

| Source policy | IC native/profile pairs | Rho native/profile pairs |
|---|---:|---:|
| `both` | 9/9 | 3/3 |
| `scaled` | 15/15 | 5/5 |
| `pairinv` | 15/15 | 5/5 |

[History and exact source/binary identities](history.json) and [replay receipts](audits/)
bind each row. The separately registered
`evidence/ic-public-linux-controls-20260925.tar.zst` preserves both the initial
superseded run and the final run, including executables, all source files
(including `.cargo/config.toml`), frozen evaluators, protocols, process reports,
profiles and certificates. Extract it into a fresh directory using the command
above. Run each extracted artifact's `ic-producer-evidence/evaluator/producer/audit.py`
with its artifact directory as the argument. The initial evaluator replays its
historical boundary; that pass does not rehabilitate the superseded rho interval.
The archive manifest supplies its digest/size and `archive-validation.json`
records the fresh extraction and replay.

The development/promotion drivers still need canonical admission
migration; the current `icx` source, instrumentation overhead and a strong rho
reference still need qualification. Only then can a frozen comparison panel
launch. This checkpoint does not satisfy the three-round improvement objective.
