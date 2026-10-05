# Offline build and full byte custody

The second validation-only freeze passed from committed source
`34304b60d85f497698cfafbc9a385200868291e1`. Both binaries were built offline
from the complete copied source and vendored dependencies. The worker's actual
compiled identity, all original build receipts and independently decoded native
asset files passed the freezer's checks. Only the worker's `build-identity`
command ran: no preparation, exporter or CMS search was invoked.

| Identity | Exact SHA-256 |
| --- | --- |
| Validation registration | `02e12b456947c8034f096b47c16dd7b73eb6c15cedb6c0b79164f887ec1406f3` |
| Source/dependency manifest | `7d5746daada759205731467b02b8cff2f7830f26ee01eed1726d07d7a42c7922` |
| Compiled build descriptor | `5723d81ab6b7f3c73f2d6e3e39f839b1ddf5d5e8682f399f66f9904fd1bf87d3` |
| Worker | `641e8b783a16498a693e323b7b64929b7bbc12ef17882df42f3d7633d5aa1223` |
| Frozen controller/auditor | `8c8acafc5e15918cff442e1cf36b2fa07a879bef813e6175d6442a741e8be33a` |
| Full published capsule archive | `78d41838ed39d1b3bc0af16563de6e0f754b5fdb00fdbed23f3a97dc6a2de345` |

[`build-validation-v2`](build-validation-v2/PUBLICATION.json) retains the four
original input/registration sidecars, every original build log and receipt, and
the complete immutable tree in `capsule.tar.gz` (49,916,133 compressed bytes).
The native data checker verifies all 6,950 archived regular files against the
sealed inventory, including source, dependencies, binaries and native assets.
It rejects links, unsafe names, omitted/duplicate members, foreign files and
changed bytes. System `shasum -a 256` independently matched the archive digest.
Build target/cache outputs are excluded. The capsule is validation-only and
unconsumed; archive custody does not authorize scientific execution.

The publishing/data-replay tool was added after the frozen source commit. Its
source is versioned in this PR; it reads the original frozen snapshot as data
and does not replace the frozen binaries or their original identity receipt.
`ordinary-control-publish-build` accepts only an unconsumed validation-only
capsule with the external seal. `ordinary-control-replay-build` requires that
same external seal and starts no archived binary or native solver. It checks
the sealed source manifest, executable pins, complete original receipt set,
actual environment pair encoding, build flags, native asset pins and every
archive member. Original absolute build paths remain receipt data; no original
installation or physical Mac is needed for data replay.

The final local replay passed with explicit source-runtime admission false.
Four build-custody controls and all 11 controller/shared-contract controls
passed. New controls reject changed or omitted original receipts even after
regenerating the outer publication inventory, wrong external seals and runtime
claims. Clippy completed with the inherited warnings described in
[VALIDATION.md](VALIDATION.md); changed-file formatting and diff checks passed.
Linux/macOS CI now runs the same data replay and retains its own receipt.
Cross-platform CI acceptance remains pending at publication time.

The [first offline freeze failure](failed-freeze-v1/FAILURE.json) remains:
vendoring exited 101 because the local cache lacked locked dependency
`zerocopy-derive 0.8.59`. The exact original log, receipt and native-asset
receipt are retained. `cargo fetch --locked` downloaded that missing dependency;
the lockfile remained `f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79`.
The second attempt used a new directory and remained fully offline. Neither
failed nor successful build launched a scientific worker.

The [first archive-publication attempt](failed-build-publication-v1/FAILURE.json)
ended with observed exit 101, no output, an incomplete gzip stream and no
completion receipt. Its cause remains unknown. Its exact 35,263,607-byte
partial archive and sidecars are retained. A host process audit confirmed no
live publication child before the separately named second attempt. The second
attempt added stage logging and completed; this does not establish why the
first attempt stopped. No failed artifact is counted as verified custody.

Exact freeze command (run through the existing busy lock from this worktree):

```sh
/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004/debug/icprog \
  ordinary-control-freeze \
  --root /Volumes/SSD990/crypto/worktrees/ic-native-ordinary-controller-20261004 \
  --out /private/tmp/ic-native-ordinary-build-validation-20261004-v2 \
  --config research/ic_candidate_tournament_20260915/goal_20260924/native-ordinary-controller-v1/build-validation-config.json \
  --host-context research/ic_candidate_tournament_20260915/goal_20260924/native-ordinary-controller-v1/build-validation-host-context.json \
  --cargo /opt/homebrew/Cellar/rust/1.93.1_1/bin/cargo \
  --rustc /opt/homebrew/Cellar/rust/1.93.1_1/bin/rustc \
  --validation-only
```

Portable replay, with a new output path:

```sh
/tmp/ic-native-busy busy -- target/debug/icprog ordinary-control-replay-build \
  --publication research/ic_candidate_tournament_20260915/goal_20260924/native-ordinary-controller-v1/build-validation-v2 \
  --validation-registration-sha256 02e12b456947c8034f096b47c16dd7b73eb6c15cedb6c0b79164f887ec1406f3 \
  --out /tmp/ordinary-build-custody-NEW.json
```

This is physical Apple M4 Pro/macOS ARM64 build correctness and data custody.
The host is uncalibrated; the busy lock is not host isolation. No preparation
yield, one-target solve, comparative cost, speedup, scientific registration or
complete-goal acceptance follows from it. Actual separately frozen F5/CMS
scientific registrations and panels, runtime audits, fresh target admission,
independently verified online intervals and paired incumbent/strong-rho
comparison remain pending. Old control and confirmation sets stay closed.
