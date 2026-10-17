# Validation and execution record

The three full graph builds and copied-database replays exited 0. Native
exports checked every conductor population against the frozen plans. No
frontier remains. Six unit and two process integration tests passed.
All delivered file hashes and ZIP CRCs passed. Drive metadata readback
verified both uploaded files. Remote whole-file hashing was not performed.

Curve catalogue cover generation and recheck both report 2,126 verified,
zero invalid and zero unsupported models. Browser and leaderboard generation
and final freshness rechecks passed. Earlier freshness attempts were made
before generation completed and reported stale files; these outputs are
retained alongside the successful rechecks. No performance figures changed.

The first metadata import exited 101 because its stdout destination was an
empty JSON file inside the scanned directory; a .txt log avoided self-scan.
The failure and empty output remain preserved. The corrected import exited 0.

The full native debug traits build asserted in binary-subfield trace code.
It was stopped with SIGTERM (143), and its queued full check was not run.
The exact failing model and clean-base reproduction remain unknown. See
[the problem record](../../problems/curve-traits-debug-subfield-trace.md).
No assertion or correctness gate was disabled.

The native incremental refresh retained all 2,296 baseline records and
computed only six new prime records. It checked all non-registry input
hashes unchanged, updated the registry digest, and preserved the normal
serialization contract. A replay from the exact baseline in a scratch
working directory produced a byte-identical updated traits file (cmp exit 0).
This does not establish full-catalogue replay. Baseline SHA-256:
`7395b479aea3e90f3fd2010dd8d4f570c60dc4fb2baf83bda97165b0e659b844`.
Updated SHA-256:
`763e70f2a2e2863c06d05ba09023e32ed776b658a57fd0baa3dff5311b3d33ef`.

The one-off native refresh source is `refresh_traits.rs`. It links the
current checkout's compiled library, using this exact build command:

```sh
rustc --edition=2021 research/full_volcano_classes_20261017/refresh_traits.rs \
  --extern crypto_lib=/tmp/volcano-debug-target/debug/deps/libcrypto_lib-038cc0b52d457b73.rlib \
  --extern serde_json=/tmp/volcano-debug-target/debug/deps/libserde_json-103d5c7fd4432996.rlib \
  -L dependency=/tmp/volcano-debug-target/debug/deps -o /tmp/volcano-refresh-traits
```

Build the library with `cargo build --bin curve_traits`; dependency artifact
suffixes vary by build. Run the refresh only in a working directory containing
the pinned baseline traits file, new registry, and unchanged files named by
the baseline's `inputs` array. It deliberately refuses another baseline.
The initial output and guarded scratch replay are retained in
`evidence/traits-incremental*.json` and associated logs.

Global naming lint retains five historical names in two unchanged frozen
reports; the added-line check is recorded separately. Conductor was attempted
and unavailable (127). No new solver timings, operation counts, relation-yield
statistics or cross-method rankings were produced. Original A/C/D/E runs and
the full ECDLP/IC campaign remain unresolved.
