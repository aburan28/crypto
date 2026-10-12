# Native full-capsule custody implementation checks

Physical macOS ARM64 local checks under the shared native busy lock:

- All 75 `icprog` tests passed, including eleven F5 controller/custody controls.
  Archive fixtures contain placeholder bytes and execute no worker. They retain
  custody separately from execution admission, require the external seal, reject
  overwritten receipts, modified/resealed inputs, missing sidecars, truncated
  archives and prohibited execution claims.
- Both shared archive controls passed: file mutation, omission, duplicates,
  unsafe paths, links, invalid trailers, and unregistered or duplicate directories.
  The existing SAT example's original rejection control also passed after the
  parser moved into the shared native module.
- The rebuilt SAT custody CLI replayed all 5,963 files of its original published
  capsule. The archive SHA-256 remains
  `150e19a1f72a3c99f289aa5ad3a94bee504611e41adf82cc24b00b5bb3e2da37`;
  the original registration, producer and frozen checker digests are unchanged.
  `local-macos-arm64-shared-parser-sat-custody-v1.json` records the initial shared
  parser replay; `v2` records the final parser with directory-uniqueness checks.
  Neither extracts files, executes a binary or reopens the consumed registration.
- Formatting and diff checks passed. Local Clippy completed with no new warning;
  the older local compiler reports existing unrelated library and unknown-lint
  warnings. CI runs its separate pinned strict lint check.

No actual F5 capsule is published or scientific job dispatched by these controls.
The compact F5 build publication's ignored-lockfile failure remains in the parent
PR. Its corrected original data also passed replay from a clean Git export of
commit `caf12485c9779c9437e11fc75484d8fd71fdbabf`, excluding local untracked files.
Final-head Linux/macOS CI and review remain acceptance gates for this follow-on.
