# Validation record

## Snapshot 03: security goals, pairing groups, and SIDH, 12 October 2026

- The source-pinned roster imports Ethereum's exact `alt_bn128`/`BN_SNARK1`
  variant and BLS12-381. `python3 validate_pairing_roster.py` passed for both:
  `ord_r(p)=12`, the published G1 cofactor relation, G1 and G2 curve
  equations, and `[r]P=O` for both generators. This independently caught and
  corrected one copied BLS12-381 G2 coordinate before export. The source
  records link [EIP-196](https://eips.ethereum.org/EIPS/eip-196),
  [EIP-197](https://eips.ethereum.org/EIPS/eip-197), and the
  [CFRG parameter draft](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/).
- `cargo test --release --locked --offline`: four passed. The Rust exporter
  `build` and `check` passed. The graph builder `--check` passed with 922 atlas
  nodes, 1,092 total nodes, 676 isogeny classes, 73 multi-model classes,
  and the same 232 certified map edges. Source/site JSON copies are byte
  identical; all 23 criterion IDs are unique; the one SIKEp434 protocol case
  is separate from graph nodes and WID1 findings. The finding count remains
  29. `cargo fmt --all --check`, `cargo clippy --release --locked --offline
  -- -D warnings`, and `node --check` for both site scripts passed; HTML IDs
  and four atlas tabs are unique and present; `git diff --check` passed.
- `lake build MinimalPairingField` and `lake env lean AxiomAudit.lean`
  passed on the pinned Lean 4.34.1/Mathlib toolchain. The three new
  declarations prove the order formula, divisibility into the pairing
  extension, and minimality of the multiplicative order. Their axiom audit
  reports only `propext`, `Classical.choice`, and `Quot.sound`; the
  finite-field roots-of-unity identification remains an explicitly open
  bridge in [PROOF_MAP.md](PROOF_MAP.md). `lake env leanchecker --fresh
  --verbose MinimalPairingField` also exited successfully with output
  `replaying MinimalPairingField with --fresh`.
- The Sites source helper pushed final commit
  `e4ffc9e3eab84c56b1cb2b977c9d9598dfadb528`, packaged
  `/Volumes/SSD990/llm/tmp/weakcurves-static-v13.tar.gz`, saved version 13
  (`appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_4381e2a6df4c81918edc52e0ff2efc1a`),
  and deployed publicly as `appgdep_6acd4c47cc38819185b09773d03fa731`.
  Production returned HTTP 200 and its `snapshot.json` contains 105 named
  curves and one protocol case at
  <https://weakcurves.aburan28.chatgpt.site>. The public audience was
  preserved. Browser pointer-interaction QA was unavailable in this run.

## Vulnerability filters and isogeny graph, 12 October 2026

- `python3 research/weak_curve_leaderboard_20261012/build_isogeny_graph.py --check`
  passed. The builder checked the P-192 walk's `VERIFY.json` digest and
  71-node/132-edge counts, all 100 large-prime map artifact SHA-256 digests,
  endpoint field/trace consistency, and unique edge IDs. The deterministic
  export has 1,090 nodes, 674 classes, 73 multi-model classes, and 232
  certified directed map edges. Of its nodes, 920 are atlas records and 170
  are additional verified-walk identities. All source digests are embedded
  in `isogeny_graph.json`.
- A separate Python integrity check verified every node/class membership,
  every edge's class equality, 319 atlas nodes in multi-model classes and
  172 graph nodes incident to certified maps. The 16 criterion IDs are
  unique. Source and site copies of `criteria.json` and `isogeny_graph.json`
  are byte-identical. An HTML parser found unique element IDs and both
  JavaScript entrypoints and filter controls.
- `cargo test --release --locked --offline` in the leaderboard crate:
  four passed. `cargo run --release --locked --offline -- check snapshot.json`
  passed. `node --check` passed for `app.js` and `graph.js`; `git diff --check`
  passed. Browser-control infrastructure was unavailable, so visual and
  pointer-interaction checks were not performed in a browser.
- `CARGO_NET_OFFLINE=true cargo test --release --lib --quiet` at repository
  root completed in 632.33 seconds with **3,991 passed, nine failed, 113
  ignored**. The nine failures are the same existing
  `cryptanalysis::pollard_collab::{cairn,net}` local-socket permission cases
  recorded in [root_lib_v3.log](validation/root_lib_v3.log) and the previous
  [root-suite baseline](validation/root_lib_v2.log). No Rust library source
  changed in this addition.
- The Sites source helper pushed commit
  `2400a588268a505818d2830f99b7ad0230d38f29`, packaged
  `/Volumes/SSD990/llm/tmp/weakcurves-static-v10.tar.gz`, saved version
  `appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_4551e69d3b248191adca19281129e313`,
  and deployed it publicly as `appgdep_6acd3b772e308191a1b89d9abf555e62`
  with status `succeeded` at
  <https://weakcurves.aburan28.chatgpt.site>. The existing public audience
  was read before publication and preserved.

## Extension-field coverage audit, 12 October 2026

- The source-linked [coverage audit](COVERAGE_AUDIT.md) adds three distinct
  routes to `criteria.json`: direct extension-field decomposition,
  Joux–Vitse cover-and-decomposition, and Tian's cubic Jacobian isogeny.
  The criterion count is 13; the exact-curve, WID1, and published-case
  counts remain 920, 29, and 36. The new entries assign no unverified
  weakness IDs.
- Primary source pages checked: [Joux–Vitse 2010/157](https://eprint.iacr.org/2010/157),
  [Joux–Vitse 2011/020](https://eprint.iacr.org/2011/020),
  [Tian 2012.07173](https://arxiv.org/abs/2012.07173),
  [Iijima–Momose–Chao 2015/805](https://eprint.iacr.org/2015/805),
  [2009/613](https://eprint.iacr.org/2009/613), and
  [2008/215](https://eprint.iacr.org/2008/215). The repository's
  `jacobian-transfers.json` reports nine extension-index-calculus and 312
  same-field-cover formula rows; the snapshot has 661 degree-six, 77
  degree-ten, and 76 degree-fourteen local models.
- `cargo test --release --locked --offline` in this crate: four passed.
  `cargo run --release --locked --offline -- check snapshot.json`: passed.
  JSON parsing, 13 unique criterion IDs, source/site byte parity, HTML
  section balance, and `node --check` for `app.js`: passed.
- `CARGO_NET_OFFLINE=true cargo test --release --locked --lib --quiet` at
  repository root: **3,991 passed, nine failed, 113 ignored**, exit 101.
  All nine failures are the same `pollard_collab::{cairn,net}` local-socket
  `Operation not permitted` cases recorded in [root_lib_v2.log](validation/root_lib_v2.log).
  No Rust library source changed in this audit.
- The Sites source helper opened the existing Site and pushed source commit
  `e481237f3017a641935f0efc1a2ff9a34af2761c`. It packaged
  `/Volumes/SSD990/llm/tmp/weakcurves-static-v7.tar.gz`, saved version
  `appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_b7e1e7c4c3cc8191948b33fad0d8f4c2`,
  and deployed it publicly as
  `appgdep_6acd2a72945c8191ba910d9e7194ec44` with status
  `succeeded` at <https://weakcurves.aburan28.chatgpt.site>.
  The existing public audience was preserved.


## Snapshot 02: named curves and WID1

- `cargo test --release --locked`: four focused tests passed, including
  source-locked Oakley records and WID1 stability across evidence grades.
- `cargo run --release --locked -- build snapshot.json` then `check`:
  passed. The exporter reads and hashes the 1,929-entry ICV1 registry and
  the Jacobian-transfer registry, imports 103 named curves, and asserts 920
  unique curve IDs and 29 unique WID1 findings. Two named rows are Oakley
  GHS candidates.
- `cargo clippy --release --locked -- -D warnings`,
  `cargo fmt --all --check`, `node --check website/weakcurves/dist/app.js`,
  three source/site byte comparisons and a JSON count check: passed.
- The three prime-field ICV1 IDs and slugs match independent calls to
  `scripts/curve_id.py prime` with their exact coefficients and orders;
  their former `EX1` labels remain searchable aliases. A separate Python
  SHA-256 replay recomputed all 29 WID1 IDs from canonical sorted-key JSON
  and found no duplicates.
- `/Volumes/SSD990/cryptanalysis/sage --runtime-info` is retained in
  [sage-runtime.json](validation/sage-runtime.json). Then
  `/Volumes/SSD990/cryptanalysis/sage research/weak_curve_leaderboard_20261012/validation/oakley.sage`
  passed for both RFC points: field polynomials irreducible, points on their
  stated curves, both large factors prime, and `[4]P` and `[2]P` of exact
  prime order. The [raw output](validation/oakley.stdout) is retained.
  The first validation attempt used a Sage field method absent from the
  accepted NTL backend (`fetch_int`); it stopped before any mathematical
  check and was corrected to polynomial-basis conversion.
- `lake build PrimeExamples OakleyOrders`,
  `lake env lean AxiomAudit.lean`, and
  `lake env leanchecker --fresh --verbose OakleyOrders`: passed. The
  [axiom audit](validation/lean_axioms.log) includes all seven Oakley
  declarations; their only foundations are the standard Lean
  `propext`, `Classical.choice` and `Quot.sound` where needed. The
  [fresh checker output](validation/oakley_lean_fresh.log) records replay.
- `CARGO_NET_OFFLINE=true cargo test --release --lib`: 3,991 passed, nine
  failed, 113 ignored. The nine failures are the same pre-existing
  `pollard_collab::{cairn,net}` socket-permission cases recorded below.
  The full [new raw log](validation/root_lib_v2.log) retains the failures.
- The final Sites source was pushed at
  `8777eaa7657b33a3f07dbd3f90b01039b356180c`, packaged as
  `weakcurves-static-v6.tar.gz`, and saved as version
  `appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_cb3b9b1abd4881919dd08d5851a54ca6`.
  The current public access policy was read before publication. Deployment
  `appgdep_6acd26c00fec8191ab129a99cde4aa05` returned `succeeded`
  at <https://weakcurves.aburan28.chatgpt.site>. Version 5 was an
  intermediate successful deployment before the prime examples' ICV1
  identities were reconciled; version 6 supersedes it.
- `python3 scripts/check_curve_names.py` stopped before linting because
  this sparse worktree lacks the tracked
  `research/boolean_affine_substitution_20260923/ARCHIVE_NOTES.md` file.
  This is a checkout-materialization limit; the new ICV1 values were
  cross-checked with the reference implementation above.

The full RFC curve orders remain source-declared in the inherited standards
registry. Sage independently verified the published generator orders, and
Lean proves the arithmetic and the conditional group-order implication.
The GHS cover, subgroup map and qualified cost remain separate open gates.

## Source and export

- `cargo test --release --locked --offline`: two focused tests passed.
- `cargo run --release --locked --offline -- build snapshot.json` and the
  matching `check`: passed after the final large-field screen join.
- `cargo fmt --all --check` and
  `cargo clippy --release --locked --offline -- -D warnings`: passed after a
  local lint correction; the frozen snapshot still passes `check`.
- The exporter asserts 814 exact source models, 24 positive model-cover
  signals, 537 large-field models, and large-field screen counts 247 / 137 /
  153. It also asserts all 36 published parameter cases are imported.
- `research/weak_curve_leaderboard_20261012/snapshot.json` and
  `website/weakcurves/dist/snapshot.json` are byte-identical; the same holds
  for `criteria.json`.

## Independent prime-field replay

`gp -q validation/pari.gp` passed. The retained [output](validation/pari.stdout)
gives `#E=101`, point order 101 for the anomalous model; `#E=100`, group
structure `[10,10]`, point order 10 for the smooth model; and `#E=104`, point
order 13 for the pairing candidate. The Rust enumerator and group law are
independent of PARI. `103 % 13=12` and `103² % 13=1` give embedding degree 2.

## Lean

- Lean 4.34.1, Mathlib commit `d13f23b723b8a846827a245b89c10fc7d3f11612`.
- `lake build PrimeExamples`: passed, 690 jobs; actual theorem module built.
- `lake env lean AxiomAudit.lean`: passed for all twelve declarations. The
  complete [axiom output](validation/lean_axioms.log) contains only `propext`,
  `Classical.choice`, and `Quot.sound` among the three point-count results;
  the other arithmetic results use subsets of those standard foundations.
- `lake env leanchecker --fresh --verbose PrimeExamples`: passed; the retained
  [fresh replay output](validation/lean_fresh.log) confirms kernel replay.

The proof map identifies the open elliptic-group and attack-reduction bridge.
No declaration uses `sorry`, `admit`, `sorryAx`, a custom research axiom, or
unverified native evaluation.

## Website and repository

- `node --check website/weakcurves/dist/app.js`: passed.
- JSON shape and final 247 / 137 / 153 counts: passed with `jq -e`.
- Snapshot 01 hosted Site source commit:
  `bd71035b2197a4c4db4193be1c5dcd121500de78`.
- Sites saved and deployed the packaged static source as version
  `appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_ce4ad015b5d08191adf58403ac5e288b`.
  The deployment returned `succeeded` at
  <https://weakcurves.aburan28.chatgpt.site>.
- A browser-control tool was unavailable in this environment. The Sites
  deployment succeeded and generated a server-side screenshot; JavaScript
  syntax and data files were checked locally.
- `CARGO_NET_OFFLINE=true cargo test --release --lib`: completed with
  **3,991 passed, 9 failed, 113 ignored**. The full [raw log](validation/root_lib.log)
  shows all nine failures in existing `cryptanalysis::pollard_collab::{cairn,net}`
  tests, each blocked by an `Operation not permitted` socket error in this
  sandbox. The same nine failures appear in the previous merged
  [case-index root-suite baseline](../weak_curve_case_index_20261012/root_cargo_test.log.gz);
  that earlier run also had two additional unrelated failures. No touched
  exporter or site path is involved in the nine failures.

The custom hostname `weakcurves.com` was pending DNS validation in
snapshot 01. The current Site URL and access mode are recorded above.

The Conductor task template also lists `go vet ./...` and `go test ./...`.
Both are inapplicable at this Rust repository root because it has no
`go.mod`; the commands report `directory prefix . does not contain main
module`. The Rust checks above are the relevant validation for this change.
