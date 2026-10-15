# Validation record

## Genus-128 GHS model-stage audit, 15 October 2026

- PR #1738 head `ae67d4a2f2` exposed a clean Linux runner mismatch in the
  pre-existing catalog snapshot check: the six Rust tests passed, then
  `weak-curve-leaderboard check snapshot.json` reported differing bytes. The
  completed job log is preserved in
  `validation/ghs128_ci_snapshot_mismatch_attempt1.log` (GitHub run
  `38710405870`, job `116340086585`). Local source-file hashes match the
  frozen snapshot and the same check passes on macOS. The exporter now reports
  the first differing byte and short contexts so the next clean run can
  isolate the field without replacing the snapshot on a guess. Its
  deliberately altered-snapshot control is in
  `validation/ghs128_mismatch_diagnostic_control.log`.
- The mismatch was traced to newer `main` copies of the curve registry and
  Jacobian-transfer registry in the PR merge checkout. Both were copied
  byte-for-byte into the feature branch. The rebuilt snapshot changes only
  their two source hashes; the native GHS audit changes only the transfer
  registry hash. The graph and decomposition exports change their pinned
  snapshot hash, with the same node, class, edge, property and finding counts.
  The [refresh run](reports/run-20261015T0929Z/REPORT.md) preserves the
  quantitative comparison and inspected PDF. Current-head clean CI remains
  the final gate.
- The next clean runner passed the snapshot, six catalog tests, native GHS
  audit, and attack-map validator, then stopped at the decomposition index's
  `covers_sha256` comparison. Its raw job log is
  `validation/ghs128_ci_decomposition_source_attempt2.log` (run
  `38713758398`, job `116351297610`). The current-main
  `docs/curves/covers.json` was copied byte-for-byte and the research and
  site decomposition indexes regenerated; the [alignment report](reports/run-20261015T0950Z/REPORT.md)
  shows the source dependency. All other declared source paths for the
  catalog, graph and GHS audit were compared with `origin/main` at
  `f2c9202ca48ed16f086815fda9c0ab73f2aabcda` and matched or were
  branch-local rosters. The current-head clean replay remains pending.

- `cargo run --release --locked --bin ghs_genus128_cost_audit -- --check` recomputes exact
  irreducible-polynomial and smoothness counts for all 84 degree bounds in
  the declared `GF(4)`, genus-128 search envelope. The optimum is bound 15
  for both `c2pnb176w1` and `c2pnb272w1`, with 48,921,052 nominal half-split
  places, smoothness log2 fraction -26.328277, and model relation-query
  log2 count 51.872229. It independently reproduces the transfer registry's
  71.336 and 71.774 log2 destination prices from the pinned source cost
  parameters. The receipt retains every candidate price and the exact
  smoothness numerator, with null fields for curve-specific arithmetic and
  measured relation stages.
- The binary GHS criterion now requires actual usable base classes,
  ordinary-query relation yield, independent rank, and complete cost after
  the type/genus and subgroup-transfer checks. The cost audit checks a
  conditional model; it does not instantiate those later witnesses.
- The database workflow now replays the catalog exporter, both property
  validators, this exact-integer cost audit, and the versioned artifact hashes
  on a clean runner. Its result is tracked separately from local validation.
- The native source, receipt, cost-sweep SVG, report source and inspected
  one-page PDF are pinned in `ARTIFACTS.sha256`. Local
  `cargo test --release --locked --offline` passed six tests (two GHS-audit
  tests and four catalog tests); the native `--check`, catalog snapshot check,
  attack-map validator, decomposition validator, site JavaScript syntax check,
  and artifact hash check passed. Logs for native build/write/replay are in
  `validation/ghs128_native_*.log`; the first replay command's wrong-working-
  directory failure is retained separately from the successful study-local
  replay.

## Oakley Type II genus reconciliation, 15 October 2026

- The current binary transfer registry classifies Oakley Groups 3 and 4 as
  Type II: each has magic number 5, nonzero relative trace of `b`, and
  projected genus 16. The earlier leaderboard snapshot asserted genus 15
  and carried stale formula costs. The regenerated snapshot now reads the
  current route, records genus 16, updates both cost projections, and gives
  the corrected GHS findings `WID1:GHS:56f762afe456` and
  `WID1:GHS:4df6967b110c`. Each finding's `corrects` field names its old
  genus-15 ID. The graph and per-curve decomposition exports were regenerated
  from that snapshot and copied into the bundled site data. The binary GHS
  criterion now asks for the type and trace certificate before assigning a
  projected genus.
- `CARGO_TARGET_DIR=/private/tmp/weakcurves-oakley-target cargo test --release
  --locked --offline` passed all four exporter tests. The same environment's
  `cargo run --release --locked --offline -- check snapshot.json` passed.
  Both graph and decomposition builders passed `--check`; the attack-map and
  decomposition validators passed with 922 atlas identities, 1,092 graph
  nodes, 232 directed certified maps, and 13 decomposition filters.
- `lake build OakleyOrders`, `lake env leanchecker --fresh --verbose
  OakleyOrders`, and `lake env lean AxiomAudit.lean` passed. The new conditional
  formula theorem reports only `propext`; its model-specific Type II premise
  comes from the binary-field transfer computation. Build, fresh replay,
  and axiom outputs are retained in `validation/oakley_genus_*.log`.
- The preceding PR head's `leaderboard-lean-proof` CI job completed its
  RelationSupport build and fresh replay, then failed because it had not built
  `PrimeExamples`, which the shared axiom audit imports. The workflow now
  builds all four imported theorem modules before either fresh replay or
  the audit. The exact four-target `lake build` command passed locally and
  its output is retained in `validation/leaderboard_all_build.log`; the new
  head's CI result must be checked separately.

## Fixed-base relation support, 15 October 2026

- `lake build RelationSupport` passed on pinned Lean 4.34.1 and Mathlib
  `d13f23b723b8a846827a245b89c10fc7d3f11612`.
  `lake env lean AxiomAudit.lean` checked the generic unordered-sum image
  and accepted-set intersection bounds. Both declarations report only
  `propext`, `Classical.choice`, and `Quot.sound`. `lake env leanchecker
  --fresh --verbose RelationSupport` passed. The build, axiom and fresh
  replay outputs are retained in `validation/relation_support_*.log`.
- `python3 validate_attack_map.py` checks 24 indexed criteria and 38 WAM1
  research checks, including byte equality of the research and site criterion
  copies. The new criterion records the P-192 fixed-base support cap alongside
  the independently verified 64 accepted divisor-pair representations,
  zero three-sum relations, and rank zero.
- `cargo run --release --locked --offline -- check snapshot.json` compiled
  but stopped at the existing Oakley source/route assertion: the current
  `docs/curves/jacobian-transfers.json` records destination genus 16 for
  Groups 3 and 4, while `src/main.rs` still expects 15. Those source and
  exporter files are unchanged in this criterion addition. The broader
  snapshot check is therefore unresolved pending reconciliation of that
  separate source drift.

## Joux–Vitse intermediate-field route, 13 October 2026

- The source-linked `2011/020` method now has separate WAM1 checks for the
  intermediate-field cover and the odd-characteristic quadratic sieve. The
  indexed criterion specifies the norm–pullback subgroup kernel and both
  sextic towers. The 922 curve identities each carry 13 property records:
  14 cubic model labels, 647 additional sextic tower screens, and 261 records
  outside these screened towers. Every lower-field cover, subgroup transport,
  saturated sieve ideal, natural relation rank and target descent remains an
  explicit open field. The genus-two `F_(p³)` Scholten branch is unassessed.
- The [published case record](published_cover_cases.json) preserves the
  paper's `p=33,554,467` field tower, alpha coefficients, 149-bit selected
  prime, factor-base size and source-reported phase measurements. The
  validator checks its field-size and order arithmetic, coefficient bounds,
  and byte-identical site copy. It has no local atlas identity or independent
  replay.
- `python3 build_decomposition_properties.py --check`,
  `python3 validate_decomposition_properties.py`,
  `python3 validate_attack_map.py`, and `node --check
  website/weakcurves/dist/app.js` passed. The attack map has 23 indexed
  criteria and 38 research checks. This addition classifies a published
  method and records observations; it introduces no new mathematical theorem
  requiring a separate Lean declaration.
- The Sites source helper pushed commit
  `51711fc6688c43011704a108d8d7a9cca43b74a4` and packaged
  `/Volumes/SSD990/llm/tmp/weakcurves-static-v17.tar.gz`. Sites saved
  version 17 and returned `succeeded` for public deployment
  `appgdep_6acd8bcceafc81919c217e1d1f686abe` at
  <https://weakcurves.aburan28.chatgpt.site>. The existing public audience
  was preserved.

## Curve-level hyperelliptic detection properties, 12 October 2026

- `python3 build_decomposition_properties.py --check` and
  `python3 validate_decomposition_properties.py` passed. The deterministic
  export binds 922 unique atlas IDs to eleven property families. It imports
  917 `verified_over_declared_field` cover certificates from
  `docs/curves/covers.json` by exact ICV1 identity, slug and full model hash.
  The other five prime-field records have separately labeled formula replays:
  the exporter checks the substituted sextic identity, squarefreeness, and
  source nonsingularity. The copied site JSON is byte-identical. Current
  property counts include 14 cubic and ten quadratic positive genus-three
  model labels, 247 cubic necessary-screen survivors, 49 curves with exact
  GHS invariant records, 37 size-admissible high-genus projections and 103
  same-field formula-cost rows. Divisor relations, Jacobian factorization,
  explicit subgroup transport and complete measured route cost remain
  individually recorded as open fields rather than inferred from geometry.
- `python3 validate_attack_map.py` passed with 23 criteria, 36 method checks,
  1,092 graph nodes and 232 directed certified maps. The isogeny graph
  builder `--check` retained 922 atlas nodes and 676 classes. Python bytecode
  compilation of both new scripts, `node --check` on the updated `app.js`,
  an HTML parser check of 51 unique atlas IDs and four guide IDs with all
  local assets present, and `git diff --check` passed. No new mathematical
  theorem is asserted by the property join; the structural Lean bridge
  remains described in [PROOF_MAP.md](PROOF_MAP.md).
- The Sites source helper pushed commit
  `3cff28af0e42580cf0862f53b14cf0da5d205177`, packaged
  `/Volumes/SSD990/llm/tmp/weakcurves-static-v16.tar.gz`, saved version 16
  (`appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_7c4a88ba0fd08191bba85a0aedcc2f12`),
  and deployed it publicly as `appgdep_6acd87c44c8c81919dd585c440ff0a37`.
  Sites returned `succeeded` at <https://weakcurves.aburan28.chatgpt.site>;
  the existing public audience was preserved. The live
  `decomposition_properties.json` matched the research artifact byte for
  byte. Browser pointer-interaction QA was unavailable.

## Hyperelliptic decomposition routes, 12 October 2026

- The [method sheet](HYPERELLIPTIC_DECOMPOSITION.md) separates divisor
  decomposition for Jacobian index calculus from splitting a Jacobian up to
  isogeny. It records a third, distinct cover-transfer gate and the field of
  definition, selected-subgroup image, and complete-cost conditions. Primary
  method sources checked: [Sarkar–Singh](https://eprint.iacr.org/2014/815),
  [Gaudry–Thomé–Thériault–Diem](https://eprint.iacr.org/2004/153),
  [Enge](https://doi.org/10.1090/S0025-5718-01-01363-1),
  [Shaska](https://arxiv.org/abs/1209.3187), and
  [Ahmadi–McGuire–Rojas-León](https://arxiv.org/abs/1409.8510).
- `python3 validate_attack_map.py` passed with 23 indexed criteria and 36
  research checks, five of them new hyperelliptic gates. The total searchable
  checks are 59; the existing 29 findings and 922-model atlas are unchanged.
  `python3 build_isogeny_graph.py --check` passed with 1,092 graph nodes, 676
  classes and 232 explicit directed maps. `node --check` passed for all three
  site JavaScript files. An HTML parser found 48 unique IDs on the atlas page
  and four unique IDs on the new method page, with every local asset link
  resolving; `git diff --check` passed. The guide indexes published results;
  the existing [Lean proof map](PROOF_MAP.md) retains the formalization status
  of transfer implications.
- The Sites source helper pushed commit
  `8441d4b8d530981507bfa2ae554493157d8a5277`, packaged
  `/Volumes/SSD990/llm/tmp/weakcurves-static-v15.tar.gz`, saved version 15
  (`appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_18a4fc856218819186bf6848b26c774d`),
  and deployed publicly as `appgdep_6acd7e0c8ebc819186deb3723ff02a6e`.
  The deployment returned `succeeded` at
  <https://weakcurves.aburan28.chatgpt.site>. The public audience was
  preserved. Live `attack_map.json` matched the source byte-for-byte, and
  the new `hyperelliptic-decomposition.html` served the expected sections.
  Browser pointer-interaction QA was unavailable.

## WAM1 attack map and focused graph, 12 October 2026

- `python3 validate_attack_map.py` passed: 23 distinct indexed criteria plus
  31 source-linked research checks, 1,092 graph nodes and 232 certified
  directed maps. Every edge has both endpoints in the same finite-field
  isogeny class; the research and site copies of `attack_map.json` match
  byte-for-byte. The WAM1 entries are an unassessed method queue and change
  neither the 29 WID1 findings nor the 922-model atlas.
- `node --check` passed for `app.js`, `graph.js` and `attack-map.js`.
  An HTML parser found 48 unique IDs and all attack-map and graph controls.
  The local static server returned HTTP 200 for `index.html`, the new JavaScript
  and JSON assets, and the graph data. The graph builder `--check` passed with
  676 classes and 232 map edges. Browser pointer-interaction QA was unavailable.
- Primary sources were checked for the new input, auxiliary-information,
  implementation and protocol families, including [RFC 7748](https://www.rfc-editor.org/rfc/rfc7748),
  [RFC 9380](https://www.rfc-editor.org/rfc/rfc9380),
  [Cheon](https://doi.org/10.1007/11761679_1),
  [Gorla–Massierer](https://eprint.iacr.org/2014/318),
  [Neves–Tibouchi](https://eprint.iacr.org/2015/1233), and
  [Galbraith–Petit–Shani–Ti](https://eprint.iacr.org/2016/859).
- The Sites source helper pushed commit `f8e9fd15cc364207222637fa770d01afae7157e1`,
  packaged `/Volumes/SSD990/llm/tmp/weakcurves-static-v14.tar.gz`, saved
  version 14 (`appgprj_6acd13c887f0819194ae65d7f61c7b33~appgver_97c5908652448191903dbda1eed35e75`),
  and deployed publicly as `appgdep_6acd542272908191aec104e80d8a4726`.
  The deployment returned `succeeded` at
  <https://weakcurves.aburan28.chatgpt.site>; the existing public audience was
  preserved. Live `attack_map.json` and `graph.js` match the packaged files
  byte-for-byte; the live HTML contains the new controls and receives a
  Cloudflare challenge script at the end of the response. This update contains
  no new mathematical theorem, so the pinned
  Lean proof map and checked declarations are unchanged.

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
