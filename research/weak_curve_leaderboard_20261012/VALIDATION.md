# Validation record

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
