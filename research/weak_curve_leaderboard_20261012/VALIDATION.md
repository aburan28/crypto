# Validation record

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
- Hosted private Site source commit:
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

The custom hostname `weakcurves.com` is pending DNS validation; the generated
Site URL above is live and private.

The Conductor task template also lists `go vet ./...` and `go test ./...`.
Both are inapplicable at this Rust repository root because it has no
`go.mod`; the commands report `directory prefix . does not contain main
module`. The Rust checks above are the relevant validation for this change.
