import TreeSplit.Verdict
import TreeSplit.Replay
import TreeSplit.Leak
import Lean.Util.CollectAxioms
import Lean.Elab.Command

/-!
Axiom audit for the `TreeSplit` theorems (sheriff case SHF-20261016-3caa6f,
finding L3). Run `lake build && lake env lean AxiomAudit.lean` in this
directory; CI does (`.github/workflows/tree-split-structured-v-lean.yml`).

- The kernel-checked theorems (`Verdict.lean`, `Replay.targets_count`,
  `Leak.lean`) may depend on `propext`, `Classical.choice` and `Quot.sound`
  only. Anything else, `sorryAx` or a compiler-trust axiom included, fails
  the run.
- `Replay.replay_ok` is the one theorem closed by `native_decide`, which
  trusts the Lean compiler. It may depend, besides the three standard
  axioms, on its own compiler-trust axioms and nothing else: the
  `TreeSplit.Replay.replay_ok._native.*` axioms Lean 4.34 generates for
  `native_decide`, or `Lean.ofReduceBool` / `Lean.trustCompiler`, which
  other Lean versions use instead. It must depend on at least one of them,
  so the log cannot present it as kernel-checked. `RESULTS.md` calls it a
  compiler-trusted replay, not a proof.

Each audited theorem logs one `AXIOM_AUDIT|<trust>|<name>|<axioms>` line,
which CI counts and keeps. The `#print axioms` lines are the same facts in
Lean's own format.
-/

open Lean Elab Command

run_cmd do
  let standard : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  let fmt := fun (axioms : Array Name) => ", ".intercalate (axioms.toList.map (·.toString))
  for name in #[`TreeSplit.cells_count,
      `TreeSplit.registered_report,
      `TreeSplit.Replay.targets_count,
      `TreeSplit.Leak.in_gf2_7,
      `TreeSplit.Leak.not_in_gf2_3] do
    let axioms ← Lean.collectAxioms name
    let forbidden := axioms.filter fun ax => !standard.contains ax
    unless forbidden.isEmpty do
      throwError "Forbidden proof dependency in {name}: {forbidden}"
    logInfo m!"AXIOM_AUDIT|kernel|{name}|{fmt axioms}"
  -- The one compiler-trusted theorem: its compiler-trust axioms only.
  let name := `TreeSplit.Replay.replay_ok
  let axioms ← Lean.collectAxioms name
  let trust := axioms.filter fun ax => !standard.contains ax
  let ownNative := fun (ax : Name) =>
    ax.toString.startsWith "TreeSplit.Replay.replay_ok._native." ||
      ax == `Lean.ofReduceBool || ax == `Lean.trustCompiler
  if trust.isEmpty then
    throwError "{name} no longer depends on a compiler-trust axiom: audit it as kernel-checked"
  unless trust.all ownNative do
    throwError "Forbidden proof dependency in {name}: {trust.filter fun ax => !ownNative ax}"
  logInfo m!"AXIOM_AUDIT|compiler|{name}|{fmt axioms}"

#print axioms TreeSplit.cells_count
#print axioms TreeSplit.registered_report
#print axioms TreeSplit.Replay.targets_count
#print axioms TreeSplit.Leak.in_gf2_7
#print axioms TreeSplit.Leak.not_in_gf2_3
#print axioms TreeSplit.Replay.replay_ok
