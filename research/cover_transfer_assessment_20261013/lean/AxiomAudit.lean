import CoverTransfer
import Lean.Util.CollectAxioms

open Lean Elab Command

run_cmd do
  let theoremName := `CoverTransfer.pull_injective_on_subgroup
  let axioms ← Lean.collectAxioms theoremName
  let permitted : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  let unexpected := axioms.filter fun ax ↦ !permitted.contains ax
  unless unexpected.isEmpty do
    throwError m!"Forbidden proof dependency in {theoremName}: {unexpected}"
  logInfo m!"AXIOM_AUDIT|{theoremName}|{axioms}"

#print axioms CoverTransfer.pull_injective_on_subgroup
