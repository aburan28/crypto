import CensusCounting
import Lean.Util.CollectAxioms
import Lean.Elab.Command

open Lean Elab Command

run_cmd do
  let allowed := #[`propext, `Classical.choice, `Quot.sound]
  for name in #[`CensusCounting.inclusion_exclusion,
      `CensusCounting.union_zero_partition, `CensusCounting.union_bounded_by_population] do
    let axioms ← Lean.collectAxioms name
    unless (axioms.filter fun ax => !allowed.contains ax).isEmpty do
      throwError "Forbidden dependency in {name}: {axioms}"
    logInfo m!"AXIOM_AUDIT|{name}|{axioms}"
  logInfo "AXIOM_AUDIT_COMPLETE|public_theorems=3"
