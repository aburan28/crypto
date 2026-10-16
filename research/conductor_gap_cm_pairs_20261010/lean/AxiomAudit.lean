import GenericInvariance
import TransferBound
import Lean.Util.CollectAxioms

/-!
Axiom audit for `GenericInvariance.lean` and `TransferBound.lean`. Run
`lake build && lake env lean AxiomAudit.lean` in this directory. The
`run_cmd` block fails the run if any listed declaration depends on an axiom
outside `propext`, `Classical.choice` and `Quot.sound`, for example
`sorryAx`. The `#print axioms` lines are the log that CI keeps
(`.github/workflows/conductor-gap-lean.yml`).
-/

open Lean Elab Command

run_cmd do
  let allowed : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  for name in #[`GenericInvariance.decide_eq_map,
      `GenericInvariance.step_map,
      `GenericInvariance.run_map,
      `GenericInvariance.dlog_invariant,
      `TransferBound.Hom.map_pow,
      `TransferBound.transfer,
      `TransferBound.transfer_cost,
      `TransferBound.transfer_answer,
      `TransferBound.within_isogeny_cost] do
    let axioms ← Lean.collectAxioms name
    let forbidden := axioms.filter fun ax => !allowed.contains ax
    unless forbidden.isEmpty do
      throwError "Forbidden proof dependency in {name}: {forbidden}"
    logInfo m!"AXIOM_AUDIT|{name}|{axioms}"

#print axioms GenericInvariance.decide_eq_map
#print axioms GenericInvariance.step_map
#print axioms GenericInvariance.run_map
#print axioms GenericInvariance.dlog_invariant
#print axioms TransferBound.Hom.map_pow
#print axioms TransferBound.transfer
#print axioms TransferBound.transfer_cost
#print axioms TransferBound.transfer_answer
#print axioms TransferBound.within_isogeny_cost
