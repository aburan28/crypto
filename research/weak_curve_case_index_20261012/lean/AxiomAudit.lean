import PublishedCases
import Lean.Util.CollectAxioms

open Lean Elab Command

run_cmd do
  let allowed : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  for name in #[`PublishedCases.cubic_frobenius_factorization,
      `PublishedCases.septic_frobenius_factorization] do
    let axioms ← Lean.collectAxioms name
    let forbidden := axioms.filter fun ax => !allowed.contains ax
    unless forbidden.isEmpty do
      throwError "Forbidden proof dependency in {name}: {forbidden}"
    logInfo m!"AXIOM_AUDIT|{name}|{axioms}"

#print axioms PublishedCases.cubic_frobenius_factorization
#print axioms PublishedCases.septic_frobenius_factorization
