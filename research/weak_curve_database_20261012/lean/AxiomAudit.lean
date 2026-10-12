import WeakCurveDatabase
import Lean.Util.CollectAxioms

open Lean Elab Command

run_cmd do
  let allowed : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  for name in #[`WeakCurveDatabase.cubic_cocycle_norm_one,
      `WeakCurveDatabase.cubic_cocycle_power,
      `WeakCurveDatabase.fourth_power_of_order_three_mod_four,
      `WeakCurveDatabase.odd_square_mod_four,
      `WeakCurveDatabase.cubic_torus_fourth_power,
      `WeakCurveDatabase.cubic_torus_fourth_power_prime_square] do
    let axioms ← Lean.collectAxioms name
    let forbidden := axioms.filter fun ax => !allowed.contains ax
    unless forbidden.isEmpty do
      throwError "Forbidden proof dependency in {name}: {forbidden}"
    logInfo m!"AXIOM_AUDIT|{name}|{axioms}"

#print axioms WeakCurveDatabase.cubic_cocycle_norm_one
#print axioms WeakCurveDatabase.cubic_cocycle_power
#print axioms WeakCurveDatabase.fourth_power_of_order_three_mod_four
#print axioms WeakCurveDatabase.odd_square_mod_four
#print axioms WeakCurveDatabase.cubic_torus_fourth_power
#print axioms WeakCurveDatabase.cubic_torus_fourth_power_prime_square
