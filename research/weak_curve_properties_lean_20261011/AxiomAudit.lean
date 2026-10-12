import WeakCurveProperties
import Lean.Util.CollectAxioms

open Lean Elab Command

/- Audit all theorems in the proof namespace, including generated helpers.
This uses the same transitive axiom collector as `#print axioms`. -/
run_cmd do
  let env ← getEnv
  let declarations := env.constants.toList.foldl (init := #[]) fun names entry =>
    if "WeakCurveProperties.".isPrefixOf entry.1.toString then
      match entry.2 with
      | .thmInfo _ => names.push entry.1
      | _ => names
    else names
  if declarations.isEmpty then
    throwError "No curve-property theorems found"
  let permitted : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  for name in declarations.qsort Name.quickLt do
    let axioms ← Lean.collectAxioms name
    let unexpected := axioms.filter fun ax => !permitted.contains ax
    unless unexpected.isEmpty do
      throwError "Forbidden proof dependency in {name}: {unexpected}"
    logInfo m!"AXIOM_AUDIT|{name}|{axioms}"
  logInfo m!"AXIOM_AUDIT_COMPLETE|declarations={declarations.size}"

#print axioms WeakCurveProperties.conductor_floor
#print axioms WeakCurveProperties.conductor_full_two
#print axioms WeakCurveProperties.conductor_full_four
#print axioms WeakCurveProperties.fixes_torsion_iff_factors_difference
#print axioms WeakCurveProperties.finite_cyclic_twist_descent
#print axioms WeakCurveProperties.finite_field_subfield_twist_test
#print axioms WeakCurveProperties.short_isomorphism_iff
#print axioms WeakCurveProperties.generic_twist_isomorphism_iff
#print axioms WeakCurveProperties.quartic_twist_isomorphism_iff
#print axioms WeakCurveProperties.sextic_twist_isomorphism_iff
#print axioms WeakCurveProperties.generic_normalization
#print axioms WeakCurveProperties.normalized_generic_descent_iff
#print axioms WeakCurveProperties.torsion_order_depth_bounds_of_bridges
