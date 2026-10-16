import S3Roots
import Lean.Util.CollectAxioms
import Lean.Elab.Command

/-!
Axiom audit for `S3Roots.lean` (Proposition 1 of
`rootfinding_decomposition.tex`). Run `lake exe cache get && lake build &&
lake env lean AxiomAudit.lean` in this directory. The `run_cmd` block fails
the run if any listed declaration depends on an axiom outside `propext`,
`Classical.choice` and `Quot.sound`, for example `sorryAx`. The
`#print axioms` lines are the log that CI keeps
(`.github/workflows/rootfinding-s3roots-lean.yml`).
-/

open Lean Elab Command

run_cmd do
  let allowed : Array Name := #[`propext, `Classical.choice, `Quot.sound]
  for name in #[`SemaevS3.d2_mul_addX,
      `SemaevS3.S3_factor,
      `SemaevS3.S3_root_iff] do
    let axioms ← Lean.collectAxioms name
    let forbidden := axioms.filter fun ax => !allowed.contains ax
    unless forbidden.isEmpty do
      throwError "Forbidden proof dependency in {name}: {forbidden}"
    logInfo m!"AXIOM_AUDIT|{name}|{axioms}"

#print axioms SemaevS3.d2_mul_addX
#print axioms SemaevS3.S3_factor
#print axioms SemaevS3.S3_root_iff
