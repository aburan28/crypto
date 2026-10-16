import RhoFormal

/-! Axiom audit: run with `lake env lean RhoFormal/Axioms.lean`.  Not part of
the default build target. -/

open RhoFormal.LeastRotation in
section
#print axioms val_lt_iff_lt
#print axioms lt_of_lz_lt
#print axioms least_rotation_has_max_lz
#print axioms cand_nonempty
#print axioms least_is_cand
#print axioms cand_min_is_least
#print axioms exists_cand_min
#print axioms least_rotation_via_candidates
#print axioms rot_const
#print axioms rot_all_false
#print axioms rot_all_true
#print axioms const_every_k_least
#print axioms val_rot_eq_rust
#print axioms getElem_rot
#print axioms zmask_spec
#print axioms rust_candidates_exact
end

open RhoFormal.BatchInv in
section
#print axioms batchInv_correct
#print axioms batchInv_eq_map
#print axioms batchInv_eq_map_inv
#print axioms combine4_eq_inv
#print axioms batchInv4_eq_map
#print axioms batchInv4_correct
end

open RhoFormal.RhoWalk in
section
#print axioms jump_step
#print axioms canon_step
#print axioms canon_compose
#print axioms canon_neg
#print axioms escape_step
#print axioms collision
#print axioms candidate_check
#print axioms collision_same_b
#print axioms denominator_zero_iff
#print axioms cursor_init
#print axioms cursor_advance
#print axioms rho_invariant
end
