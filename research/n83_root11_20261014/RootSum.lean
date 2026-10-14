import Std.Tactic

/-! Bit reconstruction and capacity obligations for one retained N83 root. -/

namespace N83RootSum

open Lean.Grind

def stateBound (k : Nat) : Nat := 83 * (k * (k - 1) / 2) + 42 * k

theorem xor_reconstruct (first second : Nat) :
    first ^^^ (first ^^^ second) = second := by
  rw [← Nat.xor_assoc, Nat.xor_self, Nat.zero_xor]

/-- Cancel the nonzero squared product in the regular S3 solver's `q`.
    The solver's `a` is `(x+y)^2` and `p` is `xy`. -/
theorem regular_root_delta {α : Type} [Field α] (a p : α)
    (hp : p ≠ 0) :
    p * (p * p) * (a * (p * p))⁻¹ = p * a⁻¹ := by
  have hpp : p * p ≠ 0 := by
    intro h
    cases Field.of_mul_eq_zero h with
    | inl hzero => exact hp hzero
    | inr hzero => exact hp hzero
  calc
    p * (p * p) * (a * (p * p))⁻¹
        = p * (p * p) * (a⁻¹ * (p * p)⁻¹) := by rw [Field.inv_mul]
    _ = (p * a⁻¹) * ((p * p) * (p * p)⁻¹) := by grind
    _ = p * a⁻¹ := by rw [Field.mul_inv_cancel hpp, Semiring.mul_one]

theorem k3400_complete_key_bound_fits_fixed_slots :
    2 * stateBound 3400 < 2 ^ 30 := by decide

theorem k3400_state_id_fits_29_bits :
    stateBound 3400 + 1 < 2 ^ 29 := by decide

theorem k3400_one_root_and_slot_payload_below_10_gibibytes :
    11 * stateBound 3400 + 5 * (2 ^ 30) < 10 * (2 ^ 30) := by decide

#print axioms xor_reconstruct
#print axioms regular_root_delta
#print axioms k3400_complete_key_bound_fits_fixed_slots
#print axioms k3400_state_id_fits_29_bits
#print axioms k3400_one_root_and_slot_payload_below_10_gibibytes

end N83RootSum
