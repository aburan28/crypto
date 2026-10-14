import Init.Grind.Ring.Basic

/-!
The third summation polynomial is unchanged by exchanging its two indexed
inputs. Any zero-preserving semiring endomorphism transports its roots.
The link from the source's N83 Frobenius rotation and canonical-root encoding
to these algebraic hypotheses is tracked separately in PROOF_MAP.md.
-/

namespace N83SymmetricRoot

open Lean.Grind

variable {F : Type} [CommSemiring F]

def s3 (x y z : F) : F :=
  (x * y + x * z + y * z) * (x * y + x * z + y * z) + x * y * z

theorem s3_swap (x y z : F) : s3 x y z = s3 y x z := by
  have hsum : x * y + x * z + y * z = y * x + y * z + x * z := by
    rw [CommSemiring.mul_comm y x]
    rw [Semiring.add_assoc, Semiring.add_comm (x * z) (y * z)]
    rw [← Semiring.add_assoc]
  simp only [s3, hsum, CommSemiring.mul_comm y x]

theorem s3_swap_root {x y z : F} (h : s3 x y z = 0) : s3 y x z = 0 := by
  rw [s3_swap] at h
  exact h

theorem s3_map (φ : F → F)
    (hadd : ∀ a b, φ (a + b) = φ a + φ b)
    (hmul : ∀ a b, φ (a * b) = φ a * φ b)
    (x y z : F) : φ (s3 x y z) = s3 (φ x) (φ y) (φ z) := by
  simp only [s3, hadd, hmul]

theorem s3_map_root (φ : F → F)
    (hzero : φ 0 = 0)
    (hadd : ∀ a b, φ (a + b) = φ a + φ b)
    (hmul : ∀ a b, φ (a * b) = φ a * φ b)
    {x y z : F} (h : s3 x y z = 0) :
    s3 (φ x) (φ y) (φ z) = 0 := by
  rw [← s3_map φ hadd hmul, h, hzero]

/- The following arithmetic uses the source-loop upper bound for retained
   states.  The quotient's actual regular-state count can be smaller. -/
def quotientStateBound (k : Nat) : Nat :=
  83 * (k * (k - 1) / 2) + 42 * k

theorem k2543_fits_five_gibibyte_bucket :
    4 * quotientStateBound 2543 ≤ 2 ^ 30 := by decide

theorem k2543_state_index_fits_twenty_nine_bits :
    quotientStateBound 2543 + 1 < 2 ^ 29 := by decide

theorem k2544_crosses_five_gibibyte_bucket :
    2 ^ 30 < 4 * quotientStateBound 2544 := by decide

theorem k2543_state_and_slot_payload_under_twelve_gibibytes :
    28 * quotientStateBound 2543 + 5 * (2 ^ 30) ≤
      12 * (2 ^ 30) := by decide

#print axioms s3_swap
#print axioms s3_swap_root
#print axioms s3_map
#print axioms s3_map_root
#print axioms k2543_fits_five_gibibyte_bucket
#print axioms k2543_state_index_fits_twenty_nine_bits
#print axioms k2544_crosses_five_gibibyte_bucket
#print axioms k2543_state_and_slot_payload_under_twelve_gibibytes

end N83SymmetricRoot
