import CrossRatio
import Mathlib.Algebra.Group.Hom.Defs
import Mathlib.Tactic.NormNum

/-!
A sign-closed base doubles the number of available point representatives,
but each signed point uses the logarithm column of its unsigned base index.
The homomorphism is the abstract discrete-log map on the relevant subgroup;
the concrete P-192 subgroup and its divisor representatives are verified
separately.
-/

namespace P192ExtensionCover

variable {G K : Type*} [AddCommGroup G] [AddCommGroup K]

def signedBase (base : Fin B → G) (index : Fin B × Bool) : G :=
  if index.2 then -(base index.1) else base index.1

theorem signedBase_map (log : G →+ K) (base : Fin B → G)
    (index : Fin B × Bool) :
    log (signedBase base index) = signedBase (fun i => log (base i)) index := by
  cases index with
  | mk i sign => cases sign <;> simp [signedBase]

theorem signed_three_log_equation (log : G →+ K) (base : Fin B → G)
    (i j k : Fin B × Bool) :
    log (signedBase base i + signedBase base j + signedBase base k) =
      signedBase (fun n => log (base n)) i +
      signedBase (fun n => log (base n)) j +
      signedBase (fun n => log (base n)) k := by
  simp only [map_add, signedBase_map]

theorem signed_p192_accepted_support_ratio : 6 * 224 = 1344 := by norm_num

theorem signed_p192_hit_integer_certificate :
    1344 * 2 ^ 179 <
      (1569275433846670190958947353987675545133315652035359680560 : ℕ) := by
  norm_num

/-- Reversing the old degree-three divisor pair makes the same two-cubic
    fixed-origin reduction represent the negative class. The section-divisor
    equalities are the concrete finite-field obligations replayed by PARI. -/
theorem negative_fixed_origin_cubic_class_identity
    (A B T U D sectionG sectionF : G)
    (hG : sectionG = A + T + U)
    (hF : sectionF = B + U + D)
    (hsame : sectionG = sectionF) :
    D - T = -(B - A) := by
  calc
    D - T = A - B :=
      fixed_origin_cubic_class_identity B A T U D sectionG sectionF hG hF hsame
    _ = -(B - A) := by abel

end P192ExtensionCover

#print axioms P192ExtensionCover.signedBase_map
#print axioms P192ExtensionCover.signed_three_log_equation
#print axioms P192ExtensionCover.signed_p192_hit_integer_certificate
#print axioms P192ExtensionCover.negative_fixed_origin_cubic_class_identity
