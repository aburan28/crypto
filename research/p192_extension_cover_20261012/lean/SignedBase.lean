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

end P192ExtensionCover

#print axioms P192ExtensionCover.signedBase_map
#print axioms P192ExtensionCover.signed_three_log_equation
#print axioms P192ExtensionCover.signed_p192_hit_integer_certificate
