/-!
Exact integer and prime-residue arithmetic certificates for the three accepted
different-j screen cells. These theorems do not define elliptic curves or prove
that the native Vélu code implements a group homomorphism.
-/

namespace OrderScreen

def JCertificate (p a b inverse result : Nat) : Prop :=
  let denominator := (4 * a ^ 3 + 27 * b ^ 2) % p
  let numerator := (1728 * 4 * a ^ 3) % p
  (denominator * inverse) % p = 1 ∧
    (numerator * inverse) % p = result

theorem destinationA :
    (-9 : Int) * (-21) - 30 * 3 ^ 2 = -81 := by decide

theorem destinationB :
    (-27 : Int) * 37 - 70 * 3 ^ 3 - 42 * (-21) * 3 = -243 := by decide

theorem p211SourceJ : JCertificate 211 190 37 112 77 := by
  unfold JCertificate
  decide
theorem p211DestinationJ : JCertificate 211 130 179 90 160 := by
  unfold JCertificate
  decide
theorem p223SourceJ : JCertificate 223 202 37 11 183 := by
  unfold JCertificate
  decide
theorem p223DestinationJ : JCertificate 223 142 203 216 222 := by
  unfold JCertificate
  decide
theorem p241SourceJ : JCertificate 241 220 37 119 33 := by
  unfold JCertificate
  decide
theorem p241DestinationJ : JCertificate 241 160 239 87 164 := by
  unfold JCertificate
  decide

theorem p211DifferentJ : (77 : Nat) % 211 ≠ 160 % 211 := by decide
theorem p223DifferentJ : (183 : Nat) % 223 ≠ 222 % 223 := by decide
theorem p241DifferentJ : (33 : Nat) % 241 ≠ 164 % 241 := by decide

theorem p211OrderFactor : (44655 : Nat) = 3 * 5 * 13 * 229 := by decide
theorem p223OrderFactor : (50175 : Nat) = 3 ^ 2 * 5 ^ 2 * 223 := by decide
theorem p241OrderFactor : (58563 : Nat) = 3 ^ 5 * 241 := by decide
theorem p211DistinctCapacity : (64 : Nat) ≤ 229 - 1 := by decide
theorem p211FiberCount : (14885 : Nat) * 3 = 44655 := by decide
theorem p211SubgroupPairCount : (229 : Nat) * 229 = 52441 := by decide
theorem p211KernelCoprime : (229 : Nat) % 3 = 1 := by decide

end OrderScreen

#print axioms OrderScreen.destinationA
#print axioms OrderScreen.destinationB
#print axioms OrderScreen.p211SourceJ
#print axioms OrderScreen.p211DestinationJ
#print axioms OrderScreen.p223SourceJ
#print axioms OrderScreen.p223DestinationJ
#print axioms OrderScreen.p241SourceJ
#print axioms OrderScreen.p241DestinationJ
#print axioms OrderScreen.p211DifferentJ
#print axioms OrderScreen.p223DifferentJ
#print axioms OrderScreen.p241DifferentJ
#print axioms OrderScreen.p211OrderFactor
#print axioms OrderScreen.p223OrderFactor
#print axioms OrderScreen.p241OrderFactor
#print axioms OrderScreen.p211DistinctCapacity
#print axioms OrderScreen.p211FiberCount
#print axioms OrderScreen.p211SubgroupPairCount
#print axioms OrderScreen.p211KernelCoprime
