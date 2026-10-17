import Mathlib.Data.Sym.Card
import Mathlib.Tactic.NormNum

/-!
The unordered three-sum support of a base with B elements is the image
of the size-three multiset type. This proves the exact combinatorial
ceiling independently of any particular elliptic curve or query law.
-/

namespace P192ExtensionCover

noncomputable def threeSumSupport {B : ℕ} {G : Type*} [AddCommMonoid G]
    (q : Fin B → G) : Finset G := by
  classical
  exact (Finset.univ : Finset (Sym (Fin B) 3)).image
    (fun s : Sym (Fin B) 3 => ((s : Multiset (Fin B)).map q).sum)

theorem three_sum_support_card_le {B : ℕ} {G : Type*} [AddCommMonoid G]
    (q : Fin B → G) :
    (threeSumSupport q).card ≤ (B + 2).choose 3 := by
  classical
  unfold threeSumSupport
  calc
    (Finset.univ.image
      (fun s : Sym (Fin B) 3 => ((s : Multiset (Fin B)).map q).sum)).card
        ≤ Fintype.card (Sym (Fin B) 3) := by
          simpa using
            (Finset.card_image_le
              (s := (Finset.univ : Finset (Sym (Fin B) 3)))
              (f := fun s : Sym (Fin B) 3 =>
                ((s : Multiset (Fin B)).map q).sum))
    _ = (B + 2).choose 3 := by
      rw [Sym.card_sym_eq_choose]
      simp

theorem p192_base_141_three_sum_capacity :
    (141 + 2).choose 3 = 477191 := by
  norm_num [Nat.choose_eq_descFactorial_div_factorial,
    Nat.descFactorial_succ, Nat.factorial]

theorem p192_base_141_three_sum_bound {G : Type*} [AddCommMonoid G]
    (q : Fin 141 → G) :
    (threeSumSupport q).card ≤ 477191 := by
  simpa [p192_base_141_three_sum_capacity] using
    (three_sum_support_card_le q)

end P192ExtensionCover

#print axioms P192ExtensionCover.three_sum_support_card_le
#print axioms P192ExtensionCover.p192_base_141_three_sum_capacity
#print axioms P192ExtensionCover.p192_base_141_three_sum_bound
