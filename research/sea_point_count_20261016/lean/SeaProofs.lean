import Std

/- An explicit lawful additive abelian group. The fields below are group laws,
   not assumptions of BSGS correctness. No arithmetic or search claim is an axiom. -/
structure GroupOps (G : Type) where
  zero : G
  add : G → G → G
  neg : G → G
  assoc : ∀ a b c, add (add a b) c = add a (add b c)
  comm : ∀ a b, add a b = add b a
  zero_add : ∀ a, add zero a = a
  add_zero : ∀ a, add a zero = a
  add_neg : ∀ a, add a (neg a) = zero

namespace SeaProofs
variable {G : Type} (g : GroupOps G)

def smul (p : G) : Nat → G
  | 0 => g.zero
  | n+1 => g.add (smul p n) p

theorem smul_add (p : G) (a b : Nat) :
    smul g p (a+b) = g.add (smul g p a) (smul g p b) := by
  induction b with
  | zero => simp [smul, g.add_zero]
  | succ b ih =>
    calc
      smul g p (a+(b+1)) = g.add (smul g p (a+b)) p := rfl
      _ = g.add (g.add (smul g p a) (smul g p b)) p := congrArg (fun x => g.add x p) ih
      _ = g.add (smul g p a) (g.add (smul g p b) p) := g.assoc _ _ _

theorem subtract_left (a b q : G) :
    g.add a b = q ↔ b = g.add q (g.neg a) := by
  constructor
  · intro h
    calc
      b = g.add b g.zero := (g.add_zero b).symm
      _ = g.add b (g.add a (g.neg a)) := by rw [g.add_neg]
      _ = g.add (g.add b a) (g.neg a) := (g.assoc _ _ _).symm
      _ = g.add q (g.neg a) := by rw [g.comm b a, h]
  · intro h
    rw [h, g.comm q (g.neg a), ← g.assoc, g.add_neg, g.zero_add]

/- Exact predicate for the complete multimap implementation. A baby match at j
   in giant block i contributes k=i*width+j only when it is in the interval. -/
def Hits (p q : G) (count width k : Nat) : Prop :=
  ∃ i j, j < width ∧ k = i*width+j ∧ k < count ∧
    smul g p j = g.add q (g.neg (smul g p (i*width)))

theorem bsgs_complete (p q : G) (count width k : Nat) (hw : 0 < width) :
    Hits g p q count width k ↔ k < count ∧ smul g p k = q := by
  constructor
  · rintro ⟨i,j,_,rfl,hk,he⟩
    exact ⟨hk, (smul_add g p (i*width) j).trans ((subtract_left g _ _ _).mpr he)⟩
  · rintro ⟨hk,he⟩
    have hd : k = (k / width)*width + k % width := by
      simpa [Nat.mul_comm, Nat.add_comm] using (Nat.mod_add_div k width).symm
    refine ⟨k / width, k % width, Nat.mod_lt _ hw, hd, hk, ?_⟩
    apply (subtract_left g _ _ _).mp
    rw [← smul_add, ← hd, he]

/- Refinement is intersection with a sampled-point annihilator predicate.
   The sample hypothesis is explicit; deriving it for genuine curve points
   requires the group-order/Lagrange bridge outside this combinatorial proof. -/
def refine (candidates test : Nat → Prop) (k : Nat) : Prop := candidates k ∧ test k

theorem refinement_preserves (candidates test : Nat → Prop) (truth : Nat)
    (hc : candidates truth) (ht : test truth) : refine candidates test truth := ⟨hc,ht⟩

theorem direct_equals_full_intersection (candidates test : Nat → Prop) :
    refine candidates test = fun k => test k ∧ candidates k := by
  funext k
  apply propext
  exact And.comm

theorem unique_sound (candidates : Nat → Prop) (truth result : Nat)
    (ht : candidates truth) (hr : candidates result)
    (hu : ∀ a b, candidates a → candidates b → a = b) : result = truth :=
  hu result truth hr ht

/- Core of Hasse uniqueness: a genuine common divisor exceeding the interval
   width admits at most one multiple there. Exact group-order divisibility is
   an explicit premise; probable-prime factor discovery is not proved here. -/
theorem interval_unique (lower upper a b d : Nat)
    (ha : lower ≤ a ∧ a ≤ upper) (hb : lower ≤ b ∧ b ≤ upper)
    (hwidth : upper-lower < d) (hda : d ∣ a) (hdb : d ∣ b) : a = b := by
  obtain ⟨x, rfl⟩ := hda
  obtain ⟨y, rfl⟩ := hdb
  rcases Nat.lt_trichotomy x y with h | h | h
  · have hm := Nat.mul_le_mul_left d (Nat.succ_le_of_lt h)
    rw [Nat.mul_succ] at hm
    omega
  · rw [h]
  · have hm := Nat.mul_le_mul_left d (Nat.succ_le_of_lt h)
    rw [Nat.mul_succ] at hm
    omega

end SeaProofs

#print axioms SeaProofs.smul_add
#print axioms SeaProofs.subtract_left
#print axioms SeaProofs.bsgs_complete
#print axioms SeaProofs.refinement_preserves
#print axioms SeaProofs.direct_equals_full_intersection
#print axioms SeaProofs.unique_sound
#print axioms SeaProofs.interval_unique
