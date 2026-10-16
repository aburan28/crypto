/-!
# The third conjugate cover and the descent cycle, without division

With `y0² = u·u''`, `y1² = u·u'` and `y2` defined by `u·y2 = y0·y1`:

* `third_cover`: `(u·y2)² = u²·(u'·u'')`, i.e. `y2² = u'·u''` once `u` is
  cancelled;
* `cycle_closes`: `u·(y1·y2) = (u·u')·y0`, i.e. `y1·y2/u' = y0` once `u` is
  cancelled, which is what makes `y0 ↦ y1 ↦ y2 ↦ y0` of order 3.

Only associativity and commutativity of multiplication are used, so both
hold in any commutative ring. Not formalized: cancelling `u` (needs a
domain and `u ≠ 0`), the function field, and the descent itself.
-/

namespace DescentCycle

structure CommMul (K : Type) where
  mul : K → K → K
  mul_assoc : ∀ a b c, mul (mul a b) c = mul a (mul b c)
  mul_comm : ∀ a b, mul a b = mul b a

variable {K : Type} (S : CommMul K)

local infixl:70 " ⬝ " => S.mul

theorem swap4 (a b c d : K) : (a ⬝ b) ⬝ (c ⬝ d) = (a ⬝ c) ⬝ (b ⬝ d) := by
  rw [S.mul_assoc a b (c ⬝ d), ← S.mul_assoc b c d, S.mul_comm b c,
    S.mul_assoc c b d, ← S.mul_assoc a c (b ⬝ d)]

/-- `(u·y2)² = (u·u)·(u'·u'')`. -/
theorem third_cover (u u1 u2 y0 y1 y2 : K)
    (h0 : y0 ⬝ y0 = u ⬝ u2) (h1 : y1 ⬝ y1 = u ⬝ u1) (h2 : u ⬝ y2 = y0 ⬝ y1) :
    (u ⬝ y2) ⬝ (u ⬝ y2) = (u ⬝ u) ⬝ (u1 ⬝ u2) := by
  rw [h2, swap4 S y0 y1 y0 y1, h0, h1, swap4 S u u2 u u1, S.mul_comm u2 u1]

/-- `u·(y1·y2) = (u·u')·y0`. -/
theorem cycle_closes (u u1 y0 y1 y2 : K)
    (h1 : y1 ⬝ y1 = u ⬝ u1) (h2 : u ⬝ y2 = y0 ⬝ y1) :
    u ⬝ (y1 ⬝ y2) = (u ⬝ u1) ⬝ y0 := by
  rw [S.mul_comm y1 y2, ← S.mul_assoc u y2 y1, h2, S.mul_assoc y0 y1 y1, h1,
    S.mul_comm y0 (u ⬝ u1)]

end DescentCycle

#print axioms DescentCycle.third_cover
#print axioms DescentCycle.cycle_closes
