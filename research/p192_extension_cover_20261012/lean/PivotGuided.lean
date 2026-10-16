import SignedBase

/-!
The geometry filter selects scalars while retaining a known-scalar
three-sum equation. Concrete P-192 classes, cubic cuts, and support counts
are checked by the finite-field producers and independent replays.
-/

namespace P192ExtensionCover

variable {G K : Type*} [AddCommGroup G] [AddCommGroup K]

/-- A three-sum hit at the scalar-derived target `a * (-2 * generator)`
    gives an exact logarithm row. Conditioning the choice of `a` on any
    geometric predicate does not change this group identity. -/
theorem guided_scalar_three_sum_log_equation (log : G →+ K)
    (generator b₀ b₁ b₂ : G) (a : ℕ)
    (hit : b₀ + b₁ + b₂ = a • (-(generator + generator))) :
    log b₀ + log b₁ + log b₂ =
      a • (-(log generator + log generator)) := by
  calc
    log b₀ + log b₁ + log b₂ = log (b₀ + b₁ + b₂) := by simp
    _ = log (a • (-(generator + generator))) := congrArg log hit
    _ = a • (-(log generator + log generator)) := by simp

/-- The exact unordered-triple capacity of the 48-class measured base. -/
theorem two_pivot_48_class_three_sum_capacity :
    48 * 49 * 50 / 6 = 19600 := by
  norm_num

/-- With the previously certified accepted-scalar lower bound, the
    measured 4,926 admitted targets have ideal-uniform per-query hit
    probability below `2^-177`. The concrete cell count is replayed
    separately over the P-192 curve. -/
theorem two_pivot_guided_hit_bound_certificate :
    4926 * (2 ^ 177 : ℕ) <
      1569275433846670190958947353987675545133315652035359680560 := by
  norm_num

end P192ExtensionCover

#print axioms P192ExtensionCover.guided_scalar_three_sum_log_equation
#print axioms P192ExtensionCover.two_pivot_48_class_three_sum_capacity
#print axioms P192ExtensionCover.two_pivot_guided_hit_bound_certificate
