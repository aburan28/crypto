/-!
# Norm of an affine involution: the parity rule behind the extra weak line

When `g` is odd about `x0` in `F_q`, the affine symmetry `nu(s) = -s + 2*x0` is
defined over `F_q`, so its Frobenius conjugates are all `nu` and the norm
`N(nu) = phi^(n-1)(nu) o ... o phi(nu) o nu` is `nu^n`. The extra weak line
exists exactly when this norm is the identity (PROTOCOL.md, Lemma A).

This file proves, with no Mathlib, the parity rule on affine maps
`(z, c) : s ↦ z*s + c` over `Int`:

* `comp_involution`: `(-1, c) o (-1, c) = id`;
* `pw_even`, `pw_odd`: `nu^(2k) = id` and `nu^(2k+1) = nu`;
* `norm_id_iff_even`: `nu^n = id` iff `n` is even.

Not formalized: Lemma A itself (Hilbert 90 for `G_a ⋊ G_m` and the bijection
between weak lines and norm-one symmetries); the computation of `Aut(g)`.
-/

namespace InvolutionNorm

/-- Affine maps `s ↦ z*s + c`, composition `(a ∘ b)(s) = a(b(s))`. -/
def comp (a b : Int × Int) : Int × Int := (a.1 * b.1, a.1 * b.2 + a.2)

/-- `nu^n` by repeated composition. -/
def pw (a : Int × Int) : Nat → Int × Int
  | 0 => (1, 0)
  | n + 1 => comp a (pw a n)

theorem comp_involution (c : Int) : comp (-1, c) (-1, c) = (1, 0) := by
  simp only [comp, Prod.mk.injEq]
  omega

theorem pw_even (c : Int) (k : Nat) : pw (-1, c) (2 * k) = (1, 0) := by
  induction k with
  | zero => rfl
  | succ k ih =>
    have h : 2 * (k + 1) = (2 * k) + 1 + 1 := by omega
    rw [h]
    simp only [pw] at ih ⊢
    rw [ih]
    simp only [comp, Prod.mk.injEq]
    omega

theorem pw_odd (c : Int) (k : Nat) : pw (-1, c) (2 * k + 1) = (-1, c) := by
  simp only [pw]
  rw [pw_even]
  simp only [comp, Prod.mk.injEq]
  omega

theorem norm_id_iff_even (c : Int) (n : Nat) : pw (-1, c) n = (1, 0) ↔ n % 2 = 0 := by
  constructor
  · intro h
    match Nat.mod_two_eq_zero_or_one n with
    | Or.inl hn => exact hn
    | Or.inr hn =>
      have hk : n = 2 * (n / 2) + 1 := by omega
      rw [hk, pw_odd] at h
      simp only [Prod.mk.injEq] at h
      omega
  · intro h
    have hk : n = 2 * (n / 2) := by omega
    rw [hk]
    exact pw_even c (n / 2)

end InvolutionNorm

#print axioms InvolutionNorm.comp_involution
#print axioms InvolutionNorm.pw_even
#print axioms InvolutionNorm.pw_odd
#print axioms InvolutionNorm.norm_id_iff_even
