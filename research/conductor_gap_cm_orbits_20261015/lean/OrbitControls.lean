import Std

-- This proves only the arithmetic bound after the norm has been identified.
-- CM classification and norm identification are separate mathematical bridges.
theorem norm_lower_bound (a2 norm ell : Nat)
    (h : norm = a2 + 19 * ell * ell) : 19 * ell * ell ≤ norm := by
  rw [h]
  exact Nat.le_add_left _ _

-- A cheap implementation and a scalar implementation of the same action
-- produce identical iterates. The actual point equality is checked natively.
def orbit {P : Type} (f : P → P) (x : P) : Nat → P
  | 0 => x
  | n+1 => f (orbit f x n)

theorem same_iterates {P : Type} (fast slow : P → P)
    (h : ∀ x, fast x = slow x) (x : P) :
    ∀ n : Nat, orbit fast x n = orbit slow x n := by
  intro n
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [orbit]
    rw [ih, h]

#print axioms norm_lower_bound
#print axioms same_iterates

-- This is the actual two-state recurrence invariant, conditional on the
-- one-step algebraic identity. No elliptic-curve classification is assumed.
def pairOrbit {P : Type} (step : P → P → P) (x y : P) : Nat → P × P
  | 0 => (x,y)
  | n+1 => let z := pairOrbit step x y n; (z.2, step z.2 z.1)

theorem recurrence_matches_orbit {P : Type} (f : P → P)
    (step : P → P → P) (h : ∀ z, step (f z) z = f (f z)) (x : P) :
    ∀ n, pairOrbit step x (f x) n = (orbit f x n, f (orbit f x n)) := by
  intro n
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [pairOrbit, ih, orbit]
    rw [h]

#print axioms recurrence_matches_orbit
