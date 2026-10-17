import Lean

/-!
The cleared chord resultant used by the different-j chi oracle.  This
module proves a ring identity conditional on the line, coordinate and
curve equations.  It does not formalize finite-field root completeness,
exceptional fibers or the Rust factor-base enumeration.
-/

namespace ChiChord

def eliminated (m d x₁ y₁ n a b : Int) : Int :=
  let l := m*x₁*d + n - y₁*d
  l*l*l + a*m*m*d*d*l + b*m*m*m*d*d*d - m*m*m*d*n*n

theorem eliminated_of_line_and_curve
    (m d x₁ y₁ x₂ y₂ n a b : Int)
    (line : m * (x₂ - x₁) = y₂ - y₁)
    (coordinate : n = d * y₂)
    (curve : y₂*y₂ = x₂*x₂*x₂ + a*x₂ + b) :
    eliminated m d x₁ y₁ n a b = 0 := by
  have hL : m*x₁*d + n - y₁*d = m*d*x₂ := by
    rw [coordinate]
    calc
      m*x₁*d + d*y₂ - y₁*d = d*(m*x₁ + (y₂-y₁)) := by
        simp [Int.mul_add, Int.mul_sub, Int.mul_comm, Int.mul_left_comm] <;> omega
      _ = d*(m*x₁ + m*(x₂-x₁)) := by rw [line]
      _ = m*d*x₂ := by
        simp [Int.mul_add, Int.mul_sub, Int.mul_comm, Int.mul_left_comm] <;> omega
  unfold eliminated
  dsimp
  rw [hL, coordinate]
  calc
    m*d*x₂*(m*d*x₂)*(m*d*x₂) + a*m*m*d*d*(m*d*x₂) +
        b*m*m*m*d*d*d - m*m*m*d*(d*y₂)*(d*y₂) =
        m*m*m*d*d*d*(x₂*x₂*x₂ + a*x₂ + b - y₂*y₂) := by
          simp [Int.mul_add, Int.mul_sub, Int.mul_comm, Int.mul_left_comm]
    _ = 0 := by rw [← curve]; simp

end ChiChord

#print axioms ChiChord.eliminated_of_line_and_curve
