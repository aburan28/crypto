import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Sum
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

/-!
Kernel-checked affine point-count certificates for the three prime-field
models in the leaderboard. The finite sum counts (x,y) solutions modulo p;
adding the unique point at infinity gives the elliptic-curve order when the
short Weierstrass discriminant is nonzero. The curve-group construction and
published attack reductions are separate formalization obligations.
-/

namespace WeakCurveLeaderboard

set_option maxRecDepth 16384

def affineCount (p a b : Nat) : Nat :=
  (Finset.range p).sum fun x =>
    ((Finset.range p).filter fun y =>
      (y * y) % p = (x * x * x + a * x + b) % p).card

set_option maxHeartbeats 10000000 in
theorem anomalous_F101_affine : affineCount 101 1 32 = 100 := by
  decide

set_option maxHeartbeats 10000000 in
theorem smooth_F101_affine : affineCount 101 1 0 = 99 := by
  decide

set_option maxHeartbeats 10000000 in
theorem pairing_F103_affine : affineCount 103 1 0 = 103 := by
  decide

theorem anomalous_trace_identity : 101 + 1 - (affineCount 101 1 32 + 1) = 1 := by
  rw [anomalous_F101_affine]

theorem pairing_embedding_degree_two : 103 % 13 ≠ 1 ∧ (103 * 103) % 13 = 1 := by
  decide

theorem anomalous_nonsingular : (4 * 1 ^ 3 + 27 * 32 ^ 2) % 101 ≠ 0 := by
  decide

theorem smooth_nonsingular : (4 * 1 ^ 3 + 27 * 0 ^ 2) % 101 ≠ 0 := by
  decide

theorem pairing_nonsingular : (4 * 1 ^ 3 + 27 * 0 ^ 2) % 103 ≠ 0 := by
  decide

theorem anomalous_point_on_curve : (10 ^ 2) % 101 = (4 ^ 3 + 4 + 32) % 101 := by
  decide

theorem smooth_point_on_curve : (38 ^ 2) % 101 = (3 ^ 3 + 3) % 101 := by
  decide

theorem pairing_point_on_curve : (44 ^ 2) % 103 = (18 ^ 3 + 18) % 103 := by
  decide

theorem pairing_selected_order_divides_group_order : 13 ∣ 104 := by
  decide

end WeakCurveLeaderboard
