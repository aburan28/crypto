import Std.Tactic

/-! Exact 166-bit root pairs and a bounded fixed-slot N83 index. -/

namespace N83RootPair

def stateBound (k : Nat) : Nat := 83 * (k * (k - 1) / 2) + 42 * k
def pack (root0 root1 : Nat) : Nat := root0 + 2 ^ 83 * root1

theorem root_pair_roundtrip (root0 root1 : Nat)
    (h0 : root0 < 2 ^ 83) :
    pack root0 root1 % (2 ^ 83) = root0 ∧
    pack root0 root1 / (2 ^ 83) = root1 := by
  unfold pack
  constructor <;> omega

theorem root_pair_fits_21_bytes (root0 root1 : Nat)
    (h0 : root0 < 2 ^ 83) (h1 : root1 < 2 ^ 83) :
    pack root0 root1 < 2 ^ (21 * 8) := by
  unfold pack
  omega

theorem k3240_root_keys_fit_fixed_2pow30_slots :
    2 * stateBound 3240 < 2 ^ 30 := by decide

theorem k3240_state_id_fits_29_bits :
    stateBound 3240 + 1 < 2 ^ 29 := by decide

theorem k3240_payload_below_13point6_gibibytes :
    10 * (21 * stateBound 3240 + 5 * (2 ^ 30)) < 136 * (2 ^ 30) := by decide

#print axioms root_pair_roundtrip
#print axioms root_pair_fits_21_bytes
#print axioms k3240_root_keys_fit_fixed_2pow30_slots
#print axioms k3240_state_id_fits_29_bits
#print axioms k3240_payload_below_13point6_gibibytes

end N83RootPair
