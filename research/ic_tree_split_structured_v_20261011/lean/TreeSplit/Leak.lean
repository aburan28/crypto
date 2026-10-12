import TreeSplit.Replay

/-!
# Which `frob` subspaces sit inside the subfield `GF(2^7)`

Post-hoc (not registered). The registered `frob` arm excluded `g = t^d − 1`
for `d | n`, i.e. the subfields themselves, but not proper Frobenius-stable
subspaces **of** a subfield. `V = ker g(σ)` lies in `GF(2^7) = ker(σ^7 − 1)`
exactly when `g | t^7 − 1`. This file records which of the ten `g` the run
used do: `23`, `29` and `127` (`t^4+t^2+t+1`, `t^4+t^3+t^2+1`,
`t^6+…+1`). Those are the `n = 21` cells at `ℓ = 4` (both seeds) and
`ℓ = 6` (seed 20261012), and they are the only `frob` cells whose gap is
`−1`. None divides `t^3 − 1`.
-/

namespace TreeSplit.Leak
open TreeSplit.Replay

/-- `g` values of the registered `frob` cells at `n = 21`, in run order. -/
def gs : List Nat := [23, 35, 101, 159, 421, 29, 49, 127, 249, 331]

theorem in_gf2_7 : gs.filter (fun g => polyMod 129 g == 0) = [23, 29, 127] := by native_decide
theorem not_in_gf2_3 : gs.all (fun g => polyMod 9 g != 0) = true := by native_decide

end TreeSplit.Leak
