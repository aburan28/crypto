set_option exponentiation.threshold 1024

/-!
# Break-even arithmetic (README, "Break-even")

The cost model per relation, for `m = 4` summands on `n = 131`:
`cost_ω(ℓ, D) = C(4ℓ, D)^ω · 2^{131 − 4ℓ}`. The first factor is the
Macaulay-matrix size to the linear-algebra exponent `ω`. The second is the
inverse yield (one more dimension per summand buys `m = 4` bits). These are
kernel-checked integer instances **of that model**. They say nothing about
whether the model is right. That is what the experiments measure.

* slope 1 (`D = ℓ + 4`, the measured `m = 4` fit) and slope 1/2: growing
  the base always loses, for `ω ∈ {2, 3}` and `ℓ ∈ [8, 32]`;
* slope 0 (`D = 8` fixed): growing the base always wins (`ω = 2`: `ℓ ≥ 8`;
  `ω = 3`: `ℓ ≥ 10`);
* slope 1/4: wins at `ω = 2` and loses at `ω = 3`. The break-even slope lies
  between, which is why the preregistered threshold is `σ ≤ 0.25`.
-/

namespace FactorBaseDesigns.BreakEven

/-- `lo · (lo+1) ⋯ (lo+len−1)`. -/
def rise (lo : Nat) : Nat → Nat
  | 0 => 1
  | k + 1 => rise lo k * (lo + k)

/-- Binomial coefficient `C(n, k)` as `n^{\underline k} / k!`. -/
def binom (n k : Nat) : Nat := if k ≤ n then rise (n - k + 1) k / rise 1 k else 0

theorem binom_sanity : binom 5 2 = 10 ∧ binom 10 3 = 120 ∧ binom 52 5 = 2598960 ∧ binom 7 9 = 0 := by
  decide

def cost (ω ℓ D : Nat) : Nat := binom (4 * ℓ) D ^ ω * 2 ^ (131 - 4 * ℓ)

theorem slope_one_loses :
    (List.range 25).all (fun i => let ℓ := i + 8;
      decide (cost 2 ℓ (ℓ + 4) < cost 2 (ℓ + 1) (ℓ + 5)) &&
      decide (cost 3 ℓ (ℓ + 4) < cost 3 (ℓ + 1) (ℓ + 5))) = true := by decide +kernel

theorem slope_half_loses :
    (List.range 12).all (fun i => let ℓ := 2 * i + 8;
      decide (cost 2 ℓ (4 + ℓ / 2) ≤ cost 2 (ℓ + 2) (4 + (ℓ + 2) / 2)) &&
      decide (cost 3 ℓ (4 + ℓ / 2) ≤ cost 3 (ℓ + 2) (4 + (ℓ + 2) / 2))) = true := by
  decide +kernel

theorem slope_zero_wins :
    (List.range 25).all (fun i => let ℓ := i + 8; decide (cost 2 (ℓ + 1) 8 < cost 2 ℓ 8)) = true ∧
    (List.range 23).all (fun i => let ℓ := i + 10; decide (cost 3 (ℓ + 1) 8 < cost 3 ℓ 8)) = true := by
  decide +kernel

theorem slope_quarter_straddles :
    (List.range 6).all (fun i => let ℓ := 4 * i + 8;
      decide (cost 2 (ℓ + 4) (4 + (ℓ + 4) / 4) < cost 2 ℓ (4 + ℓ / 4)) &&
      decide (cost 3 ℓ (4 + ℓ / 4) ≤ cost 3 (ℓ + 4) (4 + (ℓ + 4) / 4))) = true := by
  decide +kernel

end FactorBaseDesigns.BreakEven
