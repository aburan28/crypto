/-!
# Least rotation of an n-bit word

Justifies `NormalBasis::least_rotation` in
`src/cryptanalysis/koblitz_strong_rho.rs`.

An `n`-bit word is a `List Bool` read most-significant bit first, so list
index `i` is Rust bit `n - 1 - i`.  Words are compared with core Lean's
lexicographic order `<` on `List Bool` (with `false < true`).
-/

namespace RhoFormal.LeastRotation

/-! ## Numeric value and the bridge to lexicographic order -/

/-- Numeric value of a word read most-significant bit first. -/
def val : List Bool → Nat
  | [] => 0
  | b :: bs => (if b then 2 ^ bs.length else 0) + val bs

theorem val_lt (w : List Bool) : val w < 2 ^ w.length := by
  induction w with
  | nil => simp [val]
  | cons b bs ih => cases b <;> simp [val, Nat.pow_succ] <;> omega

theorem not_lt_nil (u : List Bool) : ¬ u < [] := by
  intro h; cases u <;> cases h

/-- **(1a) Bridge lemma.** For words of equal length, the MSB-first numeric
value of `u` is below that of `v` iff `u` is lexicographically below `v`. -/
theorem val_lt_iff_lt : ∀ (u v : List Bool), u.length = v.length →
    (val u < val v ↔ u < v)
  | [], [], _ => by simp [val]
  | [], _ :: _, h => by simp at h
  | _ :: _, [], h => by simp at h
  | a :: u, b :: v, h => by
    have hl : u.length = v.length := by simpa using h
    have ih := val_lt_iff_lt u v hl
    have hu := val_lt u
    rw [hl] at hu
    have hft : (false < true) := by decide
    have htf : ¬ (true < false) := by decide
    have hv := val_lt v
    rw [List.cons_lt_cons_iff]
    cases a <;> cases b <;> simp [val, hl, hft, htf] at * <;> first | omega | exact ih

theorem val_inj : ∀ (u v : List Bool), u.length = v.length → val u = val v → u = v
  | [], [], _, _ => rfl
  | [], _ :: _, h, _ => by simp at h
  | _ :: _, [], h, _ => by simp at h
  | a :: u, b :: v, h, hv => by
    have hl : u.length = v.length := by simpa using h
    have h1 := val_lt u
    have h2 := val_lt v
    rw [hl] at h1
    cases a <;> cases b <;> simp [val, hl] at hv ⊢ <;>
      first | omega | exact val_inj u v hl (by omega)

theorem lt_asymm' {u v : List Bool} (hl : u.length = v.length) (h : u < v) : ¬ v < u := by
  rw [← val_lt_iff_lt _ _ hl] at h
  rw [← val_lt_iff_lt _ _ hl.symm]
  omega

theorem eq_of_not_lt_of_not_lt {u v : List Bool} (hl : u.length = v.length)
    (h1 : ¬ u < v) (h2 : ¬ v < u) : u = v := by
  rw [← val_lt_iff_lt _ _ hl] at h1
  rw [← val_lt_iff_lt _ _ hl.symm] at h2
  exact val_inj u v hl (by omega)

/-! ## Leading zeros -/

/-- Number of leading `false`s. -/
def lz : List Bool → Nat
  | [] => 0
  | false :: w => lz w + 1
  | true :: _ => 0

theorem lz_le_length : ∀ w : List Bool, lz w ≤ w.length
  | [] => by simp [lz]
  | false :: w => by simp [lz]; exact lz_le_length w
  | true :: _ => by simp [lz]

/-- **(1b)** Among words of equal length, more leading zeros is smaller. -/
theorem lt_of_lz_lt : ∀ (u v : List Bool), u.length = v.length → lz v < lz u → u < v
  | [], [], _, h => by simp [lz] at h
  | [], _ :: _, h, _ => by simp at h
  | _ :: _, [], h, _ => by simp at h
  | a :: u, b :: v, hl, h => by
    rw [List.cons_lt_cons_iff]
    cases a with
    | true => cases b <;> simp [lz] at h
    | false =>
      cases b with
      | true => left; decide
      | false =>
        right
        refine ⟨rfl, lt_of_lz_lt u v (by simpa using hl) ?_⟩
        simp [lz] at h; omega

/-! ## Rotations -/

/-- Left rotation by `k` (for `k ≤ n`), MSB first.  With list index `i` being
Rust bit `n-1-i`, this is `((v << k) | (v >> (n - k))) & mask`; see
`val_rot_eq_rust`. -/
def rot (k : Nat) (w : List Bool) : List Bool := w.drop k ++ w.take k

@[simp] theorem length_rot (k : Nat) (w : List Bool) : (rot k w).length = w.length := by
  simp [rot]; omega

/-- `k` is the offset of a lexicographically least rotation of `w`. -/
def IsLeastRot (w : List Bool) (k : Nat) : Prop :=
  k < w.length ∧ ∀ j, j < w.length → ¬ rot j w < rot k w

/-- **(1c) Main theorem.** Every rotation that is lexicographically least among
all rotations has the maximal number of leading zeros among all rotations,
i.e. it starts at the top of a longest cyclic run of zeros. -/
theorem least_rotation_has_max_lz (w : List Bool) (k : Nat) (hk : IsLeastRot w k) :
    ∀ j, j < w.length → lz (rot j w) ≤ lz (rot k w) := by
  intro j hj
  apply Classical.byContradiction
  intro h
  exact hk.2 j hj (lt_of_lz_lt _ _ (by simp) (by omega))

/-! ## Maximal leading-zero count and the candidate set -/

def maxLzUpTo (w : List Bool) : Nat → Nat
  | 0 => 0
  | m + 1 => max (maxLzUpTo w m) (lz (rot m w))

/-- `L`: the maximal `lz` over all `n` rotations. -/
def maxLz (w : List Bool) : Nat := maxLzUpTo w w.length

theorem le_maxLzUpTo (w : List Bool) : ∀ m k, k < m → lz (rot k w) ≤ maxLzUpTo w m
  | 0, _, h => by omega
  | m + 1, k, h => by
    simp only [maxLzUpTo]
    by_cases hk : k = m
    · subst hk; omega
    · have := le_maxLzUpTo w m k (by omega); omega

theorem exists_eq_maxLzUpTo (w : List Bool) :
    ∀ m, 0 < m → ∃ k, k < m ∧ lz (rot k w) = maxLzUpTo w m
  | 0, h => by omega
  | m + 1, _ => by
    simp only [maxLzUpTo]
    by_cases hm : m = 0
    · subst hm; exact ⟨0, by omega, by simp [maxLzUpTo]⟩
    · obtain ⟨k, hk, he⟩ := exists_eq_maxLzUpTo w m (by omega)
      by_cases hc : lz (rot m w) ≤ maxLzUpTo w m
      · exact ⟨k, by omega, by omega⟩
      · exact ⟨m, by omega, by omega⟩

theorem le_maxLz (w : List Bool) (k : Nat) (hk : k < w.length) : lz (rot k w) ≤ maxLz w :=
  le_maxLzUpTo w _ k hk

/-- Candidate offsets: `{k < n | lz (rot k w) = L}`. -/
def IsCand (w : List Bool) (k : Nat) : Prop := k < w.length ∧ lz (rot k w) = maxLz w

/-- **(1d-i)** The candidate set is nonempty for nonempty `w`. -/
theorem cand_nonempty (w : List Bool) (hw : w ≠ []) : ∃ k, IsCand w k := by
  have : 0 < w.length := by cases w <;> simp_all
  obtain ⟨k, hk, he⟩ := exists_eq_maxLzUpTo w w.length this
  exact ⟨k, hk, he⟩

/-- **(1d-ii)** Every least rotation (including every tie) is a candidate. -/
theorem least_is_cand (w : List Bool) (k : Nat) (hk : IsLeastRot w k) : IsCand w k := by
  refine ⟨hk.1, ?_⟩
  have h1 := le_maxLz w k hk.1
  have : w.length ≠ 0 := by have := hk.1; omega
  obtain ⟨j, hj, he⟩ := exists_eq_maxLzUpTo w w.length (by omega)
  have h2 := least_rotation_has_max_lz w k hk j hj
  simp only [maxLz] at *
  omega

/-- **(1d-iii)** Any `k` minimal within the candidates is minimal overall. -/
theorem cand_min_is_least (w : List Bool) (k : Nat) (hk : IsCand w k)
    (hmin : ∀ j, IsCand w j → ¬ rot j w < rot k w) : IsLeastRot w k := by
  refine ⟨hk.1, fun j hj => ?_⟩
  by_cases hc : lz (rot j w) = maxLz w
  · exact hmin j ⟨hj, hc⟩
  · have := le_maxLz w j hj
    have hlt : rot k w < rot j w := lt_of_lz_lt _ _ (by simp) (by rw [hk.2]; omega)
    exact lt_asymm' (by simp) hlt

theorem exists_argmin (f : Nat → Nat) (P : Nat → Prop) :
    ∀ n, (∃ k, k < n ∧ P k) → ∃ k, k < n ∧ P k ∧ ∀ j, j < n → P j → f k ≤ f j
  | 0, ⟨_, h, _⟩ => by omega
  | n + 1, ⟨k0, hk0, hp0⟩ => by
    by_cases hex : ∃ k, k < n ∧ P k
    · obtain ⟨k, hk, hp, hmin⟩ := exists_argmin f P n hex
      by_cases hn : P n ∧ f n < f k
      · refine ⟨n, by omega, hn.1, fun j hj hpj => ?_⟩
        by_cases hjn : j = n
        · subst hjn; omega
        · have := hmin j (by omega) hpj; omega
      · refine ⟨k, by omega, hp, fun j hj hpj => ?_⟩
        by_cases hjn : j = n
        · subst hjn
          apply Classical.byContradiction; intro hc; exact hn ⟨hpj, by omega⟩
        · exact hmin j (by omega) hpj
    · have hpn : P n := by
        by_cases hk : k0 = n
        · subst hk; exact hp0
        · exact absurd ⟨k0, by omega, hp0⟩ hex
      refine ⟨n, by omega, hpn, fun j hj hpj => ?_⟩
      by_cases hjn : j = n
      · subst hjn; omega
      · exact absurd ⟨j, by omega, hpj⟩ hex

/-- A rotation that is least within the candidate set exists. -/
theorem exists_cand_min (w : List Bool) (hw : w ≠ []) :
    ∃ k, IsCand w k ∧ ∀ j, IsCand w j → ¬ rot j w < rot k w := by
  obtain ⟨k0, hk0⟩ := cand_nonempty w hw
  obtain ⟨k, hk, hp, hmin⟩ :=
    exists_argmin (fun k => val (rot k w)) (IsCand w) w.length ⟨k0, hk0.1, hk0⟩
  refine ⟨k, hp, fun j hj hlt => ?_⟩
  rw [← val_lt_iff_lt _ _ (by simp)] at hlt
  have := hmin j hj.1 hj
  omega

/-- **(1d) Corollary justifying the algorithm.**  For nonempty `w`:
* the candidate set `{k < n | lz (rot k w) = L}` is nonempty;
* `k` is a least rotation overall iff it is a candidate and least among the
  candidates (so restricting the scan to candidates loses no minimiser and
  no tie);
* a least rotation exists, and all least rotations are the same word, so the
  minimum over the candidates equals the minimum over all rotations. -/
theorem least_rotation_via_candidates (w : List Bool) (hw : w ≠ []) :
    (∃ k, IsCand w k) ∧
    (∀ k, IsLeastRot w k ↔ (IsCand w k ∧ ∀ j, IsCand w j → ¬ rot j w < rot k w)) ∧
    (∃ k, IsLeastRot w k) ∧
    (∀ k k', IsLeastRot w k → IsLeastRot w k' → rot k w = rot k' w) := by
  refine ⟨cand_nonempty w hw, fun k => ⟨fun h => ⟨least_is_cand w k h, fun j hj => h.2 j hj.1⟩,
    fun ⟨hc, hm⟩ => cand_min_is_least w k hc hm⟩, ?_, ?_⟩
  · obtain ⟨k, hc, hm⟩ := exists_cand_min w hw
    exact ⟨k, cand_min_is_least w k hc hm⟩
  · intro k k' hk hk'
    exact eq_of_not_lt_of_not_lt (by simp) (hk'.2 k hk.1) (hk.2 k' hk'.1)

/-! ## Edge case: constant words -/

theorem rot_replicate (n k : Nat) (b : Bool) :
    rot k (List.replicate n b) = List.replicate n b := by
  simp only [rot, List.drop_replicate, List.take_replicate, List.replicate_append_replicate]
  congr 1; omega

/-- **(1e)** If every bit of `w` is the same (all `false` or all `true`),
every rotation equals `w`. -/
theorem rot_const (w : List Bool) (b : Bool) (hb : ∀ x ∈ w, x = b) (k : Nat) :
    rot k w = w := by
  have : w = List.replicate w.length b := List.eq_replicate_iff.mpr ⟨rfl, hb⟩
  rw [this, rot_replicate]

theorem rot_all_false (w : List Bool) (h : ∀ x ∈ w, x = false) (k : Nat) : rot k w = w :=
  rot_const w false h k

theorem rot_all_true (w : List Bool) (h : ∀ x ∈ w, x = true) (k : Nat) : rot k w = w :=
  rot_const w true h k

/-- For a constant word every offset `k < n` is a least rotation, so for
`n > 1` there is a tie (the Rust returns `tie = n > 1`). -/
theorem const_every_k_least (w : List Bool) (b : Bool) (hb : ∀ x ∈ w, x = b) (k : Nat)
    (hk : k < w.length) : IsLeastRot w k := by
  refine ⟨hk, fun j _ => ?_⟩
  rw [rot_const w b hb j, rot_const w b hb k]
  intro h; exact lt_asymm' rfl h h

end RhoFormal.LeastRotation
