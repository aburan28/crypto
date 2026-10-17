import RhoFormal.LeastRotation

/-!
# Bridges from the list model of `least_rotation` to the Rust bit operations

* `val_rot_eq_rust`: `rot` on MSB-first lists is the Rust
  `((v << k) | (v >> (n - k))) & mask` on the numeric value, for `k ≤ n`.
* `zmask`: the zero-run masks `Z_1 = ¬c`, `Z_{j+1} = Z_j ∧ rotate(¬c, j)`
  of the Rust loop, with `zmask_spec` (bit set iff the rotation starting there
  has `≥ j` leading zeros) and `rust_candidates_exact` (the loop stops at
  `len = L` and the set bits of `run` are exactly the candidates `IsCand`).
-/

namespace RhoFormal.LeastRotation

/-! ## Numeric form of `rot` -/

theorem val_append : ∀ (u v : List Bool), val (u ++ v) = val u * 2 ^ v.length + val v
  | [], v => by simp [val]
  | b :: u, v => by
    have ih := val_append u v
    cases b <;> simp [val, ih, Nat.pow_add] <;> grind

/-- **Bridge to `NormalBasis::rotate`.**  For `k ≤ n` (the Rust calls it with
`k < n`, and returns `v` itself for `k = 0`, which the formula also gives),
the value of `rot k w` is `((v <<< k) ||| (v >>> (n - k))) &&& (2^n - 1)`. -/
theorem val_rot_eq_rust (w : List Bool) (k : Nat) (hk : k ≤ w.length) :
    val (rot k w) =
      ((val w <<< k) ||| (val w >>> (w.length - k))) &&& (2 ^ w.length - 1) := by
  have hsplit : val w = val (w.take k) * 2 ^ (w.length - k) + val (w.drop k) := by
    conv => lhs; rw [← List.take_append_drop k w]
    rw [val_append]; simp
  have hT := val_lt (w.take k)
  have hD := val_lt (w.drop k)
  simp only [List.length_take, List.length_drop] at hT hD
  have hmin : min k w.length = k := by omega
  rw [hmin] at hT
  have hn : w.length = k + (w.length - k) := by omega
  -- v >>> (n - k) = T
  have hshr : val w >>> (w.length - k) = val (w.take k) := by
    rw [Nat.shiftRight_eq_div_pow, hsplit, Nat.add_comm, Nat.add_mul_div_right _ _ (by
      exact Nat.two_pow_pos _), Nat.div_eq_of_lt hD, Nat.zero_add]
  -- (v <<< k) % 2^n = D * 2^k
  have hshl : (val w <<< k) % 2 ^ w.length = val (w.drop k) * 2 ^ k := by
    rw [Nat.shiftLeft_eq, hsplit]
    have e : (val (w.take k) * 2 ^ (w.length - k) + val (w.drop k)) * 2 ^ k
        = 2 ^ w.length * val (w.take k) + val (w.drop k) * 2 ^ k := by
      conv => rhs; rw [hn, Nat.pow_add]
      grind
    rw [e, Nat.mul_add_mod]
    apply Nat.mod_eq_of_lt
    have : val (w.drop k) * 2 ^ k < 2 ^ (w.length - k) * 2 ^ k :=
      Nat.mul_lt_mul_of_pos_right hD (Nat.two_pow_pos _)
    rw [← Nat.pow_add] at this
    have e2 : w.length - k + k = w.length := by omega
    rw [e2] at this; exact this
  have hTn : val (w.take k) < 2 ^ w.length :=
    Nat.lt_of_lt_of_le hT (Nat.pow_le_pow_right (by decide) hk)
  rw [Nat.and_two_pow_sub_one_eq_mod, Nat.or_mod_two_pow, hshl, hshr,
    Nat.mod_eq_of_lt hTn, Nat.mul_comm, ← Nat.two_pow_add_eq_or_of_lt hT]
  simp only [rot, val_append, List.length_take, hmin]
  grind

/-! ## Indexing lemmas -/

theorem getElem_rot (w : List Bool) (k i : Nat) (hk : k ≤ w.length) (hi : i < w.length) :
    (rot k w)[i]'(by simp; exact hi) = w[(k + i) % w.length]'(Nat.mod_lt _ (by omega)) := by
  simp only [rot, List.getElem_append, List.length_drop, List.getElem_drop, List.getElem_take]
  split
  · next h =>
    have : (k + i) % w.length = k + i := Nat.mod_eq_of_lt (by omega)
    simp only [this]
  · next h =>
    have : (k + i) % w.length = i - (w.length - k) := by
      have e : k + i = (i - (w.length - k)) + w.length := by omega
      rw [e, Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
    simp only [this]

theorem le_lz_iff : ∀ (u : List Bool) (j : Nat) (hj : j ≤ u.length),
    (j ≤ lz u ↔ ∀ t (h : t < j), u[t]'(by omega) = false)
  | [], j, hj => by simp at hj; subst hj; simp [lz]
  | b :: u, 0, _ => by simp
  | b :: u, j + 1, hj => by
    have ih := le_lz_iff u j (by simp at hj; omega)
    cases b with
    | true => simp [lz]; exact ⟨0, by omega, rfl⟩
    | false =>
      simp only [lz]
      constructor
      · intro h t ht
        cases t with
        | zero => rfl
        | succ t => simpa using (ih.mp (by omega)) t (by omega)
      · intro h
        have := ih.mpr (fun t ht => by simpa using h (t + 1) (by omega))
        omega

/-! ## The zero-run masks of the Rust loop -/

/-- `zmask c 0` is all ones; `zmask c (j+1) = zmask c j ∧ rotate(¬c, j)`.
So `zmask c 1 = ¬c` is the Rust `zeros`, and the Rust's
`next = run & rotate(zeros, len)` with `run = zmask c len` is
`zmask c (len + 1)`. -/
def zmask (c : List Bool) : Nat → List Bool
  | 0 => List.replicate c.length true
  | j + 1 => List.zipWith (· && ·) (zmask c j) (rot j (c.map (!·)))

@[simp] theorem length_zmask (c : List Bool) : ∀ j, (zmask c j).length = c.length
  | 0 => by simp [zmask]
  | j + 1 => by simp [zmask, length_zmask c j]

/-- Bit `i` (MSB-first index, i.e. Rust bit `n-1-i`, offset `k = i`) of
`zmask c j` is set iff `c` is zero at the cyclic positions `i, …, i+j-1`. -/
theorem zmask_cyclic (c : List Bool) : ∀ j, j ≤ c.length → ∀ i (hi : i < c.length),
    ((zmask c j)[i]'(by simp; exact hi) = true ↔
      ∀ t, t < j → c[(i + t) % c.length]'(Nat.mod_lt _ (by omega)) = false)
  | 0, _, i, hi => by simp [zmask]
  | j + 1, hj, i, hi => by
    have ih := zmask_cyclic c j (by omega) i hi
    have hr := getElem_rot (c.map (!·)) j i (by simp; omega) (by simp; exact hi)
    simp only [List.length_map] at hr
    simp only [Nat.add_comm j i] at hr
    simp only [zmask, List.getElem_zipWith, Bool.and_eq_true, hr, List.getElem_map,
      Bool.not_eq_eq_eq_not, Bool.not_true, ih]
    constructor
    · intro ⟨h1, h2⟩ t ht
      by_cases e : t = j
      · subst e; exact h2
      · exact h1 t (by omega)
    · intro h
      exact ⟨fun t ht => h t (by omega), h j (by omega)⟩

/-- **`zmask_spec`.** For `j ≤ n`: bit `i` of `zmask c j` is set iff
`rot i c` has at least `j` leading zeros. -/
theorem zmask_spec (c : List Bool) (j : Nat) (hj : j ≤ c.length) (i : Nat)
    (hi : i < c.length) :
    (zmask c j)[i]'(by simp; exact hi) = true ↔ j ≤ lz (rot i c) := by
  rw [zmask_cyclic c j hj i hi, le_lz_iff (rot i c) j (by simp; exact hj)]
  constructor
  · intro h t ht
    rw [getElem_rot c i t (by omega) (by omega)]; exact h t ht
  · intro h t ht
    have := h t ht
    rw [getElem_rot c i t (by omega) (by omega)] at this; exact this

theorem lz_rot_lt_of_true (c : List Bool) (i : Nat) (hi : i < c.length)
    (htrue : ∃ x ∈ c, x = true) : lz (rot i c) < c.length := by
  apply Classical.byContradiction
  intro hlt
  have hall := (le_lz_iff (rot i c) c.length (by simp)).mp (by have := lz_le_length (rot i c); simp at this; omega)
  obtain ⟨x, hx, hxt⟩ := htrue
  obtain ⟨idx, hidx, he⟩ := List.getElem_of_mem hx
  have ht : ∃ t, t < c.length ∧ (i + t) % c.length = idx := by
    by_cases h : i ≤ idx
    · exact ⟨idx - i, by omega, by rw [Nat.mod_eq_of_lt (by omega)]; omega⟩
    · refine ⟨idx + c.length - i, by omega, ?_⟩
      have e : i + (idx + c.length - i) = idx + c.length := by omega
      rw [e, Nat.add_mod_right, Nat.mod_eq_of_lt hidx]
  obtain ⟨t, htl, hte⟩ := ht
  have := hall t htl
  rw [getElem_rot c i t (by omega) (by omega)] at this
  simp only [hte] at this
  rw [he, hxt] at this
  exact Bool.noConfusion this

/-- **Rust candidate set is exact.**  For a word `c` that is neither all-zero
nor all-one (the Rust handles those first; see `rot_const`):
* `1 ≤ L < n` where `L = maxLz c`;
* for `1 ≤ j ≤ L` the mask `zmask c j` is nonzero, and `zmask c (L+1)` is
  zero, so the loop `while run & rotate(zeros, len) ≠ 0 { run = …; len += 1 }`
  started at `run = zeros, len = 1` stops with `len = L` and `run = zmask c L`;
* the set bits of `zmask c L` are exactly the candidates `IsCand c i`. -/
theorem rust_candidates_exact (c : List Bool) (hf : ∃ x ∈ c, x = false)
    (ht : ∃ x ∈ c, x = true) :
    1 ≤ maxLz c ∧ maxLz c < c.length ∧
    (∀ j, 1 ≤ j → j ≤ maxLz c → ∃ i, ∃ hi : i < c.length,
        (zmask c j)[i]'(by simp; exact hi) = true) ∧
    (∀ i (hi : i < c.length), (zmask c (maxLz c + 1))[i]'(by simp; exact hi) = false) ∧
    (∀ i (hi : i < c.length), (zmask c (maxLz c))[i]'(by simp; exact hi) = true ↔ IsCand c i) := by
  have hne : c ≠ [] := by intro h; subst h; simp at hf
  have hpos : 0 < c.length := by cases c <;> simp_all
  obtain ⟨k, hk, hkL⟩ := exists_eq_maxLzUpTo c c.length hpos
  have hLlt : maxLz c < c.length := by
    have := lz_rot_lt_of_true c k hk ht; simp only [maxLz] at *; omega
  have hL1 : 1 ≤ maxLz c := by
    obtain ⟨x, hx, hxf⟩ := hf
    obtain ⟨idx, hidx, he⟩ := List.getElem_of_mem hx
    have h1 : 1 ≤ lz (rot idx c) := by
      rw [← zmask_spec c 1 (by omega) idx hidx, zmask_cyclic c 1 (by omega) idx hidx]
      intro t ht
      have : t = 0 := by omega
      subst this
      simp only [Nat.add_zero, Nat.mod_eq_of_lt hidx]; rw [he, hxf]
    have := le_maxLz c idx hidx; omega
  refine ⟨hL1, hLlt, ?_, ?_, ?_⟩
  · intro j _ hj
    refine ⟨k, hk, ?_⟩
    rw [zmask_spec c j (by omega) k hk]; simp only [maxLz] at *; omega
  · intro i hi
    have hspec := zmask_spec c (maxLz c + 1) (by omega) i hi
    have := le_maxLz c i hi
    cases h : (zmask c (maxLz c + 1))[i]'(by simp; exact hi) with
    | false => rfl
    | true => rw [h] at hspec; have := hspec.mp rfl; omega
  · intro i hi
    rw [zmask_spec c (maxLz c) (by omega) i hi]
    have := le_maxLz c i hi
    simp only [IsCand]
    constructor
    · intro h; exact ⟨hi, by omega⟩
    · intro ⟨_, h⟩; omega

end RhoFormal.LeastRotation
