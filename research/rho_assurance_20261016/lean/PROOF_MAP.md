# Rho assurance: Lean proof map

Machine-checked lemmas behind the strong Pollard-rho reference
(`src/cryptanalysis/koblitz_strong_rho.rs`) and its Montgomery batch inversion.

- Toolchain: `leanprover/lean4:v4.34.1`, **core Lean only, no Mathlib**.
- Build and audit, one line: `./check.sh` (or `lake build && lake env lean RhoFormal/Axioms.lean`).
- Status: every statement below is **proved**. There is no `sorry`, `admit` or `native_decide`
  and no custom `axiom`. `axioms.log` shows that each reported theorem depends only on a
  subset of `propext`, `Classical.choice` and `Quot.sound`.
- Algebraic structure: core `Lean.Grind.Field` (part of Lean core in v4.34.1; it holds the
  usual commutative-ring axioms plus `inv_zero : 0⁻¹ = 0` and
  `mul_inv_cancel : a ≠ 0 → a * a⁻¹ = 1`). No class or axiom was defined here. The field
  is a typeclass parameter, so every result holds for every instance.
  `batch_inv` additionally takes `[DecidableEq F]`, which models the Rust `x != 0` test.

Files: `RhoFormal/LeastRotation.lean`, `RhoFormal/RotationBits.lean`,
`RhoFormal/BatchInv.lean`, `RhoFormal/RhoWalk.lean`, audit `RhoFormal/Axioms.lean`
(not in the default target), logs `build.log`, `axioms.log`.

## 1. Least rotation (`NormalBasis::least_rotation`, commit 23ed8d3a)

Model: a word is `List Bool`, most significant bit first. List index `i` is Rust bit
`n-1-i`, and offset `k` corresponds to the Rust `k = n-1-p` for top bit `p`. Words are
ordered by core Lean's lexicographic `<` on `List Bool` (`false < true`).
`val` is the MSB-first `Nat` value. `rot k w := w.drop k ++ w.take k`. `lz` counts
leading `false`s. `IsLeastRot w k := k < n ∧ ∀ j < n, ¬ rot j w < rot k w`.
`maxLz w` (= L) is the maximum of `lz (rot k w)` over `k < n`.
`IsCand w k := k < n ∧ lz (rot k w) = maxLz w`.

| Stmt | Lean name | Hypotheses | Status | Gap to Rust |
|---|---|---|---|---|
| 1a bridge | `val_lt_iff_lt` | `u.length = v.length`; concludes `val u < val v ↔ u < v` | proved | Rust compares `u64`/`u128` words. Words are masked to `n` bits, so `val < 2^n` holds automatically (`val_lt`). |
| (aux) | `val_inj`, `lt_asymm'`, `eq_of_not_lt_of_not_lt` | equal lengths | proved | Show the order on equal-length words is strict and total. |
| 1b | `lt_of_lz_lt` | `u.length = v.length`, `lz v < lz u`; concludes `u < v` | proved | none |
| 1c main | `least_rotation_has_max_lz` | `IsLeastRot w k`; concludes `∀ j < n, lz (rot j w) ≤ lz (rot k w)` | proved | "Top of a longest cyclic zero run" is expressed as "maximal `lz` over all rotations". `zmask_spec` (below) ties `lz (rot i c)` to the cyclic run at `i`. |
| 1d nonempty | `cand_nonempty` | `w ≠ []` | proved | none |
| 1d every least, and every tie, is a candidate | `least_is_cand` | `IsLeastRot w k` | proved | This is why the Rust `tie` flag is exact: no minimiser lies outside `tops`. |
| 1d min within candidates is min overall | `cand_min_is_least` | `IsCand w k`, `∀ j, IsCand w j → ¬ rot j w < rot k w` | proved | none |
| 1d existence | `exists_cand_min` | `w ≠ []` | proved | none |
| 1d packaged | `least_rotation_via_candidates` | `w ≠ []`; concludes: candidates nonempty; `IsLeastRot w k ↔ (IsCand w k ∧ least among candidates)`; a least rotation exists; all least rotations are the same word | proved | The Rust's choice of the smallest `k` among equal values is not formalized; it affects only which offset is reported, not the value or the tie flag. |
| 1e | `rot_const`, `rot_all_false`, `rot_all_true` | every bit equals `b` | proved | Matches the Rust early return `c == 0 \|\| c ^ mask == 0`. |
| 1e tie | `const_every_k_least` | constant `w`, `k < n` | proved | So for `n > 1` there is a tie, as the Rust's `tie = n > 1` says. |
| rotate bridge | `val_rot_eq_rust` | `k ≤ n`; concludes `val (rot k w) = ((v <<< k) \|\|\| (v >>> (n-k))) &&& (2^n-1)` with `v = val w` | proved | Stated on unbounded `Nat`. Rust `<<` on `u64`/`u128` drops bits above the word width, which the mask `2^n-1` (`n ≤` width) removes anyway; this step is argued, not proved. Rust returns `v` for `k = 0`, which is what the formula gives. |
| index form | `getElem_rot` | `k ≤ n`, `i < n`; concludes `(rot k w)[i] = w[(k+i) % n]` | proved | none |
| run-mask loop | `zmask_spec` | `zmask c 0 = 1…1`, `zmask c (j+1) = zmask c j ∧ rot j (¬c)`; for `j ≤ n`, `i < n`: bit `i` set ↔ `j ≤ lz (rot i c)` | proved | `zmask c 1 = ¬c` is the Rust `zeros`, and `zmask c (len+1)` is `run & rotate(zeros, len)`. AND on `u64` is modelled as `zipWith (&&)`. |
| candidate set exact | `rust_candidates_exact` | `c` has a `false` and a `true`; concludes `1 ≤ L < n`; `zmask c j ≠ 0` for `1 ≤ j ≤ L`; `zmask c (L+1) = 0`; bit `i` of `zmask c L` ↔ `IsCand c i` | proved | So the Rust loop stops at `len = L` with `run = zmask c L`, whose set bits are exactly the candidates. The loop control flow itself (the `loop`/`break`, and `trailing_zeros` iteration over `tops`) is not modelled as a program; the theorem states the facts it relies on. |

Not covered for (1): `to_normal`/`to_poly` (the normal-basis change of coordinates and
the claim that squaring is a one-bit rotation in it), and `choice_exhaustive`
(the `(x, y)` / negation tie-break). These still rely on the Rust tests
(`least_rotation_matches_the_exhaustive_scan`,
`canonical_form_is_orbit_invariant_and_multiplier_is_correct`).

## 2. Montgomery batch inversion

| Stmt | Lean name | Hypotheses | Status | Gap to Rust |
|---|---|---|---|---|
| model, 1 accumulator | `fwd`, `bwd`, `batchInv` | any `Field F`, `DecidableEq F` | definition | Mirrors `WideGf2::batch_inv` (`wide_gf2m.rs:216`, the field the rho uses for `n > 64`) and the identical `Gf2_128::batch_inv` (`semaev_decomp.rs:1116`). The forward pass pushes `acc` and then multiplies when `x ≠ 0`. There is one `inv` of the total. The reverse pass processes the highest index first and, when `x ≠ 0`, sets `out = inv·prefix` and `inv = inv·x`. |
| 2 | `batchInv_correct` | any list (including `[]` and lists with zeros); concludes the length is preserved and `∀ i < len, out[i] = if xs[i] = 0 then xs[i] else xs[i]⁻¹` | proved | Field elements are abstract. The GF(2^m) `mul`/`inv` implementations (Itoh–Tsujii, CLMUL folding) are assumed to be correct field operations and are not verified. |
| 2 (forms) | `batchInv_eq_map`, `batchInv_eq_map_inv` | none | proved | `out = xs.map (·⁻¹)`, using `0⁻¹ = 0`. |
| model, 4 accumulators | `fwdL`, `bwdL`, `combine4`, `batchInv4` | lane of index `i` is `i % 4` | definition | Mirrors `Gf2::batch_inv_inline` (and the test copy at `semaev_decomp.rs:3063`). `fwdL`/`bwdL` are proved for **any** lane assignment. `combine4` is exactly the `p01, p23, inv_all, i01, i23` combination. |
| 4-lane combine | `combine4_eq_inv` | every accumulator is nonzero (proved for the actual run by `fwdL_ne_zero`) | proved | none |
| 2 (4 lanes) | `batchInv4_correct`, `batchInv4_eq_map` | any list | proved | The Rust reverse loop visits chunks in descending order but the lanes *within* a chunk in ascending order. The model visits indices in strictly descending order. The four elements of a chunk use four different accumulators, so the two orders compute the same thing, but that reordering is **argued, not proved**. `batch_inv_clmul` keeps accumulators in a shifted vector representation (`fold_shifted`); its equivalence to `mul` is not modelled. |

## 3. Rho walk algebra

**Modelling assumption.** The prime-order subgroup `⟨G⟩` is identified with the scalar
field `F`: `[x]G ↦ x`, `G ↦ 1`, `Q ↦ d`, point addition ↦ `+`, and `[m]P ↦ m·P`. When
`F = Z/rZ` with `r` prime this is the group isomorphism `Z/rZ ≅ ⟨G⟩`. The following
are **not** formalized: that the binary-curve formulas in `StrongRho` (`add`, `double`,
`add_with_slope`) implement that group; that Frobenius^k and negation act on `⟨G⟩` as
`[±λ^k]` with the `lam_pow` table indexed consistently with the rotation direction; and
that `S::*_mod` arithmetic is the field `Z/rZ`. Scalars are a general field, not
`Z/rZ`. In particular, `inverse_mod` returning `None` exactly when the denominator is 0
relies on `r` being prime.

State `St = (P, a, b)`. `Rep d s := s.P = s.a + s.b * d` (Rust: `point = [a]G + [b]Q`).

| Stmt | Lean name | Hypotheses | Status | Gap to Rust |
|---|---|---|---|---|
| 3a jump | `jump_step` | `P = a + b·d`; concludes `P + c = (a+c) + b·d` | proved | A jump is `JumpG { point: [c]G, a: c }`, and the step keeps `b` (code: `a: add_mod(a, jump.a)`, `b: b`). |
| 3b canonicalize | `canon_step` | `P = a + b·d`; concludes `m·P = (m·a) + (m·b)·d` for **every** `m` | proved | Holds whichever `m` `canonicalize` picks (`±λ^k`, or `m = 1` for the point at infinity). |
| 3b compose | `canon_compose` | none | proved | `m₁(m₂x) = (m₁m₂)x` on point, `a` and `b`. |
| 3b negation | `canon_neg` | `P = a + b·d` | proved | `−P = (−1)·P = (−a) + (−b)·d`. |
| 3c escape | `escape_step` | `P = a + b·d`; concludes `P+P = (a+a) + (b+b)·d` | proved | Code: `add(next.point, next.point)`, `add_mod(a, a)`, `add_mod(b, b)`. |
| 3d collision | `collision` | `a + b·d = a' + b'·d`, `b ≠ b'`; concludes `d = (a'−a)·(b−b')⁻¹` | proved | Sign convention checked against `solve`: lane `(a, b) = lane.state`, hit `(a', b') = (hit_a, hit_b)`, code `candidate = (hit_a − a)·(b − hit_b)⁻¹`. They match. |
| 3d check never fails | `candidate_check` | as above | proved | `[candidate]G = Q` becomes `candidate·1 = d`. The step from this to "the Rust check never increments `failed_collisions` when the invariant holds" uses the modelling assumption above. |
| 3d wasted merge | `collision_same_b`, `denominator_zero_iff` | `b = b'`; concludes `a = a'`. Also `b − b' = 0 ↔ b = b'` | proved | Wasted, not wrong. |
| 3e cursor | `cursor_init`, `cursor_advance` | `CursorRep d C Ca := C = Ca + d` | proved | `cursor = [c]G + Q`, then `cursor += stride`, `cursor_a += stride_a`. |
| 3e one theorem | `rho_invariant` | `Reachable d`: the inductive closure of `start` (from a `CursorRep` cursor, `b = 1`), `step`, `canon` (any `m`) and `escape`; concludes (i) `Reachable → Rep`, (ii) each operation preserves `Rep`, (iii) two reachable states at the same point give `d = (a'−a)(b−b')⁻¹` if `b ≠ b'`, and `a = a'` if `b = b'` | proved | DP-table bookkeeping, `partition`, the fruitless-cycle detector, step caps and RNG are not modelled. They do not affect correctness. Skipping infinity starts is consistent with the model, since such a start is never reached. |

## Summary of what is not covered

1. Curve group ↔ scalar field (above), including the multiplier `±λ^k` of each automorphism.
2. Correctness of the GF(2^m) `mul`/`inv` code and the normal-basis maps.
3. Fixed-width machine words: `u64`/`u128` overflow in `rotate` is argued, not proved;
   lists stand in for words.
4. Imperative control flow: the run-mask loop, the `trailing_zeros` scan over `tops`, and
   the within-chunk order of the four-lane reverse pass are reduced to the facts they rely
   on, not modelled as programs.
5. Scalars form an abstract field rather than `Z/rZ` with `r` prime.
