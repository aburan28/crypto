# Correction: GHS magic number and genus type (2026-10-13)

`src/cryptanalysis/ec_trapdoor.rs` defined the GHS magic number as the plain
orbit dimension `dim_F2 span{σ^i(√b)}` and decided the genus type by "orbit
length of `√b` equals `m`".  Both are wrong against Maurer–Menezes–Teske,
*Analysis of the GHS Weil descent attack on the ECDLP over characteristic two
finite fields of composite degree*, LMS J. Comput. Math. 5 (2002):

| | before | now |
|---|---|---|
| magic number | `deg Ord_b = dim span{σ^i(√b)}` | `dim span{(1, σ^i(√b))}` (MMT eq. (5)) `= deg Ord_b + c`, `c = 1` iff `(x − 1) ∤ Ord_b` (Theorem 5(i)) |
| `magic_number_full` | added `(σ − 1)σ^i(a)` as generators | equals `magic_number`; `a` does not enter `m` |
| type I (genus `2^(m−1) − 1`) | orbit length of `√b` `== m` | `(x − 1) ∤ Ord_b`, i.e. `(1, 0)` lies in the MMT span |
| genus at `m = 1` | `0` (type I by the orbit rule) | `1` (`E` itself; type II) |

`Ord_b ∈ F_2[x]` is the `σ`-annihilator of `√b`, `σ: x ↦ x^(2^l)`, `n = N/l`.
For `n` odd, type I is exactly `Tr_{K/k}(b) = 0`, `b ∉ k`; for `n` a power of
two, `x^n − 1 = (x − 1)^n` and every curve is type II with `m = deg Ord_b`.

## Why (derivation, checked by tests)

The Artin–Schreier functions of the GHS compositum are
`x + σ^i(a) + σ^i(√b)/x`.  Constants `σ^i(a) − a` have absolute trace zero, so
they lie in `℘(K) = {h² + h}` and drop out; a relation modulo `℘(K(x))` needs
`Σ c_i = 0` (the `x` term) and `Σ c_i σ^i(√b) = 0` (the `1/x` term).  The span
therefore has dimension `dim span{(1, σ^i(√b))}`.  The genus of an elementary
abelian 2-extension of `K(x)` is the sum of the genera of its quadratic
subfields; `w² + w = c·x + const + β/x` has genus 1 iff `c = 1` and `β ≠ 0`,
else 0.  So the genus is the number of `(1, β)` in the span with `β ≠ 0`:
`2^(m−1)`, less one exactly when `(1, 0)` is in the span.

Evidence (all in this PR):

- `ec_trapdoor::tests::mmt_theorem5_counts_n12` — `N = 12`, `n = 3`:
  `#{m = 1} = 15`, `#{m = 3} = 4080`, `m = 2` not attained; `n = 2`:
  `63` / `4032`.  These are the counts in the `ghs_class_coverage` test
  `mmt_theorem5_counts` (branch `claude/isogeny-walks-genus3-d4c9g1`).
- `ec_trapdoor::tests::plain_orbit_undercounts_on_relative_trace_zero` —
  over every `b ∈ F_{2^12}^*`, `n = 3`: MMT `=` plain `+ 1` exactly on the 255
  `b ∉ k` with `Tr_{K/k}(b) = 0`, and those are exactly the type-I curves.
- `ec_trapdoor::tests::genus_type_n5` — `N = 10`, `n = 5`: generic `b` has
  genus 16, the 255 trace-zero `b ∉ k` genus 15.  This is the 15/16 split
  known from the `F_{2^155}`, `n = 5` GHS literature; the old rules gave 15 and
  8.
- `isogeny_degree_search::ball::tests::achievable_magic_numbers_match_exhaustive_mmt`
  — the divisor formula for attainable `m` equals an exhaustive scan over
  `F_{2^12}/F_2` and `F_{2^10}/F_2`.
- An independent pure-Python rank computation (not using the Rust code)
  reproduced the Oakley Group 3/4, sect113r1, c2pnb272w1, c2pnb304w1 and
  c2pnb368w1 values quoted below.

## Conclusions that change

Interpretation is labelled; every number is from the code in this PR.

### 1. Every curve over `F_2^131` (ECC2K-130's field)

| `b` | was `m` / genus | now `m` / type / genus |
|---|---|---|
| `b = 1` | 1 / 1 | 1 / II / 1 |
| `Tr(b) = 0`, `b ≠ 1` | **130 / `2^129`** | **131 / I / `2^130 − 1`** |
| `Tr(b) = 1` | 131 / `2^130` | 131 / II / `2^130` |

- The attainable magic numbers over `F_2^131` are `{1, 131}`, not
  `{0, 1, 130, 131}`.  `{0, 1, 130, 131}` is the set of Frobenius-invariant
  **subspace** dimensions; statements about invariant subspaces and
  quasi-subfield factor bases that cite it are unaffected.
- **Unchanged:** the GHS window `2 ≤ m ≤ 6` is empty over `F_2^131`; the
  `m = 1` transfer is the zero map on `⟨G⟩` (independent of the genus formula).
- **Shifted by one bit:** the smallest non-trivial GHS genus over `F_2` is
  `2^130 − 1`, not `2^129`.  In
  `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` the gap from
  the paying window `[130, 290…300]` becomes `130 − 8.23 = 121.77` bits (was
  120.77), and the §3.1 description-length floor becomes `≥ 2^130 − 1`
  operations, about `2^69.19×` the `2^60.81` rho reference (was `2^68.19×`).
  Interpretation: the GHS route for ECC2K-130 is closed by one more bit than
  stated; the verdict does not change.
- **§7 counterfactual (`N = 130`) changes and is not recomputed here.**  With
  MMT values, the least admissible `(m, genus)` for `l · g ≥ N − 2` is
  `l = 1`: `(9, 255)` (was `(8, 128)`); `l = 2`: `(13, 4095)` (was
  `(7, 64)`); `l = 13`: `(5, 15)` (was `(5, 16)`); `l = 26`: `(5, 15)` (was
  `(4, 8)`).  The row "genus 128, cost `2^36.41`, breaks: yes" therefore
  does not stand as written; the cost at genus 255 under the note's §5 model
  has not been computed, so whether the counterfactual still "breaks" is
  **open**.

Frozen artefacts that record the old values are not rewritten:
`experiments/ecc2k130_hyperelliptic_cover_boundary.json`
(`boundary_C_ghs_genus_trichotomy`), `experiments/isogeny_degree_search.json`
(`magic_numbers_131`), and
`research/notes/ecc2k130/curve_construction_20261006/results/curve_construction.json`
(`ghs.census`, `trichotomy_holds`).  Re-running the native producers gives
the values above; `curve_construction::legacy::ghs_census` now reports
`mmt_classes_hold` and keys `genus_by_class` by `m` and trace.

### 2. `docs/curves/jacobian-transfers.json` (regenerated)

`jacobian_transfer::binary::ghs_magic` already used the MMT magic number; only
its type rule was wrong.  Comparing a build of `main` (ecb69d15) with this
branch on the same registry and `--budget 200000`: 0 magic numbers change,
252 GHS rows change genus, 16 rows change class, and 4 curve verdicts change.
These are formula-level prices from the frozen cost model, not attack results.

| curve | best route (branch) | genus was → now | ratio_log2 was → now | verdict was → now |
|---|---|---|---|---|
| Oakley Group 3 (`N = 5·31`) | GHS to `GF(2^31)` | 15 → 16 | +1.15 → −3.04 | formula_marginal → boundary |
| Oakley Group 4 (`N = 5·37`) | GHS to `GF(2^37)` | 15 → 16 | +10.94 → +6.76 | formula_advantage → formula_marginal |
| c2pnb368w1 | (was GHS to `GF(2^92)`, `m = 4`) now the heuristic isogeny-extension row to `GF(2^46)` | 7 → 8 | +10.05 → +7.92 | formula_advantage → formula_marginal |
| c2pnb272w1 | GHS to `GF(2^2)`, `m = 8` "lower bound" row | 127 → 128 | (was unresolved at +41.32) → +72.20 | unresolved → formula_advantage |

Every `n = 2`, `m = 2` row moves from genus 1 / `not_carried` to genus 2 /
`boundary` (11 rows: c2pnb176w1, c2pnb208w1, c2pnb272w1, c2pnb304w1,
c2pnb368w1 and six toy-field registry curves).

Caveat on c2pnb272w1 (not changed here): the screen brackets `Tr(a) = 1`,
`n` even, as "`m` or `m + 1`" and marks both rows heuristic.  Under MMT Lemma 6
`a` does not change `m`; it decides applicability, condition (8): `n` odd or
`Tr_{K/k}(b) ≠ 0`.  For c2pnb272w1 and c2pnb304w1 at `l = 2`,
`Tr_{K/k}(b) ≠ 0` (independent check), so MMT gives `m = 8` exactly.  Replacing
the bracket by condition (8) is a separate modelling change.

The regenerated JSON and `docs/curves/JACOBIAN_TRANSFERS.md` table also pick
up drift that was already present on `main` (`jacobian_transfers check`
reported the committed file stale before this change); the four verdicts above
are the ones attributable to this correction.

### 3. sect113r1 (`N = 113`, prime)

Magic number 113 (unchanged), plain orbit rank 113, so `Ord_b = x^113 − 1` and
the curve is **type II, genus `2^112`**, not type I, `2^112 − 1`.  No row is
within the genus bound either way; the screen's conclusion is unchanged.  The
preserved, hash-pinned outputs under
`research/sect113r1-weak-curve-20261006/addenda/ghs-three-node-20261007/` and
`docs/ghs-transport-evidence/sect113r1-screen.json` keep the old type and
genus.  `tests/sect113r1_ghs_evidence.rs` now expects the corrected live
values; that test was already failing on `main` at its producer-source hash
pin (`src/cryptanalysis/ghs_screen.rs`), before any GHS value is compared, and
still fails there.

### 4. Other screens and code paths

- `docs/ghs-transport-evidence/binary192-screen.json` (`a = 0`, `b = 1`): all
  13 rows have `m = 1`; the recorded genus 0 / type I is genus 1 / type II.
  All 13 rows remain inside the genus-4 bound.  Frozen; not rewritten.
- `ghs_screen` / `audit_curve` with `a ∉ F_2`: `m` no longer depends on `a`
  (AES field, `b = 1`, `a = z`: was `m = 2, 3, 6`, now `1, 1, 1`).
- `ghs_descent::descend_m2_abstract`: the "type II, genus 2" 3-cycle case is
  MMT `m = 3`, type I, genus 3; the "type I, genus 1" 2-cycle case is `m = 2`,
  type II, genus 2.  The field `is_type_ii` is replaced by `three_cycle`,
  `magic`, `genus`.  `ghs_full_attack` dispatches both `m = 2` and `m = 3`.
- `binary_isogeny::walk_to_magic`, `ghs_full_attack`, `isogeny_degree_search`
  `screen` and the `ghs_attack_demo` example call `magic_number_full` and now
  get MMT values.  Over the `(N, n, l) = (6, 3, 2)` demo tower `m ∈ {1, 3}`;
  the demo's "m = 2" curve is `m = 3`, genus 3.
- `teske161_select` already used the MMT rank of pairs and the trace-zero ⇔
  genus 7 rule; unchanged.

## Not settled here

- MMT Lemma 6 applicability (`Tr_{K/F_2}(a) = 1`) is not checked by
  `audit_curve`.
- The `ghs_class_coverage` module on branch
  `claude/isogeny-walks-genus3-d4c9g1` asserts in `plain_matches_ec_trapdoor`
  that `ec_trapdoor::magic_number` equals its `magic_plain`; after this change
  that test must compare against `ec_trapdoor::plain_orbit_dim` instead.
- Conductor: `conductor check` could not run in the authoring session (no
  control-plane endpoint or login available in that container); no scope
  reservation was taken.
