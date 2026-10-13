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
- **§7 counterfactual (`N = 130`), recomputed** in
  [`research/ghs_counterfactual_n130_20261013/`](../../research/ghs_counterfactual_n130_20261013/README.md)
  (native `curve_construction::counterfactual`, the note's cost model and
  admissibility rule unchanged; checked against the legacy Python pricer).
  The least admissible `(m, genus)` for `l · g ≥ N − 2`:

  | `l` | was | now | cost was → now |
  |---:|---|---|---|
  | 1 | `(8, 128)` | `(9, type I, 255)` | `2^36.41` → **`2^55.60`** (rho `2^60.31`) |
  | 2 | `(7, 64)` | `(13, I, 4095)` | `2^299.93` → `2^43242.05` † |
  | 13 | `(5, 16)` | `(5, I, 15)` | `2^68.63` → `2^64.52` † |
  | 26 | `(4, 8)` | `(5, I, 15)` | `2^60.80` → `2^88.78` † |

  († fixed-genus model, pessimistic for large `g`.)  At `l = 1` the row still
  prices below rho in the note's model, but by **4.71 bits, not 23.91**.  The
  note's own §8 estimate of the cost unit (one relation trial is *perhaps*
  `2^7` curve additions) is larger than that margin; the earlier margin
  absorbed it.  Read that as: the "breaks: yes" entry survives the
  magic-number fix only inside the unconverted model.  A concrete `b` with
  `m = 9`, type I, genus 255 and condition (8) satisfied exists over `F_{2^130}`
  (constructed and checked); no curve with the subgroup in its Jacobian, and no
  transfer map, is constructed.
- **A second, independent error in the §7 table.**  Its `l ≥ 2` rows took the
  available magic numbers from the factorisation of `x^n − 1` over `F_{2^l}`
  (`achievable_dims(n, 1 << l)`); the magic number is an `F_2`-dimension, so the
  factorisation is over `F_2` at every `l`.  Exhaustive scans of `F_{2^12}` at
  `(n, l) = (6, 2), (4, 3), (3, 4), (2, 6)` find the `F_2` set (the `F_q` set
  would allow `m = 2` at `n = 3, l = 4`; no `b` has it).  This changes the `l = 2`
  and `l = 26` rows (and `l = 5, 10`, which the note's table omitted) in
  addition to the `+1` and the type; it does not touch `l = 1`.  The legacy
  script and its frozen JSON are not rewritten.

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
its type rule was wrong.  Comparing a build of `main` at ecb69d15 with this
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

### Lemma 6 replaces the `Tr(a)` bracket (same PR, later commit)

The screen used to bracket `Tr(a) = 1` with `n` even as "`m` or `m + 1`",
price both as `heuristic`, and mark disagreeing pairs `unresolved`.  Under MMT
Lemma 6, `a` does not change `m`; it decides whether the descent exists.
Condition (8): `n` odd, or `Tr(a) = 0`, or `Tr_{K/k}(b) ≠ 0`.

Why (8) is needed: if `n` is even and `Σ σ^i(√b) = 0`, the sum of all `n`
Artin–Schreier conjugates of `z² + z = x + a + b/x²` reduces to the constant
`Tr_{K/k}(a)`.  Its `℘`-root `h₀` lies in `F_{q²}`, and lies in `k = F_q`
exactly when `Tr_K(a) = 0`.  When `Tr_K(a) = 1`, `σ(h₀) = h₀ + 1`, so every
lift of `σ` to the compositum flips an odd number of the `z_i` and has order
`2n`.  There is then no Weil descent to `k`.

The screen now emits one `exact_invariant` row per tower with
`trace_a`, `relative_trace_b_zero` and `mmt_condition_8`.  A row where (8)
fails gets transport `descent_obstructed`, which counts as not carried.
`Tr(a)` is an isogeny invariant, so for such towers the isogeny-extension
rows count only targets satisfying (8).  Those have `(x − 1)^{2^e}` in the
annihilator and magic number exactly `d`.  The `unresolved` class and
`summary.unresolved` are removed (schema `jacobian-transfers/2`).

New test `mmt_condition_8_matches_the_relative_trace_of_b` runs exhaustively
over `F_{2^12}` on all six towers.  It compares the condition against
`Σ b^{q^i}` computed directly and checks the kernel size `q^{n−1} − 1`.  It
also checks that the magic number and type agree for `Tr(a) = 0` and `1`.

Effect on the regenerated file (2091 curves): two verdicts change, nothing
else changes class.

| curve | was | now |
|---|---|---|
| c2pnb176w1 | unresolved (`m = 8` lower-bound row, +23.64) | formula_advantage: same row, now exact, `m = 8`, genus 128 over `GF(4)`, +23.64 |
| c2pnb304w1 | unresolved (`m = 9` upper-bound row, +57.56) | formula_marginal: GHS to `GF(2^76)`, `m = 4`, genus 8, +2.00; its `GF(4)` row (`m = 8`, genus 128, 256 bits) cannot carry the 288-bit subgroup |

c2pnb272w1 keeps formula_advantage at +72.20; its `GF(4)` row is now exact
rather than a heuristic lower bound.  All three curves satisfy (8) on every
tower.  The isogeny-extension rows these three curves had are gone: each
target they priced would have been obstructed, and every row was
`not_carried` or `boundary`.  In the registry, (8) fails only on one 18-bit
curve (`icv1-f2m18-t999-40751283`, `b ∈ F_2`) at `l = 9, 3, 1`.  Those rows
were already `degenerate` / `not_carried`.

These are formula-level prices on exact invariants, not attacks: the cover
is not constructed and the map's action on the subgroup is not verified.
Graphs, the before/after table and a PDF are in
[`research/jacobian_transfer_lemma6_20261013/`](../../research/jacobian_transfer_lemma6_20261013/README.md).

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
- `ghs_class_coverage` (merged to `main` while this change was in
  progress): its `plain_matches_ec_trapdoor` test now compares `magic_plain`
  with `ec_trapdoor::plain_orbit_dim` and `magic_mmt` with
  `ec_trapdoor::magic_number`.  Its committed coverage data already used the
  MMT value and is unaffected.
- `teske161_select` already used the MMT rank of pairs and the trace-zero ⇔
  genus 7 rule; unchanged.

### Lemma 6 in `ec_trapdoor` and its callers (same PR, later commit)

- `ec_trapdoor::mmt_condition_8` (with `absolute_trace` and
  `relative_trace`) is the shared check.  `audit_curve` records it in
  `DescentRow::descends` and lists descending rows first.
- `construct_trapdoor_curve` keeps `a = 1`, which always satisfies (8):
  `Tr(1) = N mod 2`, and an even `n` forces an even `N`.  Its row records
  `descends`, and a debug assertion guards the invariant.
- `ghs_screen` reports `mmt_condition_8` per tower (schema
  `ghs-screen/v2`) and never counts an obstructed tower as within the
  bound.  `binary_isogeny::walk_to_magic` only stops on a curve whose
  cover descends.  `isogeny_degree_search::ball::ghs_hits` requires it
  (always true at `n = 131`, which is odd).  `ghs_full_attack` records it
  per neighbour, and `examples/ghs_attack_demo` prints it.
- New tests: `mmt_condition_8_counts_n12` (exhaustive over `F_{2^12}`, five
  towers: with `Tr(a) = 1` exactly `q^{n−1} − 1` failures for even `n`, none
  for odd `n`; none with `Tr(a) = 0`), `audit_curve_puts_descending_towers_first`
  and `trapdoor_curves_satisfy_condition_8`.  The `jacobian_transfer`
  exhaustive test also checks its own implementation against
  `ec_trapdoor::mmt_condition_8`.
- Effect on `jacobian-transfers.json`: only the 18-bit curve where (8)
  fails lists its routes in a different order.  No route, cost or verdict
  changes.
- The frozen `docs/ghs-transport-evidence/*-screen.json` and the sect113r1
  three-node addendum keep schema `ghs-screen/v1`.  sect113r1 has `N = 113`
  prime, so `n = 113` is odd and (8) holds; its live recompute is unchanged.

## Conductor

The authoring session for the MMT correction had no control plane, so no
scope reservation was taken.  The Lemma 6 follow-up started the in-container
control plane (`conductor up` on the image's Postgres 16, as
`.claude/hooks/session-start.sh` does).  `conductor check` for its scope was
clear, and it filed and claimed task T-1.  That control plane sees only
this container.
