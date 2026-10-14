# Protocol: root-finding decomposition oracles on both tracks, measured

Registered 2026-10-12, before any timed run of the grid below.  Builds and
correctness runs (verdict and root-set digests, brute-force checks) may
precede it; every timed row is produced after this file is committed.
Companion: the requirements document `research/pdp_root_finder_requirements_20261012/`
(gates G0–G7, §6) and the PS program's baseline `research/pdp_oracle_speedup_20261011/baseline_v1/`.

## Question

How much cheaper, per target and in counted field operations, is a
root-finding decomposition oracle (group law plus set membership) than the
Gröbner engines on the same instances, on the binary PS-0 rungs and on the
prime Kummer-tower cells the PKM rounds measured with F4; how does its cost
grow with the base dimension (`ℓ`, `t`) against the enumeration and
pair-table floors; and does any arm move the exponent?

## Arms

Binary (`examples/pdp_rootfind`, the PS harness's rung, corpus and digests):

| arm | what | memory |
|:--|:--|:--|
| `baseline` | PS-0's `semaev_decomp::decompose` (pairs-and-solve, quartic + `gcd(q, L_V mod q)`) | `O(2^ℓ)` |
| `group-scalar` | PS-2's scalar group-law row kernel | `O(2^ℓ)` |
| `rf-pairs` | this study's pairs-and-solve: both row points from one batched inversion, branch-free candidate test | `O(2^ℓ)` |
| `rf-mitm` | this study's 2-sum table: `x(P_i ± P_j)` sorted with a presence filter; `|F|` batched additions and `2|F|` lookups a target | `≈ 8 B` per stored sum (`|F|²`) |
| Gröbner reference | `pdp_split_dreg` unsplit `S₄`, `g = 0` (the split-dreg engine) on the PS-0 rung `(n, ℓ) = (24, 8)`, budget 1,800 s per target, 2 targets; censored if it does not finish | — |

Prime (`examples/pdp_rootfind`, Kummer cells rebuilt from the committed PKM
rows' `p, a, b, g, ζ, t, x_r`, and a fresh ladder):

| arm | what | memory |
|:--|:--|:--|
| `pairs3` | memoryless `m = 3`: row points `R ± P₁`, candidates `x(U ± P₂)` tested in the abscissa set of `V` | `O(2^t)` |
| `mitm3` | the 2-sum table, `2|F|` lookups a target | `≈ 8 B` per stored sum |
| `subtract2`, `mitm4` | `m = 2` and `m = 4` through the same primitives | as above |
| `f4_fp_tower` (reference) | the committed rows' multiply-adds, width, `D`, time (rounds 2–7); not re-run | — |
| exhaustive reference | `brute_rational` (general group law over every signed triple) and `algebraic_count` (the PKM checker's resultant count) | — |

## Units

- Counted: field multiplications (`mul`, the `3(k−1)` of a batched inversion
  included), squarings (`sqr`), inversions (`inv`, one per batched inversion),
  table lookups, candidate abscissae tested.  Deterministic; the metric.
- The F4 reference's unit is `F_p` multiply-adds (one multiplication and one
  addition each); a multiply-add and a field multiplication are compared one
  to one, which flatters F4 (its additions are free in that comparison).
- Wall: ns per target on a reserved core through `tools/isolated_bench.py run`,
  A/A (the binary against a copy of itself) then ABAB, at least five rounds,
  median and minimum reported; contended runs reported and not pooled.
- Memory: bytes of the table and filter.

## Cells and targets

- Binary: the PS-0 corpus, `ℓ = 8…12` (`n = 3ℓ`), the first 16 refuted
  uniform targets of each rung for timing (as PS-0), every target of the
  corpus for the digests; `ℓ = 13` and `14` (`n = 39, 42`) added with fresh
  corpora from the same generator and seed rule (`20261011 + ℓ`), 16 uniform
  + 4 planted targets, digests produced by the corpus's exhaustive enumerator.
- Prime, committed cells: every Kummer row of rounds 2, 3 and 7 at `p₁ =
  2013265921` (`m = 2`, `N = 12…22`; `m = 3`, `N = 9…18`; `m = 4`,
  `N = 8…16`), both targets each.
- Prime, ladder: `m = 3` at `p₁`, `t = 4…13`, seed `20261012`, 16 uniform +
  8 planted targets per `t`, the first 8 refuted uniform targets timed;
  `m = 2` and `m = 4` at `t = 4…10` for the `subtract2` / `mitm4` rows.

## Exactness gates (G0)

- Binary: verdict digest and root-set digest equal to the corpus's on every
  rung for every arm; a differing root set voids the arm's rows.
- Prime: every planted target found; every arm's root set equal to
  `brute_rational` where it runs (`t ≤ 7`); verdicts of all arms equal at
  every `t`; the committed row's `inconsistent` flag equal to
  `algebraic_count == 0` where the count runs (`t ≤ 8`, `m ≤ 3`); any
  algebraic solution with no rational decomposition reported.

## Predictions, before the runs

P1. `rf-mitm` decides a PS-0 target in fewer counted multiplications than
    `group-scalar` by a factor within 2× of `|F|/16` (the `2|F|`
    candidates against `≈ 2|F|²` of the pair loop, at 5–6 multiplications
    per base point against 10 per pair), i.e. ≈ 8× at `ℓ = 8` and ≈ 128× at
    `ℓ = 12`, and its wall ratio to PS-0's baseline exceeds 50× at `ℓ = 8`
    and 1,000× at `ℓ = 12`.
P2. `rf-pairs` is within 0.7–1.3× of `group-scalar` in wall: the row
    inversion sharing is a few per cent of the pair loop.
P3. Per-rung growth: `rf-mitm` ≈ 2 bits per rung of `ℓ` in counted work
    (`|F|` doubles, `2|F|` lookups), `rf-pairs` ≈ 2 bits: the pair table's
    slope at `m = 3` is 1 bit per rung and the memoryless floor 2 bits; no
    arm goes below the pair table's 1 bit per rung (no G3 pass is expected).
P4. Prime, committed cells: at `m = 3, N = 18` (`t = 6`) `mitm3` costs under
    `10³` multiplications and `pairs3` under `10⁴` against F4's `8.99·10¹³`
    multiply-adds; at `m = 4, N = 16` `mitm4` under `10⁴` against `2.7·10¹³`;
    the F4 multiply-add ÷ root-finder multiplication ratio is above `10⁹` on
    every `m ≥ 3` row with `N ≥ 15`.
P5. Prime ladder slopes: `pairs3` ≈ 2 bits per unit of `t` (`0.67` per tower
    variable), `mitm3` ≈ 1 bit per unit of `t` (`0.33` per tower variable):
    the enumeration and pair-table floors exactly; no arm below the pair
    table.
P6. Exactness: zero mismatches on every corpus and cell; the committed
    rows' `inconsistent` flags agree with the algebraic count; the `m = 3`
    algebraic count exceeds the rational decompositions on some uniform
    targets (twist-point solutions) and never on planted ones.

## Decision rule

- G0 failed on an arm → its rows are void until the cause is an instrument
  defect that is fixed and the run repeated whole.
- G1/G2 read from the one table: ratio columns against `T_enum`, the pair
  table, and the best Gröbner arm on the same cell; slopes fitted over the
  rungs each arm finishes (≥ 4 sizes required for a slope claim; two sizes are
  a ratio).
- G3 is passed only if a fitted slope with its 95 % bootstrap interval lies
  below the pair table's `⌊m/2⌋/m` per tower variable at polynomial memory;
  P3/P5 predict no pass.
- Class per AGENTS.md §3: every row expected to be **engineering** (a constant
  on the generic line); an **advance** only on a G3 pass.

## Inadmissible

Changing the base per target; counting a `BudgetExceeded` Gröbner run as a
refutation; choosing targets after seeing results; pooling contended wall
runs; quoting a per-call cost as a speed; comparing a table arm's time to a
memoryless arm's without its table in the row.
