# Finding points versus finding relations: residual walks over partial decompositions

**Module:** `src/cryptanalysis/residual_walk.rs`
**Bench:**  `cargo run --release --example residual_walk_bench -- --panel --json experiments/20_residual_walk_panel.json`
**Data:**   `experiments/20_residual_walk_panel.json`, `experiments/20_residual_walk_panel.log`
**Tables:** `python3 scripts/summarize_residual_walk_panel.py experiments/20_residual_walk_panel.json`
**Baseline:** `experiments/20_residual_walk_baseline.json` (plain) and `experiments/20_residual_walk_tuned.json` (levers on), scored by `scripts/residual_walk_scoreboard.py` (§9)
**Round 3:** `experiments/20_residual_walk_seeded.json`, `experiments/20_residual_walk_structure.json` (§10)

> **Result in one line.**  Every residual-collision hybrid tested here
> recovers the planted logarithm correctly, and every one of them needs
> `≈ √(2n(B+1))` residuals to do it — a factor `√(2(B+1))` *more* than
> plain Pollard rho needs steps.  The walk changes how cheaply a residual
> is *constructed* (from `≈ 3·log₂n` group operations down to one or
> two), it does not change how many residuals a collision *needs*, and
> nothing in the construction makes residuals land in a smaller set for
> free.  No non-generic advantage was found, and the module documents
> exactly where one would have to come from.

## 1. The distinction

Let `Q = dG` in a subgroup of prime order `n`, and fix a factor base
`F = {P_0, …, P_{B−1}}` with unknown logarithms `ℓ_i = log_G P_i`.

*Finding points* is what a random walk is good at: run a rho-style walk
and every visited point is a candidate for anything you care to test —
"is it in `F`?", "is its `x` small?".  *Finding relations* is different:
a relation is an equation

```
  aG + bQ = Σ_i c_i P_i        ⟺        a + b·d ≡ Σ_i c_i ℓ_i   (mod n),
```

and `d` is recovered only once the collected equations span the unit
vector on `d` over `ℤ/nℤ`.  Enumerating many points with a property
does not produce such equations; independently chosen targets have to
*decompose* over the base, which is a property of the base relative to
random targets, not of the points in it.

The hybrid tested here keeps the walk but changes what it walks over.
A state is a partial decomposition

```
  s = (a, b, i_1, …, i_k),          L(s) = aG + bQ − Σ_j P_{i_j},
```

and the walk looks for two states with the same residual `L`.  Then

```
  (a − a')G + (b − b')Q = Σ_j P_{i_j} − Σ_j P_{i'_j}
```

is a relation entirely over `F`: the leftover cancels, exactly as in
the large-prime variation of index calculus.  A residual that lands in
`±F ∪ {O}` is a complete decomposition and is collected too.

## 2. What was built

One module, `cryptanalysis::residual_walk`, containing:

- `u64` affine arithmetic for prime-field curves with `p < 2^60`
  (cross-checked against the crate's `BigUint` point arithmetic in a
  unit test), a prime-order curve generator (random `p`, random `a, b`,
  BSGS over the Hasse interval, Miller–Rabin on the order), and a
  known-answer target `Q = dG`.
- A factor base of the `B` points with smallest `x`, indexed by `x` so
  a residual can be tested for membership in `±F` in `O(1)`.
- Five relation generators sharing one operation counter, one exact
  hash table keyed on the affine residual, one verifier, and one
  incremental reduced-row-echelon tracker over `ℤ/nℤ`:

  | tag | strategy                    | state                                  | cost of the next residual              | collision-preserving |
  |-----|-----------------------------|----------------------------------------|----------------------------------------|----------------------|
  | A   | independent partial sums    | fresh random `(a, b, k-tuple)`         | 2 scalar mults + `k` adds ≈ `3·log₂n`  | no                   |
  | B   | local-mutation walk         | swap one tuple slot; `L ← L + P_old − P_new`, sometimes `L ← L + G` or `L + Q` | 1–2 adds | no |
  | C1  | r-adding residual walk      | `L ← L + α_jG + β_jQ − P_j`, `j = h(L)` | 1 add (+ replay on collision)          | yes                  |
  | C2  | fresh-hash walk             | `s_{t+1} = H(L(s_t))`, residual recomputed | 2 scalar mults + `k` adds           | yes                  |
  | R   | plain Pollard rho           | `L ← L + α_jG + β_jQ`                  | 1 add                                  | yes (reference)      |

  *Collision-preserving* means equal residuals have equal successors,
  which is what makes distinguished-point storage legitimate.  B is the
  literal "change one point at a time and update the residual cheaply"
  experiment; C1 is the collision-preserving construction; C2 is the
  literal `s_{t+1} = H(L(s_t))` construction.  Four generic levers
  (negation map on the residual table, a `P_i − P_j` table for B,
  segmented walks for C1 and R, no restart after a collision for B)
  are available as options and measured in §9.7.

- Every collision is classified before it counts: **trivial** (same
  canonical state — tuples are stored as sorted multisets, so two
  orderings of the same summands are one state), **direct** (same
  multiset, different `(a, b)`: this is a plain rho collision and gives
  `d` at once), **factor-base-only** (`a = a'`, `b = b'`) or **mixed**.
  Every relation is re-verified by scalar multiplication before it is
  inserted, and the tracker reports whether it was independent of the
  rows already present.  A run stops at the first moment the row space
  contains the unit vector on `d` (or when the budget is spent), and
  the recovered value is scored against the planted one.
- The rejection filter of the "trap" argument (`filter_bound`: only
  residuals with `x < M` are looked up), distinguished points for the
  collision-preserving walks, and a meet-in-the-middle 4-decomposition
  search (`mitm_four_decomposition`: `P_i + P_j = R − P_k − P_l`).

Costs are counted in affine group operations (one field inversion
each).  Setup (multiplier tables), the walk itself, coefficient replay
after a collision, and verification are reported separately; linear
algebra is reported as wall time.

## 3. What the cost model predicts

Under the independent-uniform model for residuals:

1. **Birthday count.**  `T` stored residuals give `≈ T²/(2n)`
   collisions, so `R` relations need `T ≈ √(2nR)`.  Recovering `d`
   needs `R ≈ B + 1` independent relations, hence `√(2n(B+1))`
   residuals against rho's `√(πn/2)` steps: a deficit of
   `√(4(B+1)/π) ≈ 1.13·√(B+1)` in *count*, before the per-residual cost.
2. **Per-residual cost.**  A and C2 pay `≈ 3·log₂n` operations per
   residual; B and C1 pay one or two.  Total operations are the product.
3. **Rank cap.**  For an r-adding walk, every collision relation lies in
   the span of the `r` update vectors plus the `d` column, so
   `rank ≤ r + 1` (plus one per complete decomposition, an event of
   probability `2B/n` per step).  With `r < B` the walk is Pollard rho
   whose `r` multiplier logarithms happen to be unknown: it solves `d`
   after `r + 1` collisions and never learns the individual `ℓ_i`.
4. **The filter trap.**  Accepting only residuals in a set of size `M`
   costs `n/M` samples per accepted residual and `√(2MR)` accepted
   residuals, `n·√(2R/M)` in total — never below the unfiltered
   `√(2nR)`, and worse by `√(n/M)`.
5. **Meet in the middle.**  With a pair table of `B(B+1)/2` sums and
   `B(B+1)/2` streamed differences per target, expected matches per
   target are `(B(B+1)/2)²/n`, so `B ≈ (4n)^{1/4}` gives about one
   weight-4 relation per `≈ B²/2 ≈ 2√n` operations, and `B + 1` of them
   cost `≈ n^{3/4}`.

## 4. Measurements

All numbers below are produced by the panel command at the top of this
note, on random prime-order curves of the stated size, `k = 3`, seeds
`1..3`.  "correct = yes" means the recovered `d` equalled the planted
one in every run of the row; no run in the panel recovered a wrong
value or accepted an unverifiable relation.

### P1 — operations to recover d

Mean over seeds. `samples/pred` divides the residual count by `√(2n(B+1))` (A, B, C1, C2) or by `√(πn/2)` (R). `ops/rho` divides total group operations by `√(πn/2)`; `×R` divides them by plain rho's *measured* total on the same instances.

| bits | B | tag | seeds | samples | samples/pred | ops/sample | total ops | ops/rho | ×R | full | trivial | correct |
|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 20 | 64 | A | 3 | 9,652 | 1.02 | 56.3 | 544,433 | 523.3 | 109.4 | 1.3 | 0.0 | yes |
| 20 | 64 | B | 3 | 9,888 | 1.05 | 3.3 | 32,615 | 31.6 | 6.6 | 3.0 | 155.7 | yes |
| 20 | 64 | C1 | 3 | 8,962 | 0.95 | 8.1 | 71,345 | 69.7 | 14.3 | 2.0 | 0.0 | yes |
| 20 | 64 | C2 | 3 | 9,460 | 1.00 | 56.3 | 533,414 | 514.2 | 107.2 | 1.3 | 0.0 | yes |
| 20 | 64 | R | 3 | 1,227 | 1.15 | 4.5 | 4,974 | 4.8 | 1.0 | 0.0 | 0.0 | yes |
| 20 | 256 | A | 3 | 17,896 | 0.96 | 56.8 | 1,016,818 | 980.9 | 170.3 | 14.0 | 0.0 | yes |
| 20 | 256 | B | 3 | 18,799 | 1.00 | 3.4 | 64,243 | 62.2 | 10.8 | 15.3 | 69.3 | yes |
| 20 | 256 | C1 | 3 | 18,537 | 0.99 | 11.0 | 203,055 | 196.4 | 34.0 | 9.7 | 0.0 | yes |
| 20 | 256 | C2 | 3 | 18,321 | 0.98 | 56.7 | 1,040,439 | 1,006 | 174.3 | 12.3 | 0.0 | yes |
| 20 | 256 | R | 3 | 1,227 | 1.15 | 5.3 | 5,970 | 5.7 | 1.0 | 1.0 | 0.0 | yes |
| 24 | 64 | A | 3 | 40,128 | 1.00 | 68.6 | 2,758,713 | 624.7 | 219.8 | 0.3 | 0.0 | yes |
| 24 | 64 | B | 3 | 41,507 | 1.03 | 3.0 | 124,730 | 28.2 | 9.9 | 0.0 | 660.7 | yes |
| 24 | 64 | C1 | 3 | 41,604 | 1.04 | 5.0 | 208,364 | 47.2 | 16.6 | 0.0 | 0.0 | yes |
| 24 | 64 | C2 | 3 | 39,047 | 0.98 | 68.6 | 2,683,137 | 611.0 | 213.8 | 1.0 | 0.0 | yes |
| 24 | 64 | R | 3 | 4,282 | 1.01 | 3.0 | 12,549 | 2.9 | 1.0 | 0.0 | 0.0 | yes |
| 24 | 256 | A | 3 | 78,653 | 0.99 | 68.7 | 5,411,426 | 1,234 | 431.2 | 2.3 | 0.0 | yes |
| 24 | 256 | B | 3 | 81,380 | 1.03 | 2.4 | 198,054 | 45.2 | 15.8 | 2.7 | 313.3 | yes |
| 24 | 256 | C1 | 3 | 79,802 | 1.00 | 7.0 | 557,231 | 126.4 | 44.4 | 4.7 | 0.0 | yes |
| 24 | 256 | C2 | 3 | 77,936 | 0.98 | 68.7 | 5,363,578 | 1,220 | 427.4 | 3.7 | 0.0 | yes |
| 24 | 256 | R | 3 | 4,282 | 1.01 | 3.0 | 12,549 | 2.9 | 1.0 | 0.0 | 0.0 | yes |
| 28 | 64 | A | 3 | 158,757 | 0.98 | 80.6 | 12,809,251 | 716.1 | 229.1 | 0.0 | 0.0 | yes |
| 28 | 64 | B | 3 | 173,708 | 1.07 | 3.0 | 524,279 | 29.4 | 9.4 | 0.0 | 2,709 | yes |
| 28 | 64 | C1 | 3 | 162,247 | 1.00 | 4.7 | 749,340 | 42.2 | 13.4 | 0.0 | 0.0 | yes |
| 28 | 64 | C2 | 3 | 153,736 | 0.95 | 80.6 | 12,402,987 | 693.6 | 221.9 | 0.0 | 0.0 | yes |
| 28 | 64 | R | 3 | 19,169 | 1.06 | 3.0 | 55,904 | 3.1 | 1.0 | 0.0 | 0.0 | yes |
| 28 | 256 | A | 3 | 326,929 | 1.01 | 80.6 | 26,381,868 | 1,477 | 471.9 | 0.3 | 0.0 | yes |
| 28 | 256 | B | 3 | 331,344 | 1.03 | 2.2 | 719,765 | 40.4 | 12.9 | 0.7 | 1,272 | yes |
| 28 | 256 | C1 | 3 | 326,710 | 1.01 | 6.1 | 1,976,317 | 111.2 | 35.4 | 0.0 | 0.0 | yes |
| 28 | 256 | C2 | 3 | 312,522 | 0.97 | 80.7 | 25,221,491 | 1,411 | 451.2 | 0.3 | 0.0 | yes |
| 28 | 256 | R | 3 | 19,169 | 1.06 | 3.0 | 55,904 | 3.1 | 1.0 | 0.0 | 0.0 | yes |
| 28 | 1024 | A | 1 | 628,033 | 0.99 | 80.6 | 50,628,102 | 2,881 | 558.1 | 5.0 | 0.0 | yes |
| 28 | 1024 | B | 1 | 639,939 | 1.01 | 2.1 | 1,329,128 | 75.6 | 14.7 | 6.0 | 585.0 | yes |
| 28 | 1024 | C1 | 1 | 622,226 | 0.98 | 8.1 | 5,048,026 | 287.3 | 55.7 | 8.0 | 0.0 | yes |
| 28 | 1024 | C2 | 1 | 646,114 | 1.02 | 80.6 | 52,082,780 | 2,964 | 574.2 | 8.0 | 0.0 | yes |
| 28 | 1024 | R | 1 | 30,041 | 1.71 | 3.0 | 90,709 | 5.2 | 1.0 | 0.0 | 0.0 | yes |
| 32 | 256 | A | 1 | 1,174,703 | 1.02 | 91.6 | 107,622,613 | 1,683 | 801.0 | 0.0 | 0.0 | yes |
| 32 | 256 | B | 1 | 1,110,517 | 0.96 | 2.1 | 2,372,454 | 37.1 | 17.7 | 1.0 | 4,324 | yes |
| 32 | 256 | C1 | 1 | 1,160,619 | 1.00 | 5.0 | 5,773,396 | 90.3 | 43.0 | 0.0 | 0.0 | yes |
| 32 | 256 | C2 | 1 | 1,196,278 | 1.03 | 91.6 | 109,611,521 | 1,714 | 815.8 | 0.0 | 0.0 | yes |
| 32 | 256 | R | 1 | 45,242 | 0.71 | 3.0 | 134,363 | 2.1 | 1.0 | 0.0 | 0.0 | yes |
| 32 | 1024 | A | 1 | 2,278,244 | 0.99 | 91.6 | 208,792,157 | 3,266 | 1,554 | 5.0 | 0.0 | yes |
| 32 | 1024 | B | 1 | 2,375,757 | 1.03 | 1.9 | 4,549,954 | 71.2 | 33.9 | 1.0 | 2,295 | yes |
| 32 | 1024 | C1 | 1 | 2,299,084 | 1.00 | 7.3 | 16,884,291 | 264.1 | 125.7 | 2.0 | 0.0 | yes |
| 32 | 1024 | C2 | 1 | 2,281,725 | 0.99 | 91.6 | 209,111,870 | 3,271 | 1,556 | 1.0 | 0.0 | yes |
| 32 | 1024 | R | 1 | 45,242 | 0.71 | 3.0 | 134,363 | 2.1 | 1.0 | 0.0 | 0.0 | yes |

Breakdown of the operation count (means, same runs):

| bits | B | tag | setup | walk | replay | verify | LA ms | wall ms | table |
|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|
| 20 | 64 | A | 0.0 | 540,722 | 0.0 | 3,711 | 0.3 | 124.3 | 9,590 |
| 20 | 64 | B | 0.0 | 28,886 | 0.0 | 3,730 | 0.2 | 9.6 | 9,671 |
| 20 | 64 | C1 | 3,485 | 12,293 | 39,965 | 15,602 | 0.5 | 18.3 | 8,899 |
| 20 | 64 | C2 | 0.0 | 529,660 | 0.0 | 3,753 | 0.3 | 117.8 | 9,396 |
| 20 | 64 | R | 1,698 | 1,331 | 1,896 | 49.7 | 0.0 | 1.3 | 1,226 |
| 20 | 256 | A | 0.0 | 1,002,035 | 0.0 | 14,784 | 8.2 | 235.8 | 17,653 |
| 20 | 256 | B | 0.0 | 49,571 | 0.0 | 14,672 | 8.0 | 27.9 | 18,489 |
| 20 | 256 | C1 | 13,803 | 31,413 | 95,545 | 62,295 | 18.6 | 70.5 | 18,289 |
| 20 | 256 | C2 | 0.0 | 1,025,651 | 0.0 | 14,788 | 8.3 | 238.7 | 18,077 |
| 20 | 256 | R | 1,698 | 1,331 | 2,835 | 105.3 | 0.0 | 1.5 | 1,226 |
| 24 | 64 | A | 0.0 | 2,754,251 | 0.0 | 4,462 | 0.5 | 697.0 | 40,064 |
| 24 | 64 | B | 0.0 | 120,246 | 0.0 | 4,484 | 0.3 | 48.1 | 40,783 |
| 24 | 64 | C1 | 4,264 | 45,828 | 133,399 | 24,873 | 0.5 | 61.1 | 41,539 |
| 24 | 64 | C2 | 0.0 | 2,678,672 | 0.0 | 4,465 | 0.4 | 678.0 | 38,984 |
| 24 | 64 | R | 2,089 | 4,414 | 5,984 | 62.0 | 0.0 | 3.6 | 4,281 |
| 24 | 256 | A | 0.0 | 5,393,428 | 0.0 | 17,998 | 10.5 | 1,381 | 78,399 |
| 24 | 256 | B | 0.0 | 179,967 | 0.0 | 18,087 | 7.8 | 82.2 | 80,813 |
| 24 | 256 | C1 | 16,995 | 95,960 | 301,664 | 142,612 | 22.5 | 183.7 | 79,550 |
| 24 | 256 | C2 | 0.0 | 5,345,619 | 0.0 | 17,959 | 10.4 | 1,358 | 77,683 |
| 24 | 256 | R | 2,089 | 4,414 | 5,984 | 62.0 | 0.0 | 3.6 | 4,281 |
| 28 | 64 | A | 0.0 | 12,804,016 | 0.0 | 5,235 | 0.8 | 3,586 | 158,693 |
| 28 | 64 | B | 0.0 | 518,896 | 0.0 | 5,382 | 0.3 | 202.6 | 170,935 |
| 28 | 64 | C1 | 5,042 | 167,267 | 539,402 | 37,629 | 0.7 | 256.2 | 162,182 |
| 28 | 64 | C2 | 0.0 | 12,397,718 | 0.0 | 5,269 | 0.7 | 3,453 | 153,672 |
| 28 | 64 | R | 2,489 | 19,321 | 34,020 | 74.3 | 0.0 | 18.7 | 19,168 |
| 28 | 256 | A | 0.0 | 26,360,811 | 0.0 | 21,057 | 14.6 | 7,398 | 326,674 |
| 28 | 256 | B | 0.0 | 698,484 | 0.0 | 21,281 | 8.8 | 357.7 | 329,818 |
| 28 | 256 | C1 | 20,070 | 346,466 | 1,315,755 | 294,026 | 25.7 | 717.9 | 326,453 |
| 28 | 256 | C2 | 0.0 | 25,200,436 | 0.0 | 21,055 | 14.6 | 7,113 | 312,266 |
| 28 | 256 | R | 2,489 | 19,321 | 34,020 | 74.3 | 0.0 | 37.7 | 19,168 |
| 28 | 1024 | A | 0.0 | 50,543,752 | 0.0 | 84,350 | 586.8 | 14,725 | 627,015 |
| 28 | 1024 | B | 0.0 | 1,244,166 | 0.0 | 84,962 | 510.5 | 1,305 | 638,337 |
| 28 | 1024 | C1 | 80,295 | 700,259 | 2,699,634 | 1,567,838 | 1,391 | 3,245 | 621,209 |
| 28 | 1024 | C2 | 0.0 | 51,998,517 | 0.0 | 84,263 | 553.5 | 15,060 | 645,098 |
| 28 | 1024 | R | 2,509 | 30,190 | 57,940 | 70.0 | 0.0 | 108.1 | 30,040 |
| 32 | 256 | A | 0.0 | 107,598,848 | 0.0 | 23,765 | 20.4 | 32,941 | 1,174,448 |
| 32 | 256 | B | 0.0 | 2,347,248 | 0.0 | 25,206 | 10.3 | 1,494 | 1,105,940 |
| 32 | 256 | C1 | 22,961 | 1,183,067 | 4,125,987 | 441,381 | 32.1 | 2,535 | 1,160,362 |
| 32 | 256 | C2 | 0.0 | 109,587,594 | 0.0 | 23,927 | 21.1 | 34,149 | 1,196,022 |
| 32 | 256 | R | 2,835 | 45,394 | 86,042 | 92.0 | 0.0 | 249.2 | 45,241 |
| 32 | 1024 | A | 0.0 | 208,696,628 | 0.0 | 95,529 | 694.5 | 65,306 | 2,277,227 |
| 32 | 1024 | B | 0.0 | 4,450,947 | 0.0 | 99,007 | 539.9 | 3,948 | 2,372,439 |
| 32 | 1024 | C1 | 91,727 | 2,388,850 | 11,052,226 | 3,351,488 | 1,685 | 8,673 | 2,298,061 |
| 32 | 1024 | C2 | 0.0 | 209,015,876 | 0.0 | 95,994 | 696.5 | 65,039 | 2,280,702 |
| 32 | 1024 | R | 2,835 | 45,394 | 86,042 | 92.0 | 0.0 | 496.2 | 45,241 |

### P2 — relation yield at a fixed budget

n = 34850303 (2^25.1), B = 256, budget = 1,070,757 group operations, no early stop.

| tag | samples | accepted | independent | rank | solved | ops at solve | full | trivial |
|:--|---:|---:|---:|---:|:--|---:|---:|---:|
| A | 14,657 | 14,657 | 6 | 6 | no | - | 1 | 0 |
| B | 384,545 | 384,545 | 257 | 257 | yes | 320,758 | 6 | 1507 |
| C1 | 155,015 | 155,015 | 257 | 257 | yes | 838,693 | 0 | 0 |
| C2 | 14,658 | 14,658 | 6 | 6 | no | - | 1 | 0 |
| R | 188,838 | 188,838 | 4 | 4 | yes | 24,026 | 3 | 0 |

### P3 — rank of the r-adding residual walk vs r

n = 16343879, B = 256, budget 16,777,216 ops, no early stop.

| r | samples | verified | independent | dependent | rank | full hits | rank ≤ r+1+full | solved | ops at solve |
|---:|---:|---:|---:|---:|---:|---:|:--|:--|---:|
| 8 | 1,128,677 | 32965 | 49 | 32916 | 49 | 46 | yes | yes | 40,229 |
| 32 | 1,024,415 | 30393 | 55 | 30338 | 55 | 29 | yes | yes | 130,355 |
| 128 | 970,854 | 27975 | 145 | 27830 | 145 | 38 | yes | yes | 342,396 |
| 256 | 921,463 | 25722 | 257 | 25465 | 257 | 26 | yes | yes | 663,946 |

### P4 — rejecting residuals outside x < M (strategy A)

n = 12452773, p = 12454171, B = 64. Predicted cost ratio against the unfiltered run is `√(p/M)`.

| M | samples | accepted | accepted/pred | total ops | measured ratio | predicted ratio | correct |
|:--|---:|---:|---:|---:|---:|---:|:--|
| none (M = p) | 39,541 | 39,541 | 0.98 | 2,714,913 | 1.00 | 1.00 | yes |
| p / 2^2 | 78,222 | 19,804 | 0.98 | 5,368,539 | 1.98 | 2.00 | yes |
| p / 2^4 | 153,735 | 9,607 | 0.96 | 10,548,811 | 3.89 | 4.00 | yes |
| p / 2^6 | 322,362 | 5,098 | 1.01 | 22,112,430 | 8.14 | 8.00 | yes |
| p / 2^8 | 655,791 | 2,524 | 1.00 | 44,977,003 | 16.57 | 16.00 | yes |

### P5 — distinguished points on the collision-preserving walks

n = 147966617, B = 256.

| tag | dp bits | samples | stored | stored/samples | total ops | replay ops | walks | abandoned | correct |
|:--|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| C1 | 0 | 267,157 | 266,900 | 0.9990 | 1,867,945 | 1,264,238 | 258 | 0 | yes |
| C2 | 0 | 281,844 | 281,589 | 0.9991 | 22,345,108 | 0 | 255 | 0 | yes |
| R | 0 | 7,894 | 7,893 | 0.9999 | 26,112 | 15,574 | 2 | 0 | yes |
| C1 | 4 | 271,347 | 16,841 | 0.0621 | 1,880,515 | 1,272,618 | 258 | 0 | yes |
| C2 | 4 | 286,249 | 17,695 | 0.0618 | 113,182,424 | 90,488,810 | 257 | 0 | yes |
| R | 4 | 7,913 | 512 | 0.0647 | 26,169 | 15,612 | 2 | 0 | yes |
| C1 | 8 | 329,473 | 1,036 | 0.0031 | 2,140,029 | 1,463,132 | 258 | 0 | yes |
| C2 | 8 | 387,934 | 1,236 | 0.0032 | 147,121,607 | 116,374,520 | 305 | 0 | yes |
| R | 8 | 8,129 | 32 | 0.0039 | 26,817 | 16,044 | 2 | 0 | yes |

### P6 — meet in the middle on 4-decompositions

| bits | n | B | pair table | coincidences | targets | relations | independent | total ops | ops/relation | ops/√(πn/2) | pred. matches/target | correct |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 20 | 705883 | 41 | 861 | 0 | 209 | 212 | 40 | 211,693 | 5,292 | 201.0 | 1.05 | yes |
| 20 | 705883 | 82 | 3,396 | 7 | 29 | 426 | 80 | 126,624 | 1,583 | 120.3 | 16.41 | yes |
| 24 | 13306663 | 86 | 3,740 | 1 | 360 | 490 | 87 | 1,437,677 | 16,525 | 314.5 | 1.05 | yes |
| 24 | 13306663 | 172 | 14,872 | 6 | 60 | 972 | 171 | 979,902 | 5,730 | 214.3 | 16.63 | yes |
| 28 | 199310677 | 169 | 14,365 | 0 | 902 | 980 | 169 | 13,260,581 | 78,465 | 749.4 | 1.04 | yes |
| 28 | 199310677 | 338 | 57,284 | 7 | 127 | 1953 | 334 | 7,536,513 | 22,564 | 425.9 | 16.47 | yes |

## 5. Reading the numbers

**The count is the invariant.**  In every P1 row, A, B, C1 and C2 need
`samples/pred` between `0.95` and `1.07` — the same `√(2n(B+1))`
residuals, from `n ≈ 2^19` to `n ≈ 2^32` and `B` from 64 to 1024,
whether the residual was drawn independently, mutated locally, or
produced by a collision-preserving walk.  Quadrupling `B` doubles the
count, as `√(B+1)` says.  The walks do not find useful collisions
faster than independent sampling; they only find them *cheaper per
residual*.

**Cheaper per residual is the whole gain, and it is bounded.**  B and C1
bring the per-residual cost from `56–92` operations (A and C2: two
scalar multiplications each) down to `2–3` (B) and `5–11` (C1, of which
most is coefficient replay), which is why they finish `20–45×` and
`8–15×` cheaper than A at the same count.  That gain saturates: one
group operation per residual is the floor, plain rho already sits on it,
and rho needs `√(2(B+1))` times fewer residuals.  Measured against rho's
own total on the same curve (the `×R` column, which includes rho's setup
and replay), the best hybrid costs `7–34×` more (B) and `13–126×` more
(C1), growing with `B`; A and C2 cost `100–1,500×` more.  The
fixed-budget run (P2) shows the same thing from the other side: with a
budget of `8·√(2n(B+1))` operations at `n ≈ 2^25`, B and C1 solve `d`
(at `0.32M` and `0.84M` operations), A and C2 manage six relations each,
and rho is done after `24k`.

**Full decompositions are noise.**  Complete decompositions (residual
in `±F`) occur `2B/n` of the time per residual; at 28 bits and
`B = 256` that is a handful per run.  The large-prime route (leftover
collisions) supplies essentially all relations, which is the expected
regime whenever `√(2n(B+1)) ≪ n/B`.

**Trivial collisions are real and were excluded.**  The local-mutation
walk revisits its own states (swap `i → j` then `j → i`, or return to
a state with the same sorted multiset) at a steady rate; those are
counted in the `trivial` column and never enter the linear algebra.
The collision-preserving walks produce none, as they must.

**Rank cap (P3).**  With `r` multipliers the r-adding residual walk
never exceeds rank `r + 1` plus the number of complete decompositions,
however many relations it collects: at `r = 8` it verified `33k`
relations of which `49` were independent, and solved `d` after `40k`
operations — it is rho with nine unknowns, not index calculus.  A collision-preserving walk with cheap incremental
updates can only ever generate relations inside the span of its update
vectors; to reach the full base it needs `r ≥ B` distinct updates, at
which point the coefficient vectors are dense and must be replayed (or
stored) — the `replay` column in P1.

**The filter trap (P4).**  Rejecting residuals outside `x < p/2^s`
multiplies the total work by `≈ 2^{s/2}`: measured `1.98, 3.89, 8.14,
16.57` against predicted `2, 4, 8, 16` — while the accepted-residual
count stays on the `√(2M(B+1))` birthday line (`accepted/pred` of
`0.96–1.01`).  Filtering for special-looking residuals
after the fact buys nothing, exactly as the model says.

**Distinguished points (P5).**  For C1 and R, storing only residuals
with `dp` zero hash bits cuts the table by `≈ 2^{dp}` (`6.2%` and
`0.31%` of the residuals kept at `dp = 4, 8`) for `1–15%` more
operations; the relation count and the solve are unchanged.  C2 pays
much more (`5–6×`) because each replayed step costs two scalar
multiplications and the replay has to reach back to the merge point.  For A and B the module refuses the option, because a
residual collision between non-preserving states does not propagate
to a later distinguished point.  There is a second, less obvious
requirement that the first version of this experiment got wrong: the
fresh-hash walk is *memoryless*, so once two walks have merged they
carry identical states, and the collision observed at the next
distinguished point is trivial.  The relation lives only at the merge
point, where the two predecessor states differ, and has to be located
by replaying both walks from their stored starts (the van
Oorschot–Wiener step).  Without that replay C2 with `dp = 4` needed
`16×` the residuals and produced thousands of trivial collisions; with
it, its count returns to the `√(2n(B+1))` line and the replay cost
appears in the `replay ops` column.  C1 and R avoid the issue only
because their `(walk, step)` states keep their coefficient history.

**Meet in the middle (P6).**  At `B ≈ (4n)^{1/4}` the 4-decomposition
search delivers `1.0–1.1` weight-4 relations per target (predicted
`1.05`) at `≈ B²/2` operations each; `B + 1` such relations recover
`d`, at a total of `201×, 315×, 749×` rho's `√(πn/2)` for 20, 24 and
28 bits — growing like `n^{1/4}`, i.e. a total of order `n^{3/4}`.
Doubling `B` finds `16×` more matches per target at `4×` the cost per
target, so the per-relation cost falls, but the number of relations
needed rises with `B` and the total stays hundreds of times above rho.

## 6. Where a real gain would have to come from

Everything above is consistent with the generic-group bound: a
procedure that only adds, negates and compares group elements cannot
beat `√n` for prime-order DLP, and "collide residuals of partial
decompositions" is such a procedure however the states are arranged.
The experiments pin down the two places a non-generic idea would have
to act:

1. **Cheap construction of residuals in a small set.**  The count is
   `√(2M(B+1))` when residuals are confined to a set of size `M`.  The
   trap experiment shows that *rejecting* into such a set costs `n/M`
   per accepted residual.  A gain requires *constructing* residuals in
   a set of size `M ≪ n` at a cost that does not grow like `n/M` — which
   means using coordinates or algebra, not group operations.  The
   summation-polynomial decomposition oracles elsewhere in this
   repository (`ec_index_calculus`, `semaev_*`, `koblitz_*`) are exactly
   attempts at that; they replace the walk, they do not decorate it.
2. **Relations outside the update span.**  A collision-preserving walk
   only ever produces relations in the span of its update vectors (P3).
   Any walk-based relation generator therefore needs `≥ B` distinct
   update directions, and that is where the replay/storage cost of the
   coefficient vectors enters.

What the walk *does* give is a legitimate low-memory realisation of
large-prime-style relation collection (C1 and C2 with distinguished
points): correct, verifiable, and no better than `√(2n(B+1))`.

## 7. What this does not show

- Toy sizes (`n ≤ 2^32`) and a single curve family (random prime-order
  short-Weierstrass curves over prime fields).  The measured constants
  are specific to affine `u64` arithmetic with one inversion per
  operation; the *ratios* are what the argument rests on.
- `k = 3` throughout for the explicit-state strategies.  Larger `k`
  changes the weight of the relations, not the collision count.
- No negation map, no automorphism folding, no batched inversions in
  any strategy — the comparison is like-for-like, not tuned.
- The linear algebra is dense Gaussian elimination and is not the
  bottleneck at these sizes; its time is reported but not analysed.
- Nothing here bears on binary or extension-field curves, where index
  calculus with summation polynomials has genuine asymptotic content
  (see `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, `RESEARCH_SEMAEV_DECOMPOSITION.md`).

## 8. Reproduce

```bash
cargo test --release --lib residual_walk
cargo run --release --example residual_walk_bench -- --bits 24 --fb 256 --trials 3
cargo run --release --example residual_walk_bench -- --panel --quick
cargo run --release --example residual_walk_bench -- --panel --json experiments/20_residual_walk_panel.json
python3 scripts/summarize_residual_walk_panel.py experiments/20_residual_walk_panel.json
cargo run --release --example residual_walk_bench -- --baseline --json run.json
python3 scripts/residual_walk_scoreboard.py run.json --baseline experiments/20_residual_walk_baseline.json
cargo run --release --example residual_walk_bench -- --baseline --tuned --json tuned.json
python3 scripts/residual_walk_scoreboard.py tuned.json --baseline experiments/20_residual_walk_tuned.json
cargo run --release --example residual_walk_bench -- --seeded --json seeded.json      # §10.2
cargo run --release --example residual_walk_bench -- --structure --json structure.json  # §10.3
cargo run --release --example residual_walk_bench -- --bits 28 --j0 --aut --negation --diff-table --continue --strategies B,R
```

The bench's single-instance mode also accepts `--strategies A,C1,R`,
`--dp BITS` (applied to the collision-preserving strategies only),
`--k`, `--budget OPS` and `--json FILE`.

## 9. Optimisation ledger and baseline scoreboard

Sections 4–5 report what each strategy costs.  This section states
what each *step* bought, where the remaining generic slack is, and
freezes a baseline that later work is scored against.

### 9.1 Protocol

```bash
cargo run --release --example residual_walk_bench -- --baseline --json run.json      # ≈ 50 s
python3 scripts/residual_walk_scoreboard.py run.json                                  # score it
python3 scripts/residual_walk_scoreboard.py run.json \
        --baseline experiments/20_residual_walk_baseline.json --fail-on-regression    # compare
```

The protocol is fixed: `n ≈ 2^24` and `2^28`, `B = 256`, `k = 3`, seeds
1–3, all five strategies, plus C1 and R with 8 distinguished-point
bits.  It is deterministic for a given build, so two runs of the same
code produce identical cells; the frozen reference is
`experiments/20_residual_walk_baseline.json`.  The comparison prints an
improvement factor per cell and exits non-zero on a regression beyond
the tolerance (default 10%), so it can gate a change.

### 9.2 Metrics

For every `(bits, B, dp, strategy)` cell the score is decomposed as

```
  S = total_ops / √n  ≈  κ · c  +  overheads
  κ = samples / √n            the collision *count* factor
  c = walk_ops / samples      group operations per residual
  overheads                   setup, coefficient replay, verification (fractions of total_ops)
```

with the floors `κ_floor = √(2(B+1)) = 22.67` for any residual-collision
method at `B = 256` (birthday), `κ_rho = √(π/2) = 1.25` for plain rho,
and `S_floor = κ_floor · 1` (one operation per residual, no overhead).
`ops/rel/√n` is the cost of one *independent* relation in the same
units; `stored/√n` is peak memory.

### 9.3 The frozen baseline

| bits | B | dp | tag | seeds | κ = samples/√n | κ floor | κ/floor | c ops/residual | setup | replay | verify | S = ops/√n | S floor | S/floor | ops/rel/√n | stored/√n | trivial | correct |
|---:|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 24 | 256 | 0 | A | 3 | 22.53 | 22.67 | 0.99 | 68.5 | 0.0% | 0.0% | 0.3% | 1,547.0 | 22.67 | 68.24 | 6.051 | 22.459 | 0 | yes |
| 24 | 256 | 0 | B | 3 | 23.26 | 22.67 | 1.03 | 2.2 | 0.0% | 0.0% | 9.2% | 56.7 | 22.67 | 2.50 | 0.221 | 23.097 | 313 | yes |
| 24 | 256 | 0 | C1 | 3 | 22.69 | 22.67 | 1.00 | 1.2 | 3.1% | 53.9% | 25.7% | 158.4 | 22.67 | 6.99 | 0.616 | 22.622 | 0 | yes |
| 24 | 256 | 0 | C2 | 3 | 22.26 | 22.67 | 0.98 | 68.5 | 0.0% | 0.0% | 0.3% | 1,528.8 | 22.67 | 67.43 | 5.965 | 22.188 | 0 | yes |
| 24 | 256 | 0 | R | 3 | 1.27 | 1.25 | 1.01 | 1.0 | 16.9% | 47.6% | 0.5% | 3.7 | 1.25 | 2.95 | 3.692 | 1.270 | 0 | yes |
| 24 | 256 | 8 | C1 | 3 | 91.57 | 22.67 | 4.04 | 1.1 | 1.8% | 47.1% | 17.3% | 273.3 | 22.67 | 12.06 | 1.064 | 0.097 | 0 | yes |
| 24 | 256 | 8 | R | 3 | 1.34 | 1.25 | 1.07 | 1.0 | 15.9% | 48.5% | 0.5% | 3.9 | 1.25 | 3.11 | 3.903 | 0.004 | 0 | yes |
| 28 | 256 | 0 | A | 3 | 22.95 | 22.67 | 1.01 | 80.6 | 0.0% | 0.0% | 0.1% | 1,850.6 | 22.67 | 81.63 | 7.248 | 22.929 | 0 | yes |
| 28 | 256 | 0 | B | 3 | 23.31 | 22.67 | 1.03 | 2.1 | 0.0% | 0.0% | 3.0% | 50.6 | 22.67 | 2.23 | 0.199 | 23.203 | 1,272 | yes |
| 28 | 256 | 0 | C1 | 3 | 22.87 | 22.67 | 1.01 | 1.1 | 1.0% | 66.2% | 15.0% | 139.3 | 22.67 | 6.15 | 0.542 | 22.855 | 0 | yes |
| 28 | 256 | 0 | C2 | 3 | 21.93 | 22.67 | 0.97 | 80.6 | 0.0% | 0.0% | 0.1% | 1,768.8 | 22.67 | 78.02 | 6.910 | 21.912 | 0 | yes |
| 28 | 256 | 0 | R | 3 | 1.33 | 1.25 | 1.06 | 1.0 | 6.1% | 59.7% | 0.2% | 3.9 | 1.25 | 3.10 | 3.890 | 1.327 | 0 | yes |
| 28 | 256 | 8 | C1 | 3 | 37.16 | 22.67 | 1.64 | 1.0 | 0.9% | 62.8% | 13.0% | 165.0 | 22.67 | 7.28 | 0.642 | 0.089 | 0 | yes |
| 28 | 256 | 8 | R | 3 | 1.35 | 1.25 | 1.08 | 1.0 | 5.9% | 60.0% | 0.2% | 4.0 | 1.25 | 3.16 | 3.966 | 0.005 | 0 | yes |

Optimisation ledger (dp = 0 cells; factors are ratios of S, > 1 means cheaper):

| bits | B | A → B (mutation) | A → C1 (r-adding) | A → C2 (fresh hash) | B → R (drop the base) | C1 → R | B / S floor | C1 / S floor | R / rho floor |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 24 | 256 | 27.29× | 9.77× | 1.01× | 15.36× | 42.90× | 2.50× | 6.99× | 2.95× |
| 28 | 256 | 36.55× | 13.28× | 1.05× | 13.01× | 35.82× | 2.23× | 6.15× | 3.10× |

Distinguished points (memory bought per operation spent):

| bits | B | tag | dp | stored/√n dp=0 | stored/√n dp | memory ÷ | S dp=0 | S dp | ops × |
|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|
| 24 | 256 | C1 | 8 | 22.622 | 0.0967 | 234× | 158.4 | 273.3 | 1.73× |
| 24 | 256 | R | 8 | 1.270 | 0.0042 | 301× | 3.7 | 3.9 | 1.06× |
| 28 | 256 | C1 | 8 | 22.855 | 0.0890 | 257× | 139.3 | 165.0 | 1.18× |
| 28 | 256 | R | 8 | 1.327 | 0.0054 | 247× | 3.9 | 4.0 | 1.02× |

### 9.4 What each step bought (28-bit row, `B = 256`)

| step | change | κ | c | overheads | S | gain vs previous | distance to generic floor |
|:--|:--|---:|---:|:--|---:|---:|---:|
| 0 | A — independent partial sums | 22.95 | 80.6 | verify 0.1% | 1,850.6 | — | 81.6× |
| 1 | B — mutate one slot, update `L` incrementally | 23.31 | 2.1 | verify 3.0% | 50.6 | **36.6×** | 2.23× |
| 2 | C1 — collision-preserving r-adding walk | 22.87 | 1.1 | replay 66.2%, verify 15.0%, setup 1.0% | 139.3 | 13.3× vs A, **0.36×** vs B | 6.15× |
| 2′ | C1 + 8 distinguished-point bits | 37.16 | 1.0 | replay 62.8%, verify 13.0% | 165.0 | 0.84× vs C1; memory **÷257** | 7.28× |
| 3 | C2 — `s_{t+1} = H(L(s_t))` | 21.93 | 80.6 | verify 0.1% | 1,768.8 | 1.05× vs A | 78.0× |
| 4 | R — drop the factor base (plain rho) | 1.33 | 1.0 | replay 59.7%, setup 6.1% | 3.9 | **13.0×** vs B, 35.8× vs C1 | 3.10× (rho floor 1.25) |

Reading down the table:

- **Step 1 is the only large generic win, and it is entirely in `c`.**
  The local-mutation walk keeps `κ` at the birthday value and cuts the
  cost of constructing a residual from two scalar multiplications
  (`≈ 3·log₂n`) to one swap (two additions).  The remaining `2.23×`
  above the floor is that second addition plus verification.
- **Step 2 buys memory, not operations.**  The collision-preserving walk
  walks at one operation per residual but pays two-thirds of its budget
  replaying coefficients after each collision and another 15% verifying
  the dense rows those coefficients produce.  What it makes possible is
  step 2′: with 8 distinguished-point bits the table shrinks `257×` for
  `1.18×` the operations at 28 bits.  At 24 bits the same setting costs
  `1.73×` (`κ` rises to `4×` the floor) because `2^8 = 256` is no
  longer small against the mean walk length `√n/κ ≈ 180`; the rule is
  `2^dp ≪ √n / √(2(B+1))`.
- **Step 3 is a null step.**  The fresh-hash walk pays the full
  independent-sample price per residual and gains only collision
  preservation; with distinguished points its merge-point replay makes
  it `5–6×` dearer still (panel P5).
- **Step 4 is the one that matters, and it is not an optimisation of
  the hybrid.**  Removing the factor base divides `κ` by `17` and `S`
  by `13`.  No amount of tuning in steps 1–3 can reach it: the best
  hybrid sits at `2.23×` its own floor, and that floor is `18×` rho's.

### 9.5 Remaining generic slack, quantified (predictions; measured in §9.7)

These levers are known and generic.  Their predicted effect on `S`
follows from the decomposition above; none touches `κ` except the
negation map, which lowers every strategy's `κ` by the same factor and
so leaves the gap to rho unchanged.  §9.7 applies them and reports the
measured effect next to each prediction.

| lever | applies to | mechanism | predicted S (28-bit, B = 256) | memory |
|:--|:--|:--|---:|:--|
| difference table `P_i − P_j` | B | swap becomes one addition: `c` 2.1 → 1.1 | 50.6 → ≈ 27 (−47%) | `+B²/2` points (1 MB at B = 256) |
| drop verification from the budget | B, C1 | count it as a check, not as work | B −3%, C1 −15% | — |
| store coefficient sketches instead of replaying | C1 | replay is 66% of C1 | 139 → ≈ 47 (−66%) | `× r` counters per stored residual |
| checkpoint replay at the last distinguished point | C1 + dp | replay length ÷ `2^dp` on average | 165 → ≈ 70 | — |
| negation map | all | `κ` ÷ √2 | −29% everywhere | — |
| batched inversions (Montgomery trick) | all | cheaper *operation*, same count | `S` unchanged (wall time only) | — |

Applying every lever above to B gives `S ≈ 27 / √2 ≈ 19`, against a
negation-map floor of `16.0` — and against rho at `3.9` measured, `0.89`
floor.  That is the ceiling of generic tuning: the hybrid cannot close
the last `13×` because that factor is `κ`, not `c`.

### 9.6 The target

For any change to the relation generators, the cells to beat at
`(bits = 28, B = 256, dp = 0)` are:

| tag | S baseline | κ baseline | floor |
|:--|---:|---:|---:|
| B  | 50.6  | 23.31 | 22.67 (16.0 with negation) |
| C1 | 139.3 | 22.87 | 22.67 |
| R  | 3.9   | 1.33  | 1.25 (0.89 with negation) |

- Lowering `S` while `κ/κ_floor` stays at `≈ 1` is engineering: the
  levers above bound it at `≈ 19`.
- **The result that would count as a non-generic advance is
  `κ_total/κ_floor < 0.9` on a hybrid**, where `κ_total` counts walked
  *and* precomputed residuals and the floor is `√(2(B+1)/γ)` for the
  fold order `γ` in use (§10.1), with `correct = yes` on every seed,
  zero relations failing verification, and zero trivial collisions
  counted as relations.  The scoreboard prints this as
  "count invariant beaten".
- Not admissible: changing `B` or `k` (the floor moves with them),
  changing the operation accounting, skipping verification, picking
  seeds, or counting dependent or trivial relations.

After round 2 (§9.7) the bar for *engineering* moves to the tuned
cells, frozen in `experiments/20_residual_walk_tuned.json`:

| tag | variant | S tuned | S/floor | κ/floor |
|:--|:--|---:|---:|---:|
| B  | neg+diff+cont | 19.7 | 1.23 | 1.02 |
| C1 | neg+seg512    | 39.4 | 2.46 | 1.00 |
| R  | neg+seg512    | 1.2  | 1.34 | 0.96 |

The research target is unchanged: `κ/κ_floor < 0.9` on a hybrid.

### 9.7 Round 2: the levers applied and measured

Four generic levers were implemented behind `WalkOptions` flags
(`negation_map`, `diff_table`, `segment_len`, `continue_after_collision`;
bench flags `--negation`, `--diff-table`, `--segment N`, `--continue`,
protocol `--baseline --tuned`) and measured one at a time and together
on the 28-bit, `B = 256` rows of the protocol, three seeds each.  The
fourth lever was not in the §9.5 list; the difference-table ablation
exposed it (a restart after every collision was costing the mutation
walk `0.39` operations per residual, `1,500` restarts per run).

**Single-lever ablation (28-bit, `B = 256`, `dp = 0`).**  `κ` is
`samples/√n`, `c` is walk operations per residual, fractions are of
total operations, `S` is total operations over `√n`:

| lever(s) | tag | κ | c | setup | replay | verify | S | gain vs baseline | walks | trivial | correct |
|:--|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| none (frozen baseline) | A | 22.95 | 80.58 | 0.0% | 0.0% | 0.1% | 1850.6 | 1.00× | 1 | 0 | yes |
| none (frozen baseline) | B | 23.31 | 2.11 | 0.0% | 0.0% | 3.0% | 50.6 | 1.00× | 1526 | 1272 | yes |
| none (frozen baseline) | C1 | 22.87 | 1.06 | 1.0% | 66.2% | 15.0% | 139.3 | 1.00× | 258 | 0 | yes |
| none (frozen baseline) | C2 | 21.93 | 80.59 | 0.0% | 0.0% | 0.1% | 1768.8 | 1.00× | 256 | 0 | yes |
| none (frozen baseline) | R | 1.33 | 1.01 | 6.1% | 59.7% | 0.2% | 3.9 | 1.00× | 2 | 0 | yes |
| negation map | A | 16.10 | 80.59 | 0.0% | 0.0% | 0.1% | 1299.4 | 1.42× | 1 | 0 | yes |
| negation map | B | 16.16 | 2.14 | 0.0% | 0.0% | 4.2% | 36.0 | 1.41× | 1140 | 883 | yes |
| negation map | C1 | 16.17 | 1.09 | 1.3% | 63.2% | 19.0% | 108.1 | 1.29× | 258 | 0 | yes |
| negation map | C2 | 15.64 | 80.59 | 0.0% | 0.0% | 0.1% | 1261.7 | 1.40× | 256 | 0 | yes |
| negation map | R | 0.70 | 1.02 | 11.6% | 56.2% | 0.4% | 2.1 | 1.81× | 2 | 0 | yes |
| difference table | B | 23.31 | 1.36 | 6.5% | 0.0% | 4.2% | 35.6 | 1.42× | 1526 | 1272 | yes |
| continue after collision | B | 23.02 | 1.74 | 0.0% | 0.0% | 2.0% | 41.0 | 1.24× | 1 | 1268 | yes |
| segment 128 | C1 | 22.86 | 1.62 | 3.0% | 10.9% | 7.0% | 46.9 | 2.97× | 2656 | 0 | yes |
| segment 512 | C1 | 22.08 | 1.18 | 3.1% | 25.0% | 14.3% | 45.4 | 3.07× | 751 | 0 | yes |
| segment 2048 | C1 | 22.36 | 1.08 | 2.1% | 44.8% | 18.0% | 68.6 | 2.03× | 328 | 0 | yes |
| negation + difference table + continue | B | 16.28 | 1.00 | 11.7% | 0.0% | 5.9% | 19.7 | 2.57× | 1 | 898 | yes |
| negation + segment 512 | C1 | 15.96 | 1.20 | 3.6% | 28.4% | 19.4% | 39.4 | 3.54× | 588 | 0 | yes |
| negation + segment 512 | R | 0.85 | 1.17 | 16.4% | 2.3% | 0.5% | 1.2 | 3.27× | 26 | 0 | yes |

Predicted versus measured, per lever:

| lever | predicted (§9.5) | measured | note |
|:--|:--|:--|:--|
| negation map | `κ ÷ √2`, `S` −29% on every strategy | `κ` 22.95 → 16.10 (÷1.43); A 1.42×, B 1.41×, C1 1.29×, C2 1.40×, R 1.81× | C1 gains less because its replay and verification do not scale with `κ`; rho gains more because its first collision comes sooner |
| difference table | `c` 2.1 → 1.1, `S` −47% | `c` 2.11 → 1.36, `S` 50.6 → 35.6 (1.42×) | the residual `0.36` above 1.0 was the restart cost, which became lever 4 |
| continue after collision | (not predicted) | `c` 2.11 → 1.74, `S` 1.24×, walks 1,526 → 1 | free: the RNG-driven walk has nothing to restart from |
| segments (C1) | replay −66% at best | 128: 2.97×, **512: 3.07×**, 2048: 2.03× | replay 66% → 25% at 512; shorter segments pay restarts (3%) and lose nothing else |
| all levers, B | ceiling ≈ 19 | **19.7** (2.57×), `c = 1.00` | one operation per residual reached; what is left is the difference table's setup (11.7%) and verification (5.9%) |
| all levers, C1 | — | **39.4** (3.54×) | replay 28%, verification 19%, setup 4% remain |
| all levers, R | — | **1.2** (3.27×), 1.34× its floor | reference moves too: the gap B/R is now 16×, C1/R 33× |

**The tuned protocol against the frozen baseline.**  Every cell
improves, every run recovers the planted `d`, no relation fails
verification, and — the point of the exercise — `κ/κ_floor` stays
within `0.98–1.02` on every hybrid after every lever.  The generic
levers moved `S` exactly as the `κ · c + overheads` decomposition said
they would and never touched the count:

| bits | B | dp | tag | variant | S baseline | S run | improvement (base/run) | κ/floor baseline | κ/floor run | ratio | count invariant beaten? | correct | verdict |
|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|:--|:--|:--|
| 24 | 256 | 0 | A | neg | 1,547.0 | 1,093.3 | 1.42× | 0.99 | 0.99 | 1.00 | no | yes | improvement |
| 24 | 256 | 0 | B | neg+diff+cont | 56.7 | 29.8 | 1.90× | 1.03 | 1.02 | 0.99 | no | yes | improvement |
| 24 | 256 | 0 | C1 | neg+seg512 | 158.4 | 86.6 | 1.83× | 1.00 | 1.00 | 1.00 | no | yes | improvement |
| 24 | 256 | 0 | C2 | neg | 1,528.8 | 1,094.9 | 1.40× | 0.98 | 0.99 | 1.01 | no | yes | improvement |
| 24 | 256 | 0 | R | neg+seg512 | 3.7 | 1.4 | 2.55× | 1.01 | 0.58 | 0.57 | no | yes | improvement |
| 24 | 256 | 8 | C1 | plain | 273.3 | 273.3 | 1.00× | 4.04 | 4.04 | 1.00 | no | yes | unchanged |
| 24 | 256 | 8 | R | plain | 3.9 | 3.9 | 1.00× | 1.07 | 1.07 | 1.00 | no | yes | unchanged |
| 28 | 256 | 0 | A | neg | 1,850.6 | 1,299.4 | 1.42× | 1.01 | 1.00 | 0.99 | no | yes | improvement |
| 28 | 256 | 0 | B | neg+diff+cont | 50.6 | 19.7 | 2.57× | 1.03 | 1.02 | 0.99 | no | yes | improvement |
| 28 | 256 | 0 | C1 | neg+seg512 | 139.3 | 39.4 | 3.54× | 1.01 | 1.00 | 0.99 | no | yes | improvement |
| 28 | 256 | 0 | C2 | neg | 1,768.8 | 1,261.7 | 1.40× | 0.97 | 0.98 | 1.01 | no | yes | improvement |
| 28 | 256 | 0 | R | neg+seg512 | 3.9 | 1.2 | 3.27× | 1.06 | 0.96 | 0.90 | no | yes | improvement |
| 28 | 256 | 8 | C1 | plain | 165.0 | 165.0 | 1.00× | 1.64 | 1.64 | 1.00 | no | yes | unchanged |
| 28 | 256 | 8 | R | plain | 4.0 | 4.0 | 1.00× | 1.08 | 1.08 | 1.00 | no | yes | unchanged |

Scoreboard of the tuned run (the `dp = 8` rows are unchanged because
the negation map, being table-only, is not available with distinguished
points):

| bits | B | dp | tag | variant | seeds | κ = samples/√n | κ floor | κ/floor | c ops/residual | setup | replay | verify | S = ops/√n | S floor | S/floor | ops/rel/√n | stored/√n | trivial | correct |
|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 24 | 256 | 0 | A | neg | 3 | 15.89 | 16.03 | 0.99 | 68.5 | 0.0% | 0.0% | 0.5% | 1,093.3 | 16.03 | 68.20 | 4.265 | 15.818 | 0 | yes |
| 24 | 256 | 0 | B | neg+diff+cont | 3 | 16.34 | 16.03 | 1.02 | 1.0 | 31.6% | 0.0% | 13.5% | 29.8 | 16.03 | 1.86 | 0.117 | 16.195 | 227 | yes |
| 24 | 256 | 0 | C1 | neg+seg512 | 3 | 16.03 | 16.03 | 1.00 | 1.3 | 5.7% | 40.5% | 29.1% | 86.6 | 16.03 | 5.40 | 0.337 | 15.951 | 0 | yes |
| 24 | 256 | 0 | C2 | neg | 3 | 15.92 | 16.03 | 0.99 | 68.5 | 0.0% | 0.0% | 0.5% | 1,094.9 | 16.03 | 68.30 | 4.266 | 15.843 | 0 | yes |
| 24 | 256 | 0 | R | neg+seg512 | 3 | 0.52 | 0.89 | 0.58 | 1.2 | 42.5% | 15.8% | 1.3% | 1.4 | 0.89 | 1.63 | 1.449 | 0.516 | 0 | yes |
| 24 | 256 | 8 | C1 | plain | 3 | 91.57 | 22.67 | 4.04 | 1.1 | 1.8% | 47.1% | 17.3% | 273.3 | 22.67 | 12.06 | 1.064 | 0.097 | 0 | yes |
| 24 | 256 | 8 | R | plain | 3 | 1.34 | 1.25 | 1.07 | 1.0 | 15.9% | 48.5% | 0.5% | 3.9 | 1.25 | 3.11 | 3.903 | 0.004 | 0 | yes |
| 28 | 256 | 0 | A | neg | 3 | 16.10 | 16.03 | 1.00 | 80.6 | 0.0% | 0.0% | 0.1% | 1,299.4 | 16.03 | 81.06 | 5.076 | 16.085 | 0 | yes |
| 28 | 256 | 0 | B | neg+diff+cont | 3 | 16.28 | 16.03 | 1.02 | 1.0 | 11.7% | 0.0% | 5.9% | 19.7 | 16.03 | 1.23 | 0.077 | 16.203 | 898 | yes |
| 28 | 256 | 0 | C1 | neg+seg512 | 3 | 15.96 | 16.03 | 1.00 | 1.2 | 3.6% | 28.4% | 19.4% | 39.4 | 16.03 | 2.46 | 0.153 | 15.943 | 0 | yes |
| 28 | 256 | 0 | C2 | neg | 3 | 15.64 | 16.03 | 0.98 | 80.6 | 0.0% | 0.0% | 0.1% | 1,261.7 | 16.03 | 78.70 | 4.916 | 15.623 | 0 | yes |
| 28 | 256 | 0 | R | neg+seg512 | 3 | 0.85 | 0.89 | 0.96 | 1.2 | 16.4% | 2.3% | 0.5% | 1.2 | 0.89 | 1.34 | 1.192 | 0.848 | 0 | yes |
| 28 | 256 | 8 | C1 | plain | 3 | 37.16 | 22.67 | 1.64 | 1.0 | 0.9% | 62.8% | 13.0% | 165.0 | 22.67 | 7.28 | 0.642 | 0.089 | 0 | yes |
| 28 | 256 | 8 | R | plain | 3 | 1.35 | 1.25 | 1.08 | 1.0 | 5.9% | 60.0% | 0.2% | 4.0 | 1.25 | 3.16 | 3.966 | 0.005 | 0 | yes |

**What round 2 establishes.**  The best hybrid now runs at one group
operation per residual and `1.23×` its generic floor; the r-adding walk
at `2.46×`.  The remaining slack is bookkeeping (verification, setup,
replay), worth at most another `1.2–2×`, and none of it moves `κ`.
Plain rho, tuned with the same two applicable levers, sits at `1.2`
against the hybrids' `19.7` and `39.4`: the factor between them is the
`√(B+1) / √(π/4) ≈ 18×` count ratio, unchanged since round 1.  The next
improvement that matters is not on this list; it has to lower `κ`.

## 10. Round 3: trying to lower κ

Round 2 left the count factor untouched by construction.  This round
attacks it directly, with one candidate on each side of the generic
line.

### 10.1 What "lowering κ" can and cannot mean

Every relation this module collects is a coincidence between two group
elements whose decompositions over `{G, Q} ∪ F` are known: two walked
residuals, a walked residual and a stored one, a residual and `±P_i`.
Let `P` be the number of such elements a run has paid for — walked
residuals *and* anything precomputed and stored to collide against —
and let `γ` be the size of the classes the table is keyed on (1, 2 with
the negation map, 6 with an order-3 automorphism).  In the generic
group model the elements are formal linear combinations of `G, Q, P_i`
with unknown logarithms, and two *distinct* combinations coincide with
probability `1/n` (Shoup's argument: the difference is a non-zero
linear form in the unknowns, which vanishes on at most a `1/n` fraction
of assignments).  Folding by `γ` multiplies the pairs that can coincide
by `γ`, so

```
  E[#relations]  ≤  γ · P(P−1)/(2n)      ⇒      P  ≥  √(2 n R / γ)   for R relations.
```

Hence, for `R = B + 1` independent relations,

```
  κ_total := P / √n  ≥  √(2(B+1)/γ),
```

which is exactly the floor the scoreboard prints — provided `P` counts
*every* element with a known decomposition.  Two consequences shape
what follows.

1. **Precomputation cannot lower `κ_total`.**  A residual that is
   stored before the walk collides exactly like one that is walked
   (each pair coincides with probability `1/n`), so moving work from
   the walk to a seed table changes which term of `γP²/(2n)` the
   relations come from, not their number.  What it *can* do is shrink
   the walked count `κ = samples/√n` — the quantity the scoreboard
   printed until now — which is why the scoreboard now reports
   `κ_total = (samples + seeded)/√n` and compares that to the floor.
   The Semaev `S₃` oracle ("is `L ∈ ±F ± F`?") decides membership in
   the same set of `2B²` elements algebraically instead of from memory;
   in this accounting it is the same event with the same rate and a
   higher cost per residual (`B` square roots instead of one lookup),
   so it is not run.
2. **Only a larger `γ` lowers the floor.**  Within generic operations
   the only handle on the bound is the group of automorphisms the
   fold can use — `±1` on every curve, and additionally `ζ` of order 3
   on `j = 0` curves (order 2 on `j = 1728`).  A structural fold lowers
   `κ` by `√3` for the hybrid *and* for rho, so it moves both sides of
   the comparison and leaves their ratio alone; but it is the one
   honest reduction of `κ` available, and it should be measured rather
   than asserted.

Both candidates are implemented behind flags (`seed_pairs`,
`use_automorphism` on a `generate_j0_instance`), verified by the same
scalar-multiplication check as every other relation, and measured on
the round-2 protocol sizes (`n ≈ 2^24, 2^28`, `B = 256` unknowns, seeds
1–3).

### 10.2 Candidate 1 — seeding the table with every signed pair sum

`--seeded`: the round-2 tuned mutation walk (negation map, difference
table, no restart) with the table pre-filled with the classes of
`±P_i ± P_j` for all `i < j` (`B(B−1)/2 = 32,640` pairs, two classes
each under the negation map, two additions per pair, all counted in
`setup`).  Seeds that coincide with each other or with `±P_k` yield
factor-base-only relations at seeding time, as any collision would.

| bits | B | dp | tag | variant | seeds | κ = samples/√n | seeded/√n | κ_total | κ floor | κ_total/floor | c ops/residual | setup | replay | verify | S = ops/√n | S floor | S/floor | ops/rel/√n | stored/√n | trivial | correct |
|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 24 | 256 | 0 | B | neg+diff+cont+seed2 | 3 | 3.66 | 18.84 | 22.51 | 16.03 | 1.40 | 0.9 | 81.3% | 0.0% | 8.5% | 35.0 | 16.03 | 2.19 | 0.138 | 22.453 | 48 | yes |
| 28 | 256 | 0 | B | neg+diff+cont+seed2 | 3 | 12.15 | 4.60 | 16.76 | 16.03 | 1.05 | 1.0 | 33.8% | 0.0% | 6.9% | 20.4 | 16.03 | 1.28 | 0.080 | 16.691 | 679 | yes |

At 28 bits the walked count drops from `17.09` to `12.15` — the walk
stops `29%` earlier — and the total count is `16.76` against the
unseeded `16.28`: the `4.60·√n` seeds replaced walked residuals one for
one, `S` is unchanged (`20.4` vs `19.7`), and `κ_total/κ_floor` is
`1.05` vs `1.02`.  At 24 bits the seed table alone (`18.8·√n`) already
exceeds the `16.0·√n` residuals the birthday bound needs, so the run
is over-seeded: the walk finishes after `3.7·√n` residuals but the
total is `22.5·√n` and `S` regresses to `35.0`.  Both rows are what
§10.1 predicts to the second digit.  A seed is not a cheaper residual;
it is the same residual paid for earlier.

### 10.3 Candidate 2 — folding by the `j = 0` automorphism

`--structure`: random `j = 0` curves of prime order (`y² = x³ + b`,
`p ≡ 1 mod 3`, `n ≡ 1 mod 3`), the automorphism `ζ(x, y) = (ωx, y) = λ·P`
recovered and checked at generation, a factor base of `256` orbit
representatives (smallest canonical `x`, one point per `⟨±1, ζ⟩`-orbit),
and the three one-op strategies with every round-2 lever, run twice on
the same instances: with negation folding only (`j0`, the control) and
with the six-element fold (`aut6+j0`).  Under the fold a table hit
`L(s) = f·C`, `L(s') = f'·C` gives `f'·L(s) = f·L(s')`, a relation whose
coefficients are `λ`-powers; it is stored in factored form and verified
with two full-size scalar multiplications (`Relation::scaled`).

| bits | B | dp | tag | variant | seeds | κ = samples/√n | seeded/√n | κ_total | κ floor | κ_total/floor | c ops/residual | setup | replay | verify | S = ops/√n | S floor | S/floor | ops/rel/√n | stored/√n | trivial | correct |
|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 24 | 256 | 0 | B | neg+diff+cont+aut6+j0 | 3 | 9.30 | 0.00 | 9.30 | 9.26 | 1.00 | 1.0 | 29.6% | 0.0% | 42.6% | 33.5 | 9.26 | 3.62 | 0.131 | 9.189 | 114 | yes |
| 24 | 256 | 0 | B | neg+diff+cont+j0 | 3 | 16.15 | 0.00 | 16.15 | 16.03 | 1.01 | 1.0 | 31.9% | 0.0% | 16.0% | 31.0 | 16.03 | 1.94 | 0.121 | 16.009 | 209 | yes |
| 24 | 256 | 0 | C1 | neg+seg512+aut6+j0 | 3 | 9.47 | 0.00 | 9.47 | 9.26 | 1.02 | 1.5 | 6.2% | 36.7% | 39.5% | 83.2 | 9.26 | 8.99 | 0.324 | 9.393 | 0 | yes |
| 24 | 256 | 0 | C1 | neg+seg512+j0 | 3 | 16.16 | 0.00 | 16.16 | 16.03 | 1.01 | 1.3 | 5.8% | 40.6% | 29.0% | 88.8 | 16.03 | 5.54 | 0.346 | 16.083 | 0 | yes |
| 24 | 256 | 0 | R | neg+seg512+aut6+j0 | 3 | 0.63 | 0.00 | 0.63 | 0.51 | 1.23 | 1.2 | 39.2% | 12.0% | 4.0% | 1.6 | 0.51 | 3.17 | 1.347 | 0.629 | 0 | yes |
| 24 | 256 | 0 | R | neg+seg512+j0 | 3 | 1.16 | 0.00 | 1.16 | 0.89 | 1.31 | 1.1 | 28.3% | 12.9% | 1.3% | 2.3 | 0.89 | 2.55 | 1.928 | 1.159 | 0 | yes |
| 28 | 256 | 0 | B | neg+diff+cont+aut6+j0 | 3 | 9.47 | 0.00 | 9.47 | 9.26 | 1.02 | 1.0 | 14.8% | 0.0% | 25.1% | 15.7 | 9.26 | 1.70 | 0.061 | 9.411 | 521 | yes |
| 28 | 256 | 0 | B | neg+diff+cont+j0 | 3 | 17.09 | 0.00 | 17.09 | 16.03 | 1.07 | 1.0 | 11.2% | 0.0% | 6.7% | 20.8 | 16.03 | 1.30 | 0.081 | 17.005 | 927 | yes |
| 28 | 256 | 0 | C1 | neg+seg512+aut6+j0 | 3 | 9.41 | 0.00 | 9.41 | 9.26 | 1.02 | 1.2 | 4.1% | 30.4% | 31.9% | 34.9 | 9.26 | 3.77 | 0.136 | 9.387 | 0 | yes |
| 28 | 256 | 0 | C1 | neg+seg512+j0 | 3 | 16.76 | 0.00 | 16.76 | 16.03 | 1.05 | 1.2 | 3.5% | 27.9% | 19.2% | 40.7 | 16.03 | 2.54 | 0.158 | 16.739 | 0 | yes |
| 28 | 256 | 0 | R | neg+seg512+aut6+j0 | 3 | 0.26 | 0.00 | 0.26 | 0.51 | 0.51 | 1.2 | 36.3% | 5.7% | 2.3% | 0.5 | 0.51 | 1.03 | 0.528 | 0.263 | 0 | yes |
| 28 | 256 | 0 | R | neg+seg512+j0 | 3 | 0.58 | 0.00 | 0.58 | 0.89 | 0.65 | 1.2 | 23.3% | 6.3% | 0.7% | 0.9 | 0.89 | 1.02 | 0.900 | 0.577 | 0 | yes |

Per row at 28 bits (three seeds):

| strategy | κ control | κ folded | ratio | floor ratio | κ/floor control → folded | S control → folded |
|:--|---:|---:|---:|---:|:--|:--|
| B  | 17.09 | 9.47 | 1.80 | 1.73 | 1.07 → 1.02 | 20.8 → 15.7 (1.32×) |
| C1 | 16.76 | 9.41 | 1.78 | 1.73 | 1.05 → 1.02 | 40.7 → 34.9 (1.17×) |
| R  | 0.58  | 0.26 | 2.2  | 1.73 | 0.65 → 0.51 | 0.9 → 0.5 (single-collision spread) |

The fold lowers `κ` by the predicted `√3` on every strategy and lowers
the floor by exactly the same factor, so `κ/κ_floor` stays at `1.0`:
this is the count moving *with* its bound, not below it.  `S` improves
less than `κ` because the folded relations are dearer to verify
(`25–32%` of the budget, from `7–19%`) and the difference-table setup
is a larger share of a smaller total.  Rho folds too — `0.9 → 0.5` —
so the hybrid-to-rho ratio is where it was, within rho's spread.

### 10.4 Candidate 3 — the Semaev `S₃` pair oracle

The one non-generic operation available on a prime-field curve is a
summation polynomial.  `S₃(x₁, x₂, x₃) = 0` exactly when some
`(x_i, ±y_i)` sum to `O`, so a residual `L` lies in `±F ± F` if and only
if, for some factor-base `x_i`, the quadratic `S₃(x_L, x_i, X) = 0` has
a root `X` that is itself a factor-base abscissa.  That is a membership
test for the same `2B²`-element set as the signed-pair seed table of
§10.2, done from `B` square roots per residual instead of from memory
(`s3_oracle`; `s3_in_x3` is checked against the crate's `BigUint`
implementation, and `s3_pair_oracle` against brute force over every
signed pair).  Each quadratic solved is charged as one
operation-equivalent — a square root is one field exponentiation,
about the price of an affine addition — which is generous to the
oracle.

Measured on the tuned mutation walk (`--oracle s3`; `n ≈ 2^24, 2^28`,
`B = 256`, three seeds each):

| bits | κ walked | oracle hits | oracle share of budget | S | S tuned (no oracle) | S ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 24 | 6.01 | ≈ 225 | 98.7% | 1,559.8 | 29.8 | 52× worse |
| 28 | 12.37 | ≈ 100 | 99.5% | 3,183.0 | 19.7 | 162× worse |

Two things are true at once.  The walked count does fall — to `0.77`
of the floor at 28 bits and `0.38` at 24 bits, where nearly every
residual decomposes — because each residual is now compared against
`2B²` *virtual* points that were never generated.  And the cost of
those comparisons is `B` operation-equivalents per residual against
one for a table lookup, so the total work is two orders of magnitude
above the plain walk and three above rho.  The oracle checks `2B`
potential coincidences per operation; the residual table, once it
holds `T` entries, checks `T` per operation, and `T ≈ 16√n ≈ 230,000`
at 28 bits against `2B = 512`.  The seed table of §10.2 is the same
oracle with the `2B²` checks paid once instead of per residual, and it
already showed that even at that price the total count does not move.

The scoreboard therefore withholds the "count invariant beaten" flag
on oracle runs and scores them on `S`: their walked `κ` is not a count
of points paid for, and the bound of §10.1 is a bound on points.  For
an `S₃`-style oracle to beat the walk it would have to check more than
`T` coincidences per operation, i.e. decide membership in a set larger
than the residual table at unit cost — which is what a summation
polynomial does *not* do: it trades memory for square roots, one per
factor-base element.

### 10.5 Candidate 4 — triple decompositions (`S₄`, and its meet-in-the-middle form)

Three summands reach a far larger set: signed triples with distinct
indices number `8·C(B,3) ≈ 22·10⁶` at `B = 256`, about a tenth of `n`
at 28 bits, so roughly one residual in ten decomposes outright and no
collision is needed for it.  Two oracles decide membership:

- **Algebraic `S₄`** (`s4_oracle`): `S₄ = Res(S₃, S₃)`, evaluated as
  `S₃` applied twice — for each `i` the roots `Y = x(L ∓ P_i)`, then for
  each `j > i` the roots `X` of `S₃(Y, x_j, X)`, looked up in the base.
  `B²` quadratics per residual, no table at all.  Checked against brute
  force over every signed distinct-index triple.
- **Meet in the middle** (`mitm_neighbours` with `seed_pairs`): the
  `2B` neighbours `L ∓ P_k` are computed (one group operation each) and
  looked up in the table that already holds every `±P_i ± P_j`; a hit
  on a seed is a triple, a neighbour in `±F` is a pair, a hit on a
  walked residual is an ordinary collision with one extra term.  `2B`
  operations per residual plus the `B²` seed table.

Measured on the tuned mutation walk (`--oracle s4`, `--oracle mitm3`;
`n ≈ 2^24, 2^28`, `B = 256`, three seeds each):

| oracle | bits | κ walked | decompositions found | oracle share | S | S tuned | S ratio |
|:--|---:|---:|---:|---:|---:|---:|---:|
| meet in the middle | 24 | 0.02 | ≈ 420 | 40.5% (+59.4% seed table) | 48.2 | 29.8 | 1.6× worse |
| meet in the middle | 28 | 0.15 | ≈ 705 | 92.1% | 90.1 | 19.7 | 4.6× worse |
| algebraic `S₄` | 24 | 0.04 | ≈ 251 | 99.4% | 2,504.8 | 29.8 | 84× worse |
| algebraic `S₄` | 28 | 0.17 | ≈ 262 | 100.0% | 11,153.4 | 19.7 | 566× worse |

The meet-in-the-middle oracle is the strongest count reduction in this
note by a wide margin — at 28 bits a run walks `≈ 2,000` residuals,
`0.15·√n`, and about one in three of them decomposes through a
neighbour (`≈ 700` triples for `253` independent relations; the rest
are dependent, as expected once the rank nears `B + 1`) — and it is
still `4.6×` the cost of the plain walk, `60×` rho's.  The reason is
the same ledger as before: `512` neighbour operations per residual
against one.  Per operation it checks `2B² ≈ 131,000` virtual
coincidences (each neighbour against the whole seed table), which is
close to the `T ≈ 230,000` a walked table lookup checks at these sizes
— so the two are within a small factor of each other, and the walk
wins because its checks come with a stored point that keeps paying.
As `n` grows, `T ∝ √n` outpaces `2B²` unless `B` grows like `n^{1/3}`,
at which point the seed table's `B²` cost is itself `n^{2/3}`: the
whole route is `Θ(n^{2/3})`, the textbook figure for index calculus
with 3-decompositions and a linear-algebra-sized base, and the
measured `S` of `48 → 90` from 24 to 28 bits (`1.9×` for `16×` in `n`,
i.e. `n^{0.23}`) is on its way there.

The algebraic `S₄` is the same test with the seed table replaced by a
square root per pair `(i, j)`: `B²` operation-equivalents per residual
instead of `2B`, i.e. `128×` dearer per residual at `B = 256`, with the
only advantage that nothing is stored.  Measured, it is the lowest walked count in this note — `κ = 0.17` at
28 bits, `≈ 2,400` residuals of which `262` decomposed, one in nine as
predicted — and the highest cost: `S = 11,153`, `566×` the tuned walk
and `9,000×` rho, with the oracle at `100.0%` of the budget.  The two
triple oracles decide the same membership at `2B` versus `B²`
operations per residual, and their `S` differ by that ratio (`90`
against `11,153` is `124×`, against `B/2 = 128`).

### 10.6 Where this leaves κ

| what was tried | κ (walked) | κ_total | κ_total / floor | verdict |
|:--|---:|---:|---:|:--|
| round-2 tuned B, 28 bits | 16.28 | 16.28 | 1.02 | reference |
| + signed-pair seeding | 12.15 | 16.76 | 1.05 | count relabelled, not reduced |
| `j = 0` control (negation only) | 17.09 | 17.09 | 1.07 | reference on the structured curve |
| `j = 0` with 6-fold | 9.47 | 9.47 | 1.02 | κ ÷ 1.8, floor ÷ 1.73; rho ÷ 1.8 as well |
| `S₃` pair oracle | 12.37 | (virtual) | — | walked count relabelled as `B` square roots per residual; `S` 162× worse |
| triple oracle, meet in the middle | 0.15 | (virtual) | — | one residual in three decomposes; `2B` operations each; `S` 4.6× worse, `Θ(n^{2/3})` |
| triple oracle, algebraic `S₄` | 0.17 | (virtual) | — | `B²` square roots per residual; `S` 566× worse |

Three statements now stand on measurement rather than argument:

1. **Precomputation cannot lower `κ_total`.**  Seeding traded walked
   residuals for stored ones one for one and left `S` unchanged; the
   scoreboard now counts seeds, and the only thing the old walked-only
   `κ` would have "shown" is that the walk stopped early.
2. **Structure lowers `κ` only by lowering the floor.**  The `j = 0`
   fold is the sole honest reduction of `κ` found — `√3`, as the
   automorphism group predicts — and it applies to plain rho in the
   same measure.  The generic argument of §10.1 says this is the only
   kind of reduction available to anything that adds, negates and
   compares points: the automorphism group is the whole handle.
3. **The `√(B+1)` gap to rho is not a property of the walk.**  It is
   the number of unknowns the relations have to determine, and every
   lever, structural or generic, has left `κ_total/κ_floor` between
   `1.00` and `1.07`.  Anything that lowers it below `0.9` on a hybrid
   at fixed `B` — the target of §9.6, now stated in `κ_total` with a
   fold-aware floor — has to make two *distinct* formal combinations
   coincide with probability above `γ/n`, which is to say it has to
   compute something about the coordinates that the group law does
   not.  The summation-polynomial oracles are the only candidates of
   that kind on prime fields; §10.4 and §10.5 measure `S₃`, `S₄` and the
   meet-in-the-middle triple oracle and find every one of them behind
   the walk in operations — the count they save is bought with
   per-residual work that grows with the base, and the best of them
   scales as `n^{2/3}`.

## 11. Gaudry's setting: a subspace factor base on `E(F_{p³})`

§10.4–10.5 measured the summation-polynomial oracles where they are
weakest.  On a prime field there is no proper additive subspace, so the
factor base `{x < B}` has no algebraic structure and `S₃` can only be
used one base element at a time.  Gaudry's index calculus (2009) is the
setting where the polynomials earn their keep: `E` over `F_{q^k}`, the
base `F = {P : x(P) ∈ F_q}` — an `F_q`-subspace of abscissae — and the
Weil restriction of `S_{m+1}(x_1, …, x_m, x_R) = 0`, with the `x_i`
unknown in `F_q`, a system of `k` polynomial equations over `F_q` in
`m` unknowns whose solving cost does not depend on `|F|`.  For fixed
`k ≥ 3` this beats rho asymptotically, with constants that grow fast
in `k`.

### 11.1 What was built

`cryptanalysis::gaudry_cubic`, `k = m = 3`:

- `F_{p³} = F_p[t]/(t³ − c)` with an `F_p`-multiplication counter;
  random curves `y² = x³ + ax + b`, `a, b ∈ F_{p³}`, of prime order
  `n ≈ p³` (BSGS over the Hasse interval); the subspace base of the
  `≈ p/2` points with `x ∈ F_p`.
- The **Weil-restricted `S₃` pair test**: `S₃(x_Y, X₁, X₂) = 0` with
  `X₁, X₂ ∈ F_p` unknown is three quadratics over `F_p`; `X₂` is
  eliminated by the explicit `4×4` Sylvester resultant of the first
  two, leaving a degree-`≤ 8` polynomial in `X₁` whose `F_p`-roots are
  found by Cantor–Zassenhaus, completed to `X₂` by the quadratic
  formula, checked against the third component and the base, and
  signed by group arithmetic.  Checked against brute force over every
  signed pair.  Cost `O(log p)` field multiplications, independent of
  the base size — where the prime-field `S₃` oracle needed one square
  root per base element.
- The **triple oracle** in meet-in-the-middle form: for every `P_k`
  in the base, `Y = R ∓ P_k` (one group operation) and the pair test on
  `Y`; distinct-index triples reported once.  Checked on constructed
  triples and verified by arithmetic on every decomposition.
- Relation collection by full decomposition of random `R = aG + bQ`
  until `d` is determined (the same `RelationSystem`), and **rho on the
  same group** for the reference.  Accounting: group operations are
  affine additions in `E(F_{p³})`; oracle work is counted in `F_p`
  multiplications and converted at the measured cost of one affine
  addition on that field (`≈ 63`), so `S = total/√n` is comparable
  across the whole note.

Gaudry's `O(1)` solve of the three-unknown `S₄` system, which removes
the remaining factor `2|F|` from the per-residual cost, was built
afterwards and is measured in §11.4; §11.2–11.3 are the
meet-in-the-middle numbers it is compared against.

### 11.2 Measured

`cargo run --release --example gaudry_cubic_bench -- --protocol`, two
seeds per size, `p ≡ 1 (mod 3)`:

| `p` | `n` | base | residuals | decomposition rate | pair tests / residual | `F_p` mults / pair test | ops / residual | total ops | `S` | rho steps | rho `S` | `S` / rho `S` | wall |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 271 | 2^24.2 | 132 | 815 | 0.152 | 265 | 1,512 | 6,694 | 5.5·10⁶ | 1,224 | 2,952 | 1.16 | 1,057× | 6 s |
| 523 | 2^27.1 | 248 | 1,534 | 0.152 | 496 | 1,557 | 12,835 | 19.6·10⁶ | 1,640 | 11,681 | 1.18 | 1,385× | 20 s |
| 1039 | 2^30.1 | 520 | 2,928 | 0.166 | 1,039 | 1,811 | 30,989 | 90.7·10⁶ | 2,707 | 43,086 | 1.37 | 1,978× | 91 s |
| 2083 | 2^33.1 | 1,068 | 5,555 | 0.180 | 2,136 | 1,856 | 65,170 | 361.7·10⁶ | 3,805 | 158,426 | 1.70 | 2,240× | 361 s |

Every run recovered the planted `d`, for both methods.  Least-squares
exponents over the four sizes: total operations `∝ n^{0.69}`
(prediction `n^{2/3}`), operations per residual `∝ n^{0.38}`
(`2|F| ∝ n^{1/3}` times the slow growth of the pair test, `1,512 →
1,856` `F_p` multiplications as `log p` grows), residuals `∝ n^{0.31}`
(`≈ |F|/rate`).

### 11.3 Reading it

The pieces of Gaudry's argument are all visible, and so is what is
missing:

- **The subspace makes the pair test cheap.**  `≈ 1,500–1,900` `F_p`
  multiplications, about `25–30` affine additions, to decide
  `Y ∈ ±F ± F` for a base of `130–1,070` points — the prime-field
  `S₃` oracle paid one square root *per base element* for the same
  decision (§10.4).  That is the `O(1)`-versus-`O(|F|)` gap the Weil
  restriction buys, measured.
- **The decomposition rate is what the count predicts.**
  `(2|F|)³/6n ≈ 0.15–0.18` of residuals are signed distinct-index
  triples, so `≈ 6` residuals per relation and `≈ 6|F|` residuals in
  all: `815 → 5,555` as `|F|` goes `132 → 1,068`.  The relation count
  is `∝ n^{1/3}`, not `∝ √n` — the count factor `κ` is finally
  *sub-birthday*: `815/√n = 0.20` at 24 bits and `5,555/√n = 0.06` at
  33 bits, falling as `n^{-1/6}`.
- **The cost is still `Θ(n^{2/3})`, because the triple test is still
  a loop over the base.**  Each residual runs `2|F|` pair tests, so
  operations per residual grow like `n^{1/3}` and the total like
  `n^{2/3}`; measured `n^{0.69}`.  Against rho's `n^{1/2}` the ratio
  widens with `n`: `1,057×` at 24 bits, `2,240×` at 33 bits.
- **What Gaudry's `O(1)` solve would change.**  Replacing the
  `2|F|` pair tests by one solve of the three-unknown system costs
  some constant `C₃` per residual instead of `2|F| · 1,800` `F_p`
  multiplications (`4.8·10⁵` at `p = 271`, `4.0·10⁶` at `p = 2083`);
  the total becomes `≈ 6|F| · C₃ ∝ n^{1/3}`, and the crossover with
  rho sits where `6|F| C₃ < 1.25 √n · 63`, i.e. `C₃ < 13 · n^{1/6}`
  `F_p` multiplications — `C₃ < 600` at 33 bits, `C₃ < 5,000` at 50
  bits, `C₃ < 10⁶` at 100 bits.  A resultant cascade on the symmetrised
  system is in the `10⁵–10⁶` range by the degree count of §10; a tuned
  Gröbner solve is what the literature uses.  The linear algebra
  (`|F| ∝ n^{1/3}` unknowns, `n^{2/3}` dense, `n^{1/3+ε}` sparse with
  double large primes) then decides the exponent, which is how
  Gaudry's `Õ(q^{2−2/k})` arises.

So the subspace base does what the prime-field base could not — it
makes the count sub-birthday — and at these sizes it does so at three
orders of magnitude more work than rho, with a scaling exponent that
only improves once the last loop over the base is replaced by an
algebraic solve.  Its constant `C₃` is the number that decides whether
the method beats rho at any size that fits in this module; §11.4
builds the solve and measures it.

### 11.4 Gaudry's `O(1)` solve: the three-unknown `S₄` system

The last loop over the base is removed by solving, per residual, the
symmetrised system directly.  `Solver::Groebner` in
`cryptanalysis::gaudry_cubic` does it in five stages, all in `F_p`
arithmetic with the same multiplication counter:

- **Once per curve** (`SymmetrisedS4::precompute`, `41,360` `F_p`
  multiplications): `S₄(x₁, x₂, x₃, x₄) = Res_X(S₃(x₁, x₂, X),
  S₃(x₃, x₄, X))` is expanded symbolically over `F_{p³}` (the `4×4`
  Sylvester determinant of two quadratics in `X`), then
  rewritten in the elementary symmetric polynomials `e₁, e₂, e₃` of
  `x₁, x₂, x₃` by lex-leading-term reduction.  The result
  `H(e₁, e₂, e₃, x₄)` has at most `175` terms, total degree `≤ 4` in
  the `e`s and `≤ 4` in `x₄`.
- **Weil restriction** (`15` multiplications per term): with `x₄ = x_R
  ∈ F_{p³}` substituted and `e₁, e₂, e₃ ∈ F_p` unknown, the three
  `F_p`-components of `H` are three polynomials of degree `≤ 4` in three
  unknowns over `F_p`.  Their common zeros are the `F_p`-points
  `(e₁, e₂, e₃)`; by Bézout at most `64`.
- **Macaulay matrix** at degree `10` (`252 × 286`: every equation
  times every monomial of degree `≤ 6`), reduced to row echelon form
  over `F_p`; the non-pivot columns of degree `< 10` are the standard
  monomials of the quotient.  When the quotient does not close at
  degree `10` — a normal form of `e₁ · b` is unavailable, or `1, e₁,
  e₂, e₃` is not standard — the degree is raised to `11`, `12`, `13`
  (`Θ(d⁶)` in the reduction).
- **Eigenvalues.**  The multiplication matrix `M_{e₁}` on the
  quotient (`≤ 64 × 64`), its characteristic polynomial by Hessenberg
  reduction, its `F_p`-roots by Cantor–Zassenhaus; for each root `λ` a
  left eigenvector is an evaluation functional, normalised at `1` it
  reads off `(e₂, e₃)`, and the candidate `(λ, e₂, e₃)` is checked
  against all three equations.
- **Splitting.**  `T³ − e₁T² + e₂T − e₃` is factored over `F_p`; a
  triple of roots (with multiplicity) is a triple of abscissae in the
  base; the signs are settled by group arithmetic as in the
  meet-in-the-middle oracle.

Two failure modes are handled rather than hidden.  A residual whose
affine Macaulay matrix never closes by degree `13` (solutions at
infinity in the affine coordinates `e₁, e₂, e₃` — about one residual
in a thousand) is sent to the meet-in-the-middle oracle and its cost
is counted on the same ledger (`fallback_fp_muls`).  A residual with a
repeated eigenvalue is reported (`degenerate_eigenspaces`), and the
eigenvector test still recovers every solution whose functional is a
kernel basis vector.  The solver is checked against the
meet-in-the-middle oracle on every residual of a `p = 271` run
(`--cross-check`: `0` mismatches in `728`), by end-to-end recovery of
the planted logarithm, and by a unit test that compares the two
oracles on random residuals.

`cargo run --release --example gaudry_cubic_bench -- --protocol
--groebner`, the same sizes and seeds as §11.2; `C₃` is the measured
`F_p` cost per residual (`oracle_fp_muls / residuals`), the Macaulay
share is the fraction of `C₃` spent in the row reduction, `retries`
counts residuals redone at degree `11`, `fallback` counts residuals
sent to the meet-in-the-middle oracle:

| `p` | `n` | base | residuals | decomposition rate | `C₃` (`F_p` mults / residual) | Macaulay share | retries | fallback | total ops | `S` | MITM `S` (§11.2) | rho `S` | `S` / rho `S` | wall |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 271 | 2^24.2 | 129 | 965 | 0.127 | 4.79·10⁶ | 92 % | 1 | 1 | 73·10⁶ | 16,456 | 1,224 | 0.81 | 20,321× | 38 s |
| 271 | 2^24.2 | 136 | 789 | 0.166 | 4.72·10⁶ | 92 % | 0 | 0 | 59·10⁶ | 13,276 | 1,224 | 1.51 | 8,810× | 31 s |
| 523 | 2^27.1 | 256 | 1,479 | 0.165 | 4.84·10⁶ | 92 % | 2 | 2 | 114·10⁶ | 9,500 | 1,640 | 1.38 | 6,866× | 59 s |
| 523 | 2^27.1 | 240 | 1,880 | 0.122 | 4.81·10⁶ | 92 % | 2 | 2 | 144·10⁶ | 12,016 | 1,640 | 0.98 | 12,212× | 75 s |
| 1039 | 2^30.1 | 507 | 3,045 | 0.151 | 4.85·10⁶ | 92 % | 4 | 4 | 235·10⁶ | 7,006 | 2,707 | 1.29 | 5,446× | 122 s |
| 1039 | 2^30.1 | 532 | 2,648 | 0.184 | 4.92·10⁶ | 92 % | 6 | 6 | 207·10⁶ | 6,187 | 2,707 | 1.45 | 4,264× | 108 s |
| 2083 | 2^33.1 | 1,088 | 5,555 | 0.181 | 4.81·10⁶ | 91 % | 3 | 3 | 425·10⁶ | 4,466 | 3,805 | 1.62 | 2,759× | 222 s |
| 2083 | 2^33.1 | 1,048 | 5,970 | 0.163 | 4.79·10⁶ | 92 % | 2 | 2 | 454·10⁶ | 4,779 | 3,805 | 1.78 | 2,687× | 239 s |

Every run recovered the planted `d`, for both methods.  Reading it:

- **`C₃` is a constant, as promised.**  `≈ 4.7–4.9 · 10⁶` `F_p`
  multiplications per residual at every size, `≈ 76,000` affine
  additions, of which the Macaulay reduction is `≈ 92 %` and the rest
  is the characteristic polynomial and the eigenvectors.  The
  per-residual cost no longer grows with the base — the
  meet-in-the-middle oracle paid `2|F| · 1,800`, i.e. `4.8·10⁵` at
  `p = 271` rising to `4.0·10⁶` at `p = 2083`.
- **Total work is `∝ n^{0.31}`** (prediction `n^{1/3}`, the
  count of residuals `≈ 6|F|`, measured `∝ n^{0.30}`), against
  `n^{0.69}` for the meet-in-the-middle oracle and `n^{1/2}` for rho.
  The exponent Gaudry's argument needs for the relation phase is now
  measured.
- **The constant is where the crossover was predicted to be
  decided, and it decides against.**  §11.3 gives the condition for
  beating rho as `C₃ < 13 · n^{1/6}`.  With `C₃ = 4.8·10⁶` that is
  `n^{1/6} > 3.7·10⁵`, i.e. `n > 2^{111}`: at the sizes this module
  runs the solve is `≈ 10⁴×` rho, and the ratio to rho falls only as
  `n^{-1/6}` — `≈ 12,800×` at 24 bits, `≈ 2,700×` at 33 bits (per-size means; the per-seed ratios in the table carry rho's own variance).
  Against the meet-in-the-middle oracle the solve breaks even at
  `2|F| · 1,850 ≈ 4.8·10⁶`, `|F| ≈ 1,300`, `p ≈ 2,600`, `n ≈ 2^{34}`
  — just past the largest size measured, where the two oracles cost
  the same (`S` `4,620` vs `3,805`) and the solve's flatter
  exponent takes over from there.
- **What the constant is made of.**  The Macaulay matrix at degree
  `10` has `252` rows and `286` columns because the three Weil
  components are dense quartics; the reduction is `≈ 4.4·10⁶`
  multiplications, `≈ rows · cols · rank`.  A structured solver — an
  `F₄`/`F₅`-style reduction that exploits the sparsity of the shifted
  rows, or the resultant of two of the three quartics in `e₃` first —
  would cut it by a constant factor, not change the picture: the
  crossover needs `C₃` below `600` at 33 bits, below `5,000` at 50
  bits, and the smallest conceivable dense reduction of a `64`-solution
  zero-dimensional system in three unknowns is already `≈ 64³`.
  This is the concrete reason the fixed-`k` Gaudry attack is a
  large-`n` statement: its relation-phase constant is on the order of
  `10⁶` where rho's is `1`.
- **The failures are structural, not a regularity-degree problem.**
  `20` of the `22,331` residuals did not close at degree `10`, and
  none of them closed at `11`, `12` or `13` either: every retry ended
  in the fallback.  Those are residuals whose affine system has
  solutions at infinity (or a positive-dimensional component), which
  no affine Macaulay degree resolves; a homogeneous or saturated
  formulation would, and the fallback cost they incur is `0.04 %` of
  the total, so they are noted rather than chased.  The quotient
  otherwise has Bézout's full `64` dimensions at every residual, and
  each residual yields `≈ 1` `F_p`-rational `(e₁, e₂, e₃)` on
  average, of which about one in six has a cubic that splits — the
  decomposition rate `≈ 0.16` of §11.2 again.
- **The linear algebra is not the bottleneck here, and will be
  later.**  Dense elimination on `|F| + 1 ≈ 1,090` unknowns took
  `0.4 s` of the `222 s` at `p = 2083`, but with `|F| ∝ n^{1/3}` it is
  `∝ n` dense and `∝ n^{2/3}` sparse, both worse than rho; Gaudry's
  `Õ(n^{4/9})` for `k = 3` comes from a double-large-prime variation
  that shrinks the matrix to `∝ n^{1/3+ε}` before elimination, which
  this module does not build.

So the three-unknown solve is built, correct, and base-independent —
the count is sub-birthday (`κ ∝ n^{-1/6}`) and the cost per relation
is a constant — and the constant, measured honestly in the same units
as everything else in this note, places the crossover with rho near
`2^{111}`.  That is the number to improve on: any change to the solver
is scored by `C₃`, and the target is `13 · n^{1/6}`.

## References

- J. M. Pollard, *Monte Carlo methods for index computation (mod p)*,
  Math. Comp. 32 (1978).  The walk.
- E. Teske, *Speeding up Pollard's rho method for computing discrete
  logarithms*, ANTS-III (1998).  r-adding walks.
- P. C. van Oorschot, M. J. Wiener, *Parallel collision search with
  cryptanalytic applications*, J. Cryptology 12 (1999).  Distinguished
  points.
- V. Shoup, *Lower bounds for discrete logarithms and related
  problems*, Eurocrypt 1997.  The generic-group bound.
- A. K. Lenstra, M. S. Manasse, *Factoring with two large primes*,
  Eurocrypt 1990.  Large-prime relation cancellation.
- G. Bisson, A. V. Sutherland, *A low-memory algorithm for finding short
  product representations in finite groups*, Des. Codes Cryptogr. 63
  (2012).  Low-memory meet-in-the-middle over subset sums.
- I. Semaev, *Summation polynomials and the discrete logarithm problem
  on elliptic curves*, ePrint 2004/031; C. Petit, M. Kosters,
  A. Messeng, *Algebraic approaches for the elliptic curve discrete
  logarithm problem over prime fields*, PKC 2016.  The algebraic
  decomposition oracles this note contrasts with.
- P. Gaudry, *Index calculus for abelian varieties of small dimension
  and the elliptic curve discrete logarithm problem*, J. Symbolic
  Comput. 44 (2009).  The subspace factor base on `E(F_{q^k})`, the
  symmetrised `S_{k+1}` system and its Gröbner solve; §11.
- J.-C. Faugère, *A new efficient algorithm for computing Gröbner
  bases (F₄)*, J. Pure Appl. Algebra 139 (1999); B. Mourrain,
  *Computing the isolated roots by matrix methods*, J. Symbolic
  Comput. 26 (1998).  The Macaulay-matrix and multiplication-matrix
  solve used in §11.4.
