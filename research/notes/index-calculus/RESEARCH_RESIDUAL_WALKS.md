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
  (see `research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, `research/notes/index-calculus/RESEARCH_SEMAEV_DECOMPOSITION.md`).

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
is scored by `C₃`, and the target is `13 · n^{1/6}`.  §11.5 and §11.6
are the first two rounds of that.

### 11.5 Lowering `C₃`: the Macaulay step

`C₃` was `92 %` Macaulay reduction, so that is where the round went.
Each step below was measured at `p = 271`, seed `1`, on the same
residual stream (`--cross-check` on, so the residual sequence is
identical and every step's output is compared with the
meet-in-the-middle oracle on every residual); the ledger is in `F_p`
multiplications per residual.

| step | Macaulay | of which forward / normal forms | `C₃` | ratio to baseline | cross-check |
|---|---:|---:|---:|---:|---:|
| §11.4 baseline: Gauss–Jordan to reduced echelon form, count charged per row operation | 4,445,000 | 1,100,000 / 3,350,000 (Jordan phase) | 4,817,000 | 1 | 0 / 728 |
| accounting correction: count only the multiplications performed (zero entries of the pivot row are skipped) | 2,860,000 | 620,000 / 2,240,000 | 3,230,000 | 0.67 | same algorithm |
| forward elimination only, normal forms by memoised back-substitution | 1,155,000 | 610,000 / 540,000 | 1,527,000 | 0.32 | 0 / 728 |
| learned row mask (drop the `30` Koszul-redundant rows found on the first residual) — *not kept* | 1,126,000 | 579,000 / 547,000 | 1,493,000 | 0.31 | `1 %` of residuals miss and redo (`300`-residual run; `1,515,000` without the mask on that stream) |

The first row is the §11.4 number.  The second is the same algorithm
counted honestly: the row-operation loop skipped the zero entries of
the pivot row but charged a full row length, so the §11.4 constant was
overstated by `1.5×` — the correction is recorded as a step because
the ledger has to say so, not as an improvement.  The third row is the
algorithmic change: the reduced echelon form was being built for all
`222` pivot rows (the Jordan phase, three quarters of the reduction)
when only the normal forms of the `64` products `e₁ · b` are needed;
forward elimination plus back-substitution restricted, by memoisation,
to the pivot columns those products actually reach costs a quarter of
that.  The fourth row is the classical Macaulay row selection done
empirically — the `30` rows that reduce to zero on one residual are the
same rows on the next — and it buys `2 %`, because a row that reduces
to zero is cheap to reduce, and loses part of it again to the `1 %` of
residuals where the learned set is wrong and the matrix is redone; it
is not in the code.

One correctness fix came out of the wider cross-checks (`p = 523`,
seed `2`: `1` mismatch in `400`): when two solutions share
`e₁ = λ`, the `λ`-eigenspace of `M_{e₁}ᵀ` is two-dimensional and its
kernel basis vectors are not evaluation functionals, so both solutions
were lost.  The fix restricts the commuting operators `M_{e₂}` and
then `M_{e₃}` to the eigenspace and takes their eigenvectors there; it
costs nothing measurable (the normal forms it needs are already
memoised) and the cross-checks are clean again: `0 / 728` at `p = 271`,
`0 / 400` at `p = 523`, `0 / 800` at `p = 1039`.

Protocol with the new solver (`--protocol --groebner`, the same sizes
and seeds, `experiments/21_gaudry_cubic_c3.json`):

| `p` | `n` | base | residuals | `C₃` | forward / normal forms | fallback | total ops | `S` | §11.4 `S` | MITM `S` | rho `S` | `S` / rho `S` | wall |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 271 | 2^24.2 | 129 | 941 | 1.53·10⁶ | 41 % / 35 % | 1 | 23·10⁶ | 5,139 | 14,866 | 1,224 | 0.81 | 6,346× | 18 s |
| 271 | 2^24.2 | 136 | 773 | 1.53·10⁶ | 40 % / 35 % | 0 | 19·10⁶ | 4,218 | 14,866 | 1,224 | 1.51 | 2,799× | 15 s |
| 523 | 2^27.1 | 256 | 1,461 | 1.55·10⁶ | 40 % / 35 % | 2 | 36·10⁶ | 3,022 | 10,758 | 1,640 | 1.38 | 2,184× | 28 s |
| 523 | 2^27.1 | 240 | 1,880 | 1.55·10⁶ | 40 % / 35 % | 2 | 46·10⁶ | 3,872 | 10,758 | 1,640 | 0.98 | 3,935× | 37 s |
| 1039 | 2^30.1 | 507 | 3,038 | 1.57·10⁶ | 40 % / 35 % | 4 | 76·10⁶ | 2,265 | 6,597 | 2,707 | 1.29 | 1,760× | 60 s |
| 1039 | 2^30.1 | 532 | 2,646 | 1.58·10⁶ | 40 % / 34 % | 6 | 67·10⁶ | 1,992 | 6,597 | 2,707 | 1.45 | 1,372× | 53 s |
| 2083 | 2^33.1 | 1,088 | 5,528 | 1.57·10⁶ | 39 % / 35 % | 3 | 139·10⁶ | 1,459 | 4,623 | 3,805 | 1.62 | 901× | 110 s |
| 2083 | 2^33.1 | 1,048 | 5,967 | 1.57·10⁶ | 39 % / 35 % | 2 | 149·10⁶ | 1,567 | 4,623 | 3,805 | 1.78 | 881× | 119 s |

Every run recovered the planted `d`.  What changed and what did not:

- **`C₃ ≈ 1.53–1.58 · 10⁶`, `3.1×` below §11.4**, total work `∝ n^{0.31}` again,, still a constant
  across sizes; the forward elimination is `≈ 40 %`, the normal forms
  `≈ 35 %`, the characteristic polynomial, eigenvectors and root
  finding the remaining `≈ 25 %` — the Macaulay step is no longer nine
  tenths of the cost, and the next factor of two would have to come
  from all three parts.
- **The solve now beats the meet-in-the-middle oracle inside the
  measured range.**  Break-even is `2|F| · 1,850 ≈ 1.55·10⁶`, i.e.
  `|F| ≈ 420`, `p ≈ 840`, `n ≈ 2^{29}`; at `p = 1039` the solve is
  `S ≈ 2,100` against `2,707`, at `p = 2083` `≈ 1,500` against
  `3,805`, and the gap widens as `n^{1/3}` from there.
- **The crossover with rho moves from `2^{111}` to `2^{101}`.**
  `C₃ < 13 · n^{1/6}` with `C₃ = 1.55·10⁶` gives `n^{1/6} > 1.2·10⁵`,
  `n > 2^{101}`.  A `3×` change in the constant is a `3⁶ ≈ 700×`
  change in the crossover size and ten bits of `n`; every further
  halving of `C₃` is worth six bits.  The ratio to rho at the
  measured sizes is `≈ 4,000×` at 24 bits and `≈ 890×` at 33 bits
  (per-size means).
- **What is left in `C₃`.**  The forward elimination of a `252 × 286`
  matrix of rank `222` with `35`-term rows fills in after the first
  block and costs what dense elimination costs; a structured
  elimination that keeps the degree blocks separate is worth perhaps
  a further `1.5×`, and the `30` redundant rows `≈ 5 %`.  Below that
  the only lever is the size of the system: three quartics in three
  unknowns with `64` solutions is what `S₄` on `E(F_{p³})` is, and a
  `64`-solution zero-dimensional system does not get solved in much
  under `10⁵` field operations by any dense method.  The realistic
  floor for this design is therefore `C₃ ≈ 3–5 · 10⁵`, a crossover
  near `2^{90}`; the note's conclusion stands, with the number
  sharpened.

### 11.6 `C₃`, round 2: Macaulay's row selection

One more step on the ledger of §11.5, measured the same way (`p =
271`, seed `1`, `--cross-check`, same residual stream):

| step | rows | Macaulay (forward / normal forms) | charpoly / eigenvectors / roots | `C₃` | ratio to §11.4 | cross-check |
|---|---:|---:|---:|---:|---:|---:|
| §11.5 result | 252 | 1,155,000 (610,000 / 540,000) | 372,000 | 1,527,000 | 0.32 | 0 / 728 |
| Macaulay's row selection | 226 | 507,000 (153,000 / 354,000) | 273,000 / 58,000 / 33,000 | 879,000 | 0.18 | 0 / 728 |

The learned mask of §11.5 dropped the rows that *happened* to reduce
to zero and gained `2 %`; Macaulay's rule drops the rows that are
*redundant by construction* — after triangularising the three
components so that their grevlex leading monomials are distinct, the
shift `m · f_i` is left out whenever `m` is divisible by the leading
monomial of an earlier component (`26` rows of `252`) — and gains
`1.7×`, because those are the shifts whose rows are reduced through
the longest pivot chains before they vanish.  The rank is `222` in
both cases and the cross-check is clean.  What remains is `17 %`
forward elimination, `40 %` normal forms, `31 %` the characteristic
polynomial of a `64 × 64` matrix (Hessenberg, `≈ 64³`), `10 %`
eigenvectors and roots; none of the three big parts has an obvious
factor of two left in it, and the Krylov alternatives to the
characteristic polynomial cost the same `≈ 64³`.  `C₃ ≈ 0.88 · 10⁶`
is where this design settles: `5.5×` below §11.4, crossover of the
relation phase with rho at `n ≈ 2^{96}`, break-even with the
meet-in-the-middle oracle at `|F| ≈ 240`, `n ≈ 2^{27}`.  §11.7 prices
the other half of the method and says what that `2^{96}` leaves out.

### 11.7 The other half of the method: linear algebra and large primes

§11.4–11.6 priced the relation phase, and the `2^{96}` of §11.6 is a
statement about that phase alone.  The method also has to solve the
relation system, and with `|F| ≈ p/2 ∝ n^{1/3}` unknowns that is where
the exponent really lives.  Gaudry's answer is the double-large-prime
variation: a small base `F' ⊂ F`, relations allowed to carry up to two
points of `F \ F'`, those cancelled by structured elimination, and a
system over `|F'|` unknowns only.  This round builds both halves and
measures them in the same units as the rest of the note.

**What was built** (`GaudryOptions`; bench `--sparse-la`,
`--lp-frac r`, `--lp-rule`, `--max-lp k`, `--protocol-la`):

- **Sequential Wiedemann** over `Z/nZ` — `u64` Berlekamp–Massey, the
  solution checked against every row, the logarithm confirmed by one
  scalar multiplication in the group — preceded by **filtering**: a
  column of weight one determines nothing, so its relation is dropped,
  to a fixed point, and the surplus above square is then trimmed by
  removing the relation that creates the fewest new singletons (the
  step every sieve pipeline runs between collection and the solve;
  compare `cryptanalysis::koblitz_sparse_la::filter_relations`).  The
  incremental elimination of `RelationSystem` counts its
  multiplications too, so the two are comparable.  Multiplications
  modulo `n ≈ p³` are charged `9` `F_p` multiplications (schoolbook)
  when they enter `S`.
- **Large-prime elimination** with arbitrary coefficients: one pivot
  relation per large prime; a new relation is reduced by the pivots of
  its large primes until none is left (a full relation over `F'`) or
  one is left without a pivot (it becomes that prime's pivot).  Pivots
  only ever contain primes that had no pivot when stored, so reduction
  terminates.  Checked against a hidden solution: every combined
  relation still holds.  `|F'| = ⌈|F|^{2/3}⌉` below, Gaudry's rule.

Filtering is not a refinement, it is the difference between working
and not.  A random square selection from a sparse system misses
`≈ e^{-w}` of its columns entirely — coupon collector, `w` the row
weight — and is then singular by construction.  Without it, `p = 1039`
took `445` failed solves, `66,481` residuals and `4.5 · 10⁹`
multiplications modulo `n`; with it, the same instance takes one
solve, `2,867` residuals and `2.5 · 10⁶`, filtering `533` relations to
a `387`-row core.

**Measured** (`--protocol-la --groebner`, the solver of §11.6, two
seeds per size, every run correct, `experiments/21_gaudry_cubic_la.json`;
per-size means):

| `p` | `n` | variant | `|F'|` | unknowns | row weight | residuals | LA (mults mod `n`) | LA share of `S` | `S` | `S` / rho `S` |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 271 | 2^24.2 | dense | 132 | 133 | 4 | 857 | 40,000 | 0.05 % | 2,696 | 2,328× |
| 271 | 2^24.2 | Wiedemann | 132 | 100 | 4.0 | 898 | 170,000 | 0.19 % | 2,830 | 2,443× |
| 271 | 2^24.2 | + large primes | 26 | 24 | 8.9 | 1,465 | 24,000 | 0.02 % | 4,634 | 4,000× |
| 523 | 2^27.1 | dense | 248 | 249 | 4 | 1,670 | 167,000 | 0.10 % | 1,989 | 1,680× |
| 523 | 2^27.1 | Wiedemann | 248 | 181 | 4.0 | 1,754 | 556,000 | 0.32 % | 2,094 | 1,769× |
| 523 | 2^27.1 | + large primes | 40 | 40 | 11.2 | 3,454 | 74,000 | 0.02 % | 4,119 | 3,480× |
| 1039 | 2^30.1 | dense | 520 | 521 | 4 | 2,842 | 885,000 | 0.31 % | 1,236 | 903× |
| 1039 | 2^30.1 | Wiedemann | 520 | 371 | 4.0 | 3,069 | 2,340,000 | 0.74 % | 1,340 | 979× |
| 1039 | 2^30.1 | + large primes | 65 | 64 | 14.0 | 7,860 | 238,000 | 0.03 % | 3,406 | 2,488× |
| 2083 | 2^33.1 | dense | 1,068 | 1,069 | 4 | 5,748 | 7,313,000 | 1.22 % | 898 | 528× |
| 2083 | 2^33.1 | Wiedemann | 1,068 | 783 | 4.0 | 6,215 | 10,428,000 | 1.61 % | 974 | 574× |
| 2083 | 2^33.1 | + large primes | 105 | 104 | 17.9 | 21,871 | 738,000 | 0.03 % | 3,378 | 1,989× |

Fitted exponents (least squares over the four sizes):

| variant | total ops | residuals | linear algebra, in `n` | linear algebra, in the number of unknowns `N` |
|---|---:|---:|---:|---:|
| dense (incremental elimination) | `n^{0.32}` | `n^{0.31}` | `n^{0.85}` | `N^{2.47}` |
| Wiedemann | `n^{0.32}` | `n^{0.31}` | `n^{0.68}` | `N^{2.00}` |
| + double large primes | `n^{0.44}` | `n^{0.44}` | `n^{0.56}` | — |

Reading it:

- **Wiedemann is the wrong tool at these sizes, and the right one
  later.**  Its `N^{2.00}` is exactly the `2N` matrix–vector products
  on rows of weight `4` plus Berlekamp–Massey; the incremental
  elimination is `N^{2.47}`, not `N³`, because the pivot rows of a
  weight-`4` system stay sparse.  The constants put the crossover at
  `N ≈ 11,300` unknowns — `p ≈ 22,600`, `n ≈ 2^{43}` — well past
  anything this module runs.  Filtering also shrinks the system by
  `25–30 %` (`1,069 → 783` unknowns at 33 bits), which is a real
  saving that the dense path cannot use because it consumes relations
  as they arrive.
- **The plain method cannot beat rho at any size, whatever `C₃` is.**
  Relations cost `n^{1/3}` and the linear algebra `n^{0.68}` at best,
  so `S = ops/√n` falls as `n^{-1/6}` while the linear-algebra term
  rises as `n^{+0.18}`.  Extrapolating the two measured terms, `S`
  bottoms out at `≈ 265` around `n ≈ 2^{50}` and rises after that: the
  method's best moment is still `≈ 200×` rho.  The `2^{96}` of §11.6
  was the crossover of the relation phase in isolation, and this is
  what it leaves out.
- **The `n^{4/9}` is real.**  With `|F'| = |F|^{2/3}` the total is
  `n^{0.44}` against the `4/9 = 0.444` of Gaudry's theorem for
  `k = 3`, and the residual count carries it: keeping only
  decompositions with at most two large primes keeps
  `≈ 3|F|^{-1/3}` of them (`16 % → 4 %` measured across the range), so
  residuals grow as `|F|^{4/3}`.  The system shrinks from `1,069` to
  `104` unknowns at 33 bits and its linear algebra by `10×`.  That is
  the exponent the method is famous for, measured.
- **And it is bought at `3.8×` the work at 33 bits** (`S` `3,378`
  against `898`), because `n^{4/9}` beats `n^{1/3}` only after the
  plain method's linear algebra takes over, around `2^{50}` above.
  The gap between the two variants is still widening over the measured
  range, exactly as two exponents `0.44 > 0.32` must.
- **Fill-in is the loose end.**  Row weight grows `8.9 → 17.9` as
  the eliminator combines longer chains, so the large-prime linear
  algebra measures `n^{0.56}`, above its own `4/9`.  From a `0.03 %`
  share that takes until `n ≈ 2^{98}` to matter, but it means the
  `n^{4/9}` as built is a relation-phase exponent, not an end-to-end
  one; sieve implementations cap the merge level for exactly this
  reason, and that is the piece this module does not have.
- **Where the crossover with rho actually is.**  `S` for the
  large-prime variant falls as `n^{4/9 - 1/2} = n^{-1/18}`, measured
  `n^{-0.079}` over this range.  Closing the measured factor of
  `1,989` at 33 bits at that rate takes `≈ 200` doublings of `n`: the
  crossover is somewhere past `2^{230}`.  The exponent is genuine and
  the constant is hopeless, and for a `k = 3` index calculus built out
  of a `64`-solution `S₄` system that is the whole story — which is
  why the attacks that matter in practice either cut the constant by
  orders of magnitude (Joux–Vitse's `F₄`-based variant, decompositions
  into `k − 1` points) or change the target (Weil-descent curves,
  §§ on GHS elsewhere in this repository).

The whole ledger is also drawn on one page in
`docs/index-calculus-scoreboard.html`: every variant of §§9–11 as a
multiple of rho on one logarithmic axis, the `C₃` steps of §11.5–11.6,
and the fitted exponents above against rho's one half.

So the ledger closes where it started, with numbers instead of
adjectives: on `E(F_{p³})` at the sizes this module runs, Pollard rho
costs `S ≈ 1.3` and every index-calculus variant here costs between
`528×` and `4,000×` that; the relation phase is genuinely `n^{1/3}`
and genuinely `O(1)` per residual; the linear algebra is genuinely the
thing that decides the exponent; and the double-large-prime variation
genuinely delivers `n^{4/9}` — asymptotically below rho, and out of
reach by two hundred bits.

### 11.8 The unsolved residuals: the diagnosis in §11.4 was wrong

§11.4 said of the residuals that no Macaulay degree closes:

> Those are residuals whose affine system has solutions at infinity (or
> a positive-dimensional component), which no affine Macaulay degree
> resolves; a homogeneous or saturated formulation would.

That was a guess, and it is wrong in both halves.  This round measured
it, and the class of the entry is **accounting**: the cost does not
move by much, but the reason on the record has to be the right one,
because the fix it implies is a different piece of work.

**What was measured.**  4,000 random residuals at each of `p = 271`,
`523`, `1039`, seed 7, with `GAUDRY_DEBUG_SOLVE` reporting the failing
monomial:

| `p` | residuals | unsolved | rate | failing border monomial at degree 10 |
|---:|---:|---:|---:|:--|
| 271 | 4,000 | 18 | 0.45 % | `e₁ · e₃⁹` |
| 523 | 4,000 | 7 | 0.18 % | `e₁ · e₃⁹` |
| 1039 | 4,000 | 3 | 0.075 % | `e₁ · e₃⁹` |

The rate falls roughly as `1/p`, so this is a codimension-one
degeneracy, not a structural feature of the method.

**Not solutions at infinity.**  The failing residuals' affine ideal is
zero-dimensional by the standard staircase test: every variable has a
pure power among the leading monomials, which is exactly what bounds
the staircase in that direction.  The quotient is finite, of dimension
`61` rather than Bézout's `64`, and stays `61` at degrees 10, 11, 12
and 13.

**Not a positive-dimensional component either.**  Enumerating all
`p³ = 19.9 · 10⁶` triples at `p = 271` on two failing residuals gives
**one** `F_p`-point each, against `0` for the residuals that solve
normally.  A positive-dimensional component defined over `F_p` would
carry `≈ p` points, not one.

**Not non-generic coordinates.**  Galligo says a generic linear change
of coordinates makes the initial ideal Borel-fixed, which would give a
compact staircase.  Substituting `e = A e'` for random invertible `A`
and re-solving rescued **0 of 28** failing residuals across the three
sizes.  The code for it is not kept; the negative result is.

**What it actually is.**  The multiplication matrix needs the normal
form of `e₁ · b` for every standard monomial `b`.  A column is
reducible only if it is a pivot, and `standard` is defined as the
non-pivot columns of degree `≤ degree − 1`, so a product landing on a
*non-pivot column of degree exactly `degree`* is neither reducible nor
standard and has no normal form at all.  At degree 10 the standard set
contains `e₃⁹`, and `e₁ · e₃⁹` is one of the three non-pivot columns of
degree exactly 10.  Raising the degree reproduces the same situation
one rung up — at degree 11 the standard set contains `e₃¹⁰` and
`e₁ · e₃¹⁰` is the blocker — which is why §11.4 already observed that
every retry ends in the fallback.  It is the **border of the staircase
meeting the degree cut**, a truncation artefact of the fixed-degree
Macaulay formulation, not a property of the variety.

The real fix is therefore a border-basis or a proper F4/F5 termination
criterion, not saturation.  That is a larger piece of work than the
`0.04 %` of total cost it would recover, so it stays unbuilt — but it
is now the right unbuilt thing.

**What was taken, since the retries are provably useless.**  The border
condition is detectable before the multiplication matrix is built, from
the pivot set alone.  Detecting it and going straight to the fallback
removes three `Θ(d⁶)` reductions on exactly those residuals.  Same
residual stream, same seeds, `GAUDRY_BORDER_RETRY=1` restores the old
behaviour:

| `p` | `C₃` with retries | `C₃` without | saving | unsolved, both |
|---:|---:|---:|---:|---:|
| 271 | 888,392 | 879,363 | 1.02 % | 18 |
| 523 | 896,243 | 892,719 | 0.39 % | 7 |
| 1039 | 905,207 | 903,692 | 0.17 % | 3 |

Class: **engineering**, and a small one — the ratio to the floor does
not move, and the saving falls as `1/p` because the failure rate does.
The reason it is worth committing is not the 1 %: it is that
`the_border_check_sends_only_the_residuals_no_degree_solves_to_the_fallback`
runs the same 600-residual stream both ways and asserts the *answers*
are identical, so the check is pinned as losing nothing.  `xcheck` on
the `p = 271` protocol run is `0/728` against the meet-in-the-middle
oracle, with `retries=0` and `fallback=1` unchanged.

### 11.9 Joux–Vitse: decompositions into `k − 1` points

§11.7 closed by naming the two things that make an index calculus
matter in practice, one of which was *"Joux–Vitse's `F₄`-based variant,
decompositions into `k − 1` points"*.  Built and measured here.

**What it is.**  Instead of asking whether a residual is a sum of `k = 3`
factor-base points, ask whether it is a sum of `2`.  The solve then
collapses: `S₃(x₁, x₂, x_R) = 0` Weil-restricts to three quadratics in
two `F_p` unknowns, so a resultant of two conics and Cantor–Zassenhaus
settle it — no Macaulay matrix, no eigenvalues, no characteristic
polynomial.  `weil_s3_pair_test` already existed as the inner loop of
the meet-in-the-middle oracle; `Solver::PairOnly` calls it once on the
residual itself.

**The trade.**  A random residual is a pair far less often than it is a
triple — `≈ 2|F|²/p³` against `≈ |F|³/6p³`, a factor `≈ p/12` fewer — so
the residual count goes from `Θ(|F|) = Θ(n^{1/3})` to
`Θ(p³/|F|) = Θ(p²) = Θ(n^{2/3})`.

**Measured.**  Same instances, same seed, same accounting as §11.5; every
row recovered the planted `d`:

| `p` | `n` | variant | residuals | decomp | `F_p` muls / residual | total ops | `S` | rho `S` | `S` / rho |
|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|
| 271 | 2^24.2 | three-point | 941 | 123 | 873,241 | 13.1·10⁶ | 2,939 | 0.8 | 3,674× |
| 271 | 2^24.2 | **pair-only** | 20,832 | 38 | **1,513** | 1.9·10⁶ | **427** | 0.8 | **534×** |
| 523 | 2^27.1 | three-point | 1,461 | 245 | 890,384 | 20.8·10⁶ | 1,738 | 1.4 | 1,241× |
| 523 | 2^27.1 | **pair-only** | 157,837 | 143 | 1,557 | 15.9·10⁶ | **1,331** | 1.4 | **951×** |
| 1039 | 2^30.1 | three-point | 3,038 | 461 | 903,689 | 43.9·10⁶ | **1,312** | 1.3 | **1,009×** |
| 1039 | 2^30.1 | pair-only | 486,963 | 231 | 1,810 | 55.4·10⁶ | 1,655 | 1.3 | 1,273× |
| 2083 | 2^33.1 | three-point | 5,528 | 1,004 | 919,646 | 82.3·10⁶ | **866** | 1.6 | **541×** |
| 2083 | 2^33.1 | pair-only | 1,276,421 | 323 | 1,857 | 157.8·10⁶ | 1,659 | 1.6 | 1,037× |

Fitted exponents over the four sizes:

| variant | residuals | total ops | `S` |
|---|---:|---:|---:|
| three-point (§11.6 solver) | `n^{0.29}` | `n^{0.30}` | `n^{-0.20}` |
| pair-only | `n^{0.67}` | `n^{0.72}` | `n^{+0.22}` |

Reading it:

- **The `2/3` is exact.**  The predicted residual exponent for `k − 1`
  decompositions at `k = 3` is `2/3`; measured `0.666` over four sizes.
  That is the cleanest confirmation of a predicted exponent in this
  note.
- **The constant is `580×` better and the exponent is worse.**
  Per-residual cost falls from `≈ 0.9·10⁶` `F_p` multiplications to
  `≈ 1.5·10³`.  That is the whole Joux–Vitse promise, delivered.  But
  `S` *rises* as `n^{0.22}` where the three-point solve's *falls* as
  `n^{-0.20}`, so the win is bounded.
- **They cross inside the measured range**, between `p = 523` and
  `p = 1039`, near `n ≈ 2^{28.5}`.  Below it pair-only is the best
  oracle in this module — `534×` rho at 24 bits against the three-point
  solve's `3,674×`, a factor of `6.9`.  Above it, it loses, and by 33
  bits it is `1.9×` worse.
- **Why it does not scale here.**  The variant's saving is the collapse
  of the Gröbner step, and at `k = 3` that step is a `64`-solution
  system in three unknowns — expensive, but only by a constant.  The
  variant pays an honest `Θ(p)` in decomposition rate for it.  Trading a
  constant for a factor of `p` is a good trade exactly while the
  constant dominates, which is what the crossover at `2^{28.5}` is.  The
  `k − 1` idea is built for the regime where the `k`-point Gröbner step
  is *exponentially* expensive, and `k = 3` is not that regime.

Class: **advance below `2^{28.5}`, relabelling above it** — `S` genuinely
falls against every other oracle here at the small sizes, and genuinely
rises at the large ones; the ratio to rho never approaches `1` at any
size, so this changes the constant and not the conclusion.  Rho still
costs `S ≈ 1.3` and the best row in this table is `427`.

One property worth recording because it is easy to misread: pair-only
does **not** compute the factor-base logarithms.  It stops when the `d`
column alone becomes determined — the same `solved()` criterion the
other two solvers use, so the comparison is like for like — which
happens when a combination of relations cancels every base column.  With
weight-2 relations that is a cycle, and it arrives long before the
`|F|`-row system is anywhere near full rank: `38` relations against
`130` unknowns at `p = 271`.

### 11.10 Pre-registration: does a merge-level cap restore `n^{4/9}`?

**Written before the experiment was built or run.**  §11.7 leaves one loose
end and names it: the large-prime linear algebra measures `n^{0.56}` against
its own `4/9`, because row weight grows `8.9 → 17.9` as the eliminator chains
merges, and "sieve implementations cap the merge level for exactly this reason,
and that is the piece this module does not have."  This registers the
experiment that supplies it.

**Where the excess comes from, derived.**  The large-prime variant runs
Wiedemann, so its linear algebra costs `≈ N² w` — `2N` matrix-vector products
of `N w` nonzeros each — with `N` unknowns and mean row weight `w`.  The small
base is Gaudry's rule, `N = |F|^{2/3}` (`SmallBase::Rule`), and `|F| ~ n^{1/3}`,
so `N ~ n^{2/9}` and

```text
    N² ~ n^{4/9}                    exactly the relation-phase exponent
    w  ~ n^{0.113}                  measured, 8.9 → 17.9 over 2^24.2 → 2^33.1
    N² w ~ n^{0.557}                against the fitted n^{0.56}
```

**Fill-in is therefore the entire excess.**  Not part of it — all of it.  A cap
that holds `w` to a constant puts the linear algebra exactly on `n^{4/9}`.

**What it costs.**  `LargePrimeEliminator::feed` loops until every large prime
is cancelled, subtracting one stored pivot per step and merging that pivot's
columns in.  A cap abandons a relation once it has been reduced against `k`
pivots.  Every abandoned relation is a residual that was paid for and thrown
away, so the relation phase pays for the linear algebra's saving.  **The
question this experiment answers is whether that repayment is a constant factor
or a growing one**, because only the first leaves `n^{4/9}` end to end.

**The falsifier, and the three outcomes.**  Fit the end-to-end exponent `e`
over the same four sizes, same seeds, same `--protocol-la --groebner` as
`experiments/21_gaudry_cubic_la.json`, with correctness preserved on every run:

| outcome | `e` | what it means |
|---|---|---|
| **a** | `e ≤ 0.444` | the cap restores `n^{4/9}` end to end; §11.7's loose end closes, and the exponent is the method's rather than the relation phase's |
| **b** | `0.444 < e ≤ 0.56` | the cap buys part of the gap; `n^{4/9}` stays a relation-phase exponent and the note says so as a measurement |
| **c** | `e ≥ 0.56` | the cap does not help |

**Predicted: (b), and (a) is live.**  Capping cannot raise `w` and cannot lower
`N`, so `e` cannot exceed the uncapped `0.56` except through repayment; and it
cannot fall below `4/9` at all.  Which of (a) and (b) lands depends entirely on
how the discarded-residual count scales, which is the thing being measured and
is not predicted here.

**`S` will not move, and that is not a failure.**  The linear algebra is
`0.02–0.03 %` of `S` at these four sizes, so any cap changes total cost by well
under `1 %` and the `S / rho` column stays at its published `1,989×` to
`4,000×`.  The deliverable is a fitted exponent, not a cheaper attack, and a
flat `S` column is what a correct run looks like.  By §11.7's own arithmetic
the exponent does not start paying until `n ≈ 2^{98}`.

**Inadmissible**, by §6 and the standing rules of this note: changing the four
sizes or the two seeds; changing `SmallBase::Rule`; reporting the linear-algebra
exponent alone as the end-to-end one, which is the mistake §11.7 exists to
record; quoting an improvement in `S` from a cap whose LA share is `0.03 %`;
and counting a run whose recovered logarithm was not checked.

### 11.11 The merge-level cap, measured: the phase improves, the method does not

**Runner:** `cargo run --release --example gaudry_cubic_bench -- --protocol-la
--groebner --sizes 271,523,1039,2083 --seeds 6 --merge-cap K --json …`
**Frozen:** `experiments/22_gaudry_merge_cap_uncapped.json`,
`experiments/22_gaudry_merge_cap_k12.json` (six seeds), and
`…_k10_2seed.json`, `…_k16_2seed.json` (the bracketing caps, two seeds)
**Summary:** `python3 scripts/summarize_merge_cap.py <uncapped> <capped>`
**Registered in advance:** §11.10.  Every one of the 48 runs recovered its
planted logarithm.

§11.7 named the missing piece — "sieve implementations cap the merge level for
exactly this reason, and that is the piece this module does not have".  It is
now built (`GaudryOptions::max_merge_level`) and measured.  **It does not
help.**

**Paired, because the noise is bigger than the effect.**  The fitted exponent
moves by about `0.04` between seeds, which is the size of what a cap does, so a
difference of means across two arms says almost nothing — and a two-seed
version of this table said the opposite of the truth.  Each seed is therefore
run capped and uncapped and the statistic is the mean of the per-seed
differences, which cancels the spread.

| quantity | uncapped | `k = 12` | paired `Δ` | `±` | worse/better |
|---|---:|---:|---:|---:|---:|
| **end-to-end total** | 0.424 | 0.475 | **`+0.051`** | 0.031 | 5/1 |
| linear algebra | 0.516 | 0.507 | `−0.009` | 0.019 | 2/4 |
| — of which **solve** | 0.542 | 0.498 | **`−0.045`** | 0.017 | 1/5 |
| — of which **merge** | 0.382 | 0.425 | `+0.043` | 0.030 | 5/1 |
| residuals | 0.418 | 0.468 | `+0.051` | 0.031 | 5/1 |

**The cap does exactly what it was designed to do.**  The solve exponent falls
`0.045 ± 0.017`, `2.6σ`, on five of six seeds; row weight at `p = 2083` drops
from `20.93` to `11.89`.  The sparser matrix is real.

**And the eliminator cancels it exactly.**  Merge work rises `+0.043 ± 0.030`,
so the linear algebra as a whole moves `−0.009 ± 0.019` — indistinguishable
from zero, on a 2/4 sign split.  Capping does not lower the linear-algebra
exponent.  It moves work from the solve into the eliminator.

**End to end it is worse**, `+0.051 ± 0.031` on five of six seeds, and the
residual row is identical to it: the total is the relation phase, and the cap
damages the relation phase.

### What the cap actually costs

Not the relations it discards.  Those are `0.3 %` to `2.0 %` of residuals,
which could never move an exponent.  The cost is the full relations that
**never form**, because a chain the cap truncates is a relation that would have
closed:

| `p` | full relations per residual, uncapped | capped | yield |
|---:|---:|---:|---:|
| 271 | 0.0192 | 0.0182 | `0.95×` |
| 523 | 0.0123 | 0.0114 | `0.93×` |
| 1039 | 0.0105 | 0.0067 | **`0.64×`** |
| 2083 | 0.0073 | 0.0046 | **`0.64×`** |

**And the loss grows with `n`** — `5 %` at the small sizes, `36 %` at the two
large ones.  That is the whole mechanism: a yield loss that grows in `n` is a
residual exponent that rises, which is what the table above measures.  Mean
merge depth is roughly flat across the ladder (`9.1, 9.8, 12.2, 10.6`), so
chains are not getting longer on average; **the tail past the cap thickens**,
and a fixed cap excludes a growing share of exactly the chains that close.

### Against §11.10's registered outcomes

The capped end-to-end exponent is `0.475`, which falls in band **(b)**,
`0.444 < e ≤ 0.56`.  But reading it as "the cap buys part of the gap" would be
wrong, and the registration is what is at fault: its three bands were written
assuming a cap could only help or do nothing.  **The uncapped exponent is
`0.424 ± 0.016`** — already below `4/9`, so there was no end-to-end gap to buy,
and the cap moved the number the wrong way.  Recorded here rather than quietly
re-banded: the prediction (b) was met by coincidence of arithmetic, not because
the experiment came out as expected.

**Two things §11.7 got right and one it got wrong.**  Right: the linear-algebra
exponent is genuinely above `4/9` — six seeds give `0.516 ± 0.027`, `2.7σ`
above — and fill-in is why.  Right: the module lacked the cap.  Wrong: the
implication that the cap was therefore the fix.  It is not, and the
`n^{0.56}` that motivated it was itself a two-seed number quoted without an
error bar; the six-seed value is `0.516 ± 0.027`, and `0.56` sits inside a
two-seed `±0.046`.

**The `n^{4/9}` is end-to-end here after all.**  §11.7 worried it was "a
relation-phase exponent, not an end-to-end one".  Measured over these four
sizes the total is `0.424 ± 0.016`, consistent with `4/9`, because the linear
algebra is only about `2 %` of the cost.  The higher linear-algebra exponent is
real and will eventually dominate — §11.7's own arithmetic puts that near
`n ≈ 2^{98}` — but it does not bite on this ladder, and capping the merge level
is not how to meet it when it does.

**Class: `engineering`, and negative.**  The algorithm gained a lever, the
lever was measured, and it costs more than it saves.  `S` is untouched at these
sizes, exactly as §11.10 said in advance it would be, so nothing here moves the
standing against rho.

### 11.12 Pre-registration: a border basis in place of the fixed-degree Macaulay cut

**Written before the border basis existed.**  The scoreboard has carried "the
open direction is a border basis in place of the fixed-degree Macaulay cut"
since §11.6.  This registers it, and registers first the arithmetic that says
what class of result it can possibly be — because that arithmetic is available
now, from measurements already in this note, and stating it afterwards would be
worthless.

**It cannot be an advance.  It is `engineering` by construction.**  §11.5's
protocol table measures `C₃` at `1.53, 1.53, 1.55, 1.55, 1.57, 1.58, 1.57,
1.57` (`×10⁶`) across `n = 2^{24.2}` to `2^{33.1}`.  **`C₃` is flat in `n`** —
the `S₄` system has three unknowns and 64 solutions at every size, so its cost
does not scale.  A lever on a quantity that does not scale cannot move an
exponent, and by §3 of `AGENTS.md` that is `engineering`: "legitimate, bounded,
and not a finding".

**Its ceiling is about `1.2×`, and an earlier revision of this section said
`2.3×`.**  That was wrong and the error is worth stating, because it is the
kind that makes a registration flatter itself.  §11.6 splits the
post-row-selection `C₃` into forward elimination `17 %`, normal forms `40 %`,
characteristic polynomial `31 %`, eigenvectors and roots `10 %`.  The first
revision claimed a border basis replaces the first two.  **It cannot replace
the normal forms: a border basis *is* the normal forms of the border
monomials.**  That 40 % is the thing being computed, not overhead being
removed.  Nor does the charpoly of the `64 × 64` multiplication matrix or its
eigen-solve care how the matrix was reached, so `31 %` and `10 %` survive
untouched.

What is actually available is the forward elimination — and only the part of it
spent on rows that never reach the border, since Macaulay's row selection
(§11.6) has already dropped the 26 Koszul-redundant ones.  Driving **all** of
it to zero gives `1/(1 − 0.17) ≈ 1.2×`.  The realistic expectation is `≈ 1×` or
worse, because reducing the border monomials still needs enough of the row
space to do it.  Registered at `1.2×` so that a measured `0.9×` reads as this
section having been right about the ceiling, not as the experiment failing.

**And a far larger constant has already failed to matter.**  Decomposing into
`k − 1` points cut the per-residual constant `580×`, from `0.9` million field
multiplications to `1,513`, and still landed `1,037×` above rho.  Against
§11.7's remaining `1,989×` closing as `n^{-1/18}` — about two hundred doublings
— a `1.2×` is worth roughly a third of one.  **The verdict does not move, and this section
says so in advance so that a `C₃` improvement cannot later be read as one.**

**What it might fix that is not a constant.**  `solve_at_degree` carries a
documented failure: the multiplication matrix needs the normal form of `e₁ · b`
for every standard `b`, and a product landing on a non-pivot column of degree
exactly `degree` has none.  The note records that this "recurs identically one
degree up, which is why every retry in the measured runs ended in the fallback
and none in a solution".  That is the staircase's border truncated by a
fixed-degree cut, and closing the border is exactly what a border basis does.
Whether it removes the fallbacks is a **correctness** question, separate from
the constant, and is registered as its own outcome below.

**The staircase, measured before the build** (`GAUDRY_DEBUG_SOLVE=1`,
`p = 271`, seed 1, 20 residuals).  It decides how much algorithm this needs,
so it was looked at rather than assumed:

```text
  degree 10: rows 226 cols 286 pivots 222 standard(dim) 64
             maxima [3,4,9] box=false        — identical on all 20
```

The order ideal is **the same 64 monomials on every residual**, and it is not
the box `[0,4)³` that `64 = 4³` invites you to guess.  It is
`z < 2(5 − x − y)`:

| `x` | admissible `z` by `y = 0, 1, 2, 3, 4` |
|---:|---|
| 0 | `<10`, `<8`, `<6`, `<4`, `<2` |
| 1 | `<8`, `<6`, `<4`, `<2` |
| 2 | `<6`, `<4` |
| 3 | `<4` |

Two things follow.  **Mourrain's iteration is not needed**: a residual-
independent order ideal means the staircase and its border can be computed once
per curve rather than rediscovered on each of the 941 residuals, which is the
only structural saving on offer here.  And of the 64 products `x · b`, exactly
**30** leave `O` and need a border normal form; the other 34 are shifts within
it.

**The coverage outcome is below the resolution of this cell.**  Zero
`border_unreachable` events in those 20 residuals, and §11.5 recorded one
fallback in 941 at `p = 271` seed 1 — `0.1 %`.  Whatever the border basis does
to the fallbacks cannot be established here, and this section will not claim it
was.

**The falsifier.**  Measured on §11.5's protocol — same `p ∈ {271, 523, 1039,
2083}`, same seeds, `--cross-check` on so the residual stream is identical and
every output is compared against the meet-in-the-middle oracle on every
residual:

| outcome | condition |
|---|---|
| **success (engineering)** | `C₃` below §11.6's `0.88 × 10⁶` with **zero** cross-check mismatches |
| **coverage gain** | the fallback count reaches zero where the Macaulay cut had 1–6 per run.  *Coverage*, not correctness: a fallback residual is not a wrong answer, it is one the fixed-degree cut could not solve algebraically and the meet-in-the-middle oracle solved correctly but more slowly |
| **failure** | `C₃` at or above `0.88 × 10⁶`, or any cross-check mismatch |

A mismatch is disqualifying on its own, whatever the cost column says: a
cheaper solver that returns a wrong decomposition is not a cheaper solver.

**Inadmissible**, by §6: changing the sizes, seeds or residual stream; turning
`--cross-check` off; quoting the `C₃` improvement as a change in `S / rho`
beyond the same factor; and reporting a fallback reduction as a cost result,
since the fallbacks are 1–6 residuals of thousands and cannot move `C₃` either
way.  Also inadmissible: describing the fallbacks as unsoundness in the current
solver.  They are residuals it declines and hands on, and the answer that comes
back is right.

### 11.13 The border basis, measured before it was written: a `1.06×` ceiling

> **Corrected in §11.14 — class `accounting`.**  The `5.6 %` below charges each
> elimination multiplication to the pivot row doing the subtracting.  What an
> elimination that builds only the rows it needs can skip is the work done *on*
> rows nothing reads, and the rows the normal forms depend on — the ones they
> read and, transitively, every pivot row subtracted from those — are **202 of
> 222**, not 149.  The skippable work is `1.1 %` of the elimination and `0.19 %`
> of `C₃`: the ceiling is **`1.002×`**, not `1.06×`.  The staircase is also not
> identical on every residual, only on `98 %` of them at `p = 271`.  The section
> is left as it was written; its figures are the "before" marks of §11.14.

**Instrumentation:** `GAUDRY_DEBUG_SOLVE=1`, `p = 271`, seed 1, `--cross-check`,
25 residuals — §11.5's cell, so these numbers sit beside §11.6's directly.
`solve_at_degree` now prints the staircase, how much of the echelon the normal
forms reach, and the elimination cost attributed per pivot column.

§11.12 registered this as `engineering` by construction and put its ceiling at
`1.2×`.  Three measurements taken before writing the solver put it at **`1.06×`**,
and that is the result.

**1. The staircase is fixed and is not the box.**  Identical on all 25
residuals — `rows 226, cols 286, pivots 222, dim 64` — with shape
`z < 2(5 − x − y)`, maxima `[3, 4, 9]`.  So no Mourrain iteration is needed, and
of the 64 products `x · b` exactly 30 leave the order ideal.

**2. The normal forms reach `67 %` of the echelon.**  `149` of `222` pivots,
`183` of `286` columns, again identical on every residual.  A third of the
elimination produces pivot rows that nothing afterwards consults.

**3. That third is `5.6 %` of `C₃`.**  Attributing elimination multiplications
to the pivot column that caused them:

| | per residual | of `C₃ ≈ 0.88 × 10⁶` |
|---|---:|---:|
| elimination, total | `≈ 150 300` | `17.0 %` |
| — on pivots the normal forms never reach | `≈ 49 000` | **`5.6 %`** |

`32.4 %` to `32.8 %` of the elimination across 25 residuals, a very tight band.
The `17.0 %` is worth noting on its own: it reproduces §11.6's `17 %` forward-
elimination share from an independent counter, so the split that ceiling rests
on is confirmed rather than assumed.

**The ceiling, therefore, is `1/(1 − 0.056) = 1.06×`** — and that is an
*upper* bound reached only by an elimination that skips every unreached pivot
at zero cost.  A real lazy elimination cannot: the unreached columns are zero
in the reached rows *because* the elimination zeroed them, so skipping a pivot
leaves live entries below it and the dependency has to be tracked rather than
assumed away.  The achievable figure is below `1.06×`.

**Why the solver was not then written.**  §11.12's falsifier asks for `C₃`
below `0.88 × 10⁶`, and §11.12's stop condition says not to reach for a harder
algorithm when the minimal one does not clear it.  A `1.06×` ceiling clears it
by `5.6 %` at most, on a phase that §11.5 measures as **flat in `n`** and that
§11.7 shows is `2 %` of a method sitting `1,989×` from rho.  The implementation
would confirm a number already bounded by the repo's own counters, at the price
of a lazy elimination whose correctness surface is the part of this that could
actually go wrong.

**What this closes.**  The scoreboard has carried "the open direction is a
border basis in place of the fixed-degree Macaulay cut" since §11.6.  It is not
open any more, and it did not need the build to close it: **the direction is
worth at most `1.06×`, measured.**  §11.5's own description of what it already
does — "forward elimination plus back-substitution restricted, by memoisation,
to the pivot columns those products actually reach" — is a border-basis
computation in all but name, which is why so little is left.  The `5.6 %` is
the gap between *lazy back-substitution*, which this solver has, and *lazy
elimination*, which it does not.

**Class: `engineering`, negative, and `0 / 25` cross-check mismatches** on the
instrumented runs — the diagnostics do not touch the arithmetic.

### 11.14 The border basis, written: the ceiling was `1.002×`, and the solver gets 99 % of it where it runs

**Runner:** `GAUDRY_LAZY_VERIFY=1 cargo run --release --example gaudry_cubic_bench --
--protocol --groebner --cross-check --sizes 271,523,1039,2083 --seeds 2 --json …`,
and the same without the variable for the baseline.
**Frozen:** `experiments/23_gaudry_lazy_elim_baseline.json` and
`experiments/23_gaudry_lazy_elim_verify.json` (with their `.log`),
`experiments/23_gaudry_lazy_elim_lazy.json` (the same run under
`GAUDRY_LAZY_ELIM=1`, whose operation counts equal the verified run's in every
cell — verification is uncharged), and `experiments/23_gaudry_lazy_elim_closure.log`
(the per-residual closure and lazy-cost lines).
**Summary:** `python3 scripts/summarize_lazy_elim.py <baseline> <verify>`
**Source:** `src/cryptanalysis/gaudry_cubic.rs` as of commit `55fd0dfa`; both
arms run the same binary, the switch is the environment variable.
**Registered in advance:** §11.12 (falsifier); §11.13 (ceiling — corrected below).
**Suite:** `AGENTS.md` §8's frozen WDSat regression prices a SAT-solver stage on
binary-field encodings and does not apply to this `F_p` Macaulay solve; the
matched suite is §11.5's protocol — same sizes, seeds and `--cross-check` —
run baseline and candidate, full pipeline, cold.

§11.13 declined to write the solver on the strength of a `1.06×` ceiling.  It is
now written — `GAUDRY_LAZY_ELIM=1`, off by default — and writing it showed the
ceiling was wrong by a factor of thirty.

**What was built.**  The solve keeps the Macaulay rows exactly as constructed
and reduces a row only when a normal form first reads it: `NormalForms` asks for
the pivot row of a column, and that row is cleared left of its leading column —
finding and reducing, first, each pivot row it has to be cleared with — then
normalised.  It is the forward elimination taken row by row instead of column
by column, so a row that is read costs exactly what it costs there and a row
nothing reads is never touched.  Which row leads which column comes from a plan
learned once per curve from the first full solve.  Nothing in the plan is
trusted: a planned row whose entry cancels at its column (an accident of the
values, about `1/p` per entry) is replaced by another row that leads there, and
a row that leads at a column the plan does not pivot proves this residual's
staircase is not the plan's, which sends it back to the full elimination with
the lazy work still charged.  `GAUDRY_LAZY_VERIFY=1` re-runs every solve through
the full elimination, uncharged, and counts any `M_{e₁}` entry or staircase that
differs, and any fallback on a residual whose staircase was the plan's after
all.  It is a diagnostic, not a guard: it counts, and the lazy answer is used.

**1. The ceiling was `1.002×`, not `1.06×` — class `accounting`.**  §11.13
charged each elimination multiplication to the pivot row doing the subtracting
and counted as skippable the work of every pivot the normal forms never read.
But a pivot nobody reads can still be one that was *subtracted from* a row
somebody does read, and that subtraction is what puts the read row in echelon
form; skip it and the row is wrong.  What can be skipped is the work done *on*
rows outside the **closure**: the pivot rows the normal forms read and,
transitively, every pivot row subtracted from one of those.  `echelon_traced`
now records which original row went where at what cost, so both attributions
are priced on the same residuals (`GAUDRY_DEBUG_SOLVE=1`, seed 1,
`--cross-check`; the first 300 residuals above `p = 271`):

| `p` | residuals | pivot rows read | closure | work outside the closure | §11.13's attribution |
|---:|---:|---:|---:|---:|---:|
| 271 | 727 | 148–151 of 222 | **201–202** of 222 | 1,648 muls, **1.10 %** of elimination | 48,929, 32.6 % |
| 523 | 300 | 146–151 | 202 | 1,642, 1.09 % | 49,107, 32.6 % |
| 1039 | 300 | 148–149 | 201–202 | 1,655, 1.10 % | 49,213, 32.6 % |
| 2083 | 300 | 149 | 202 | 1,643, 1.09 % | 49,260, 32.6 % |

The normal forms read two thirds of the pivot rows, as §11.13 said, but those
rows were reduced through almost all the rest: the closure is **202 of 222** at
every size.  Elimination is `17.1 %` of `C₃` (`150,099` of `875,934`), so the
skippable part is `0.19 %` of `C₃` and the ceiling is
**`1/(1 − 0.0019) = 1.002×`**.  §11.12 had registered that "the realistic
expectation is `≈ 1×` or worse, because reducing the border monomials still
needs enough of the row space to do it"; that sentence was right, and §11.13's
revision of it was not.

**2. The staircase is not the same on every residual.**  §11.12 and §11.13 saw
one staircase on 20 and 25 residuals.  Over whole runs it is the same on
`97.6–98.1 %` of residuals at `p = 271`, rising to `99.6 %` at `p = 2083`; the
rest — the `border_unreachable` ones among them — have a different order ideal
(`[0,0,8]` leaves it and the dimension drops to 61, or `[3,0,3]` leaves and
`[0,4,2]` joins).  No per-curve plan serves those, and they are exactly the lazy
solver's fallbacks: the verification mode finds **zero** fallbacks on a residual
whose staircase was the plan's.  The fraction falls roughly as `1/p`.  So
§11.12's "Mourrain's iteration is not needed" holds for the generic residual
and not for these.

**3. The measurement.**  Paired cell by cell.  `summarize_lazy_elim.py` first
asserts that residuals, decompositions, independent relations, the planted
logarithm, and the charpoly, eigenvector and root counters — none of which the
elimination touches — are identical in both arms, so the residual streams did
not diverge and every difference below is exact, not statistical.  The equal
charpoly and root counters prove more than that: a lazy attempt that got as far
as `M_{e₁}` and then fell back would have paid for a second characteristic
polynomial, so every one of the 128 fallbacks happened before the random stream
was touched.  `C₃` is
everything the solve spends per residual, including each lazy attempt that
fell back and the full elimination that followed; the speedup is `AGENTS.md`
§8's `baseline_total_operations / candidate_total_operations`, whole method,
cold.

| `p` | `log₂ n` | baseline `C₃` | lazy `C₃` | ratio | elimination / residual, base → lazy | lazy solves (all verified) | fallbacks | `M_{e₁}` mismatches | cross-check | **speedup** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 271 | 24.2 | 875,934 | 876,618 | 1.0008 | 150,099 → 150,798 | 713 | 14 | 0 | 0 / 728 | **0.9992** |
| 271 | 24.2 | 878,665 | 879,905 | 1.0014 | 150,150 → 151,476 | 699 | 17 | 0 | 0 / 719 | **0.9986** |
| 523 | 27.1 | 890,109 | 889,649 | 0.9995 | 150,623 → 150,161 | 1,587 | 15 | 0 | 0 / 1,604 | **1.0005** |
| 523 | 27.1 | 888,956 | 888,223 | 0.9992 | 150,621 → 149,931 | 1,430 | 12 | 0 | 0 / 1,444 | **1.0008** |
| 1039 | 30.1 | 904,961 | 903,920 | 0.9989 | 150,941 → 149,910 | 2,411 | 14 | 0 | 0 / 2,428 | **1.0011** |
| 1039 | 30.1 | 909,048 | 908,149 | 0.9990 | 150,967 → 150,072 | 2,625 | 18 | 0 | 0 / 2,645 | **1.0009** |
| 2083 | 33.1 | 919,672 | 918,394 | 0.9986 | 151,096 → 149,830 | 5,257 | 18 | 0 | 0 / 5,276 | **1.0012** |
| 2083 | 33.1 | 920,087 | 918,826 | 0.9986 | 151,104 → 149,851 | 5,276 | 20 | 0 | 0 / 5,302 | **1.0012** |

**Pooled speedup `1.0009`.**  Every run recovered its planted logarithm;
`S / rho` moves by the same factor in every cell, below the precision the
scoreboard draws it at.

On the 713 residuals it solves at `p = 271`, seed 1, the lazy elimination spends
`148,450` multiplications where the full one spends `150,080` on the same
residuals: `1,630` saved of the `1,648` their closure allows — **`98.9 %` of the
ceiling**.  Each fallback costs
the lazy attempt that found the staircase off the plan, about `120,000`, on top
of the full elimination.  The two cross at a fallback rate of about `1.3 %`:
above it, at `p = 271` (`1.9–2.4 %`), the lazy solver loses; below it, from
`p = 523` on, it wins — by `0.14 %` at `p = 2083` against the `0.18 %` the
closure allows there.

**Against §11.12's falsifier.**  Zero cross-check mismatches in all eight cells
(20,146 residuals per arm) and zero `M_{e₁}` mismatches over 19,998 verified
lazy solves, so the disqualifying condition is not met.  The cost condition was
written as "`C₃` below §11.6's `0.88 × 10⁶`".  Read literally it passes at
`p = 271` and fails from `p = 523` on — but only because the baseline itself
rises from `0.876` to `0.920 × 10⁶` with `p` — `28,000` of the `44,000` is root
finding, which grows with `log p` — and that says nothing about the solver.  The registration meant the matched
baseline, and against that: **failure at `p = 271` on both seeds, the
registered `engineering` success at every larger size on both seeds**, `0.14 %`
at best.  That is the size §11.12 said in advance it would be.  The **coverage**
outcome is none: `border_unreachable` and the unsolved count are identical in
the two arms, because the residuals the fixed-degree cut cannot close are among
the ones whose staircase no plan fits.  Whether Mourrain's full iteration would
close them is not answered by this build; they are 0–7 residuals per run of
thousands, and the meet-in-the-middle oracle solves them correctly.

**Two bugs the guards caught.**  The first version used the plan's row for
every column and fell back on `53` of `199` test residuals at `p = 271` —
consistent with a `1/p` chance of cancellation at each of some eighty exposed
pivot entries — which is what the row substitution is for.  A later version
passed `unwrap_or` an argument with a side effect, so a row adopted as the pivot
of an earlier column was filed under the wrong one.  In these runs both showed
up as fallbacks rather than wrong answers — a row cleared with a misfiled pivot
leaves the row space and leads at a column the plan calls standard, which the
solver reads as an off-plan staircase and refuses — and an uncharged full
elimination on each fallback, now the verification mode's needless-fallback
count, found the one "off-plan" residual whose staircase was in fact the plan's.
That, and `the_lazy_elimination_reads_the_same_normal_forms` (200 residuals in
lockstep against the full solve: identical answers, identical `M_{e₁}`, no
needless fallback), is the evidence the solver is right; the cross-check column
is the evidence its answers are.

**Class.**  The ceiling correction is `accounting`: nothing ran faster, a number
was wrong.  The solver is `engineering` — `−0.14 %` to `+0.14 %` of `C₃`,
`1.0009×` pooled, on a phase §11.5 measures as flat in `n`, in a method `1,989×`
from rho.  It does not move the verdict and §11.12 registered that it could not.
It stays in the tree behind `GAUDRY_LAZY_ELIM`, off by default so every earlier
run reproduces, because the verified negative is the result and the closure
diagnostic it brought is the instrument that would have caught §11.13's mistake
before it was published.

**What this closes.**  With the Macaulay cut, the order of the rows, the row
selection, the retries and now the elimination itself each measured, what is
left in `C₃` is the normal forms (`40 %`), which a border basis computes rather
than avoids, and the characteristic polynomial and eigen-solve (`41 %`), which
do not care how `M_{e₁}` was reached.  There is no open lever on `C₃` left in
this design, and none of the ones measured was worth more than its registration
said.

### 11.15 Close-out: what §11 established, and what would count now

§11.14 closed the last lever this section had listed.  This is the ledger of
all of them, in the note's unit, and the one question left open by the
arithmetic rather than by a missing build.

**The ledger.**  Every row was measured end to end and every run recovered its
planted logarithm.  The class column is the one each section recorded; §11.4
and §11.7 predate the practice and recorded none, so those rows say what
moved instead of being classed after the fact.

| § | lever | what moved | class, as recorded | result |
|---|---|---|---|---|
| 11.4 | `O(1)` Gröbner solve replaces the loop over the base | relation phase `n^{0.69} → n^{1/3}` | — | `C₃ = 4.82·10⁶`, flat in `n`; dearer than the loop below `n ≈ 2^{34}` |
| 11.5 | honest multiplication count | the count, not the algorithm | accounting | `C₃ = 3.23·10⁶` |
| 11.5 | echelon form with back-substitution | `C₃` | engineering | `1.53·10⁶` |
| 11.6 | Macaulay's row selection | `C₃` | engineering | `0.879·10⁶` |
| 11.8 | no retries once the border is unreachable | `C₃` | engineering | `0.876·10⁶` |
| 11.7 | Wiedemann with filtering | linear algebra `n^{0.85} → n^{0.68}` | engineering | plain method bottoms at `≈ 200×` rho near `2^{50}`, then rises |
| 11.7 | double large primes, `\|F'\| = \|F\|^{2/3}` | total exponent `0.32 → 0.44`, *below* rho's `1/2` | — | `1,989×` rho at `2^{33.1}` (`3.8×` the plain method's work), closing as `n^{-1/18}` |
| 11.9 | Joux–Vitse: decompose into `k − 1 = 2` points | `C₃ → 1,513`, `580×`; residuals `n^{1/3} → n^{2/3}` | advance below `2^{28.5}`, relabelling above | `1,037×` rho at `2^{33.1}`, exponent worse |
| 11.11 | merge-level cap on the large-prime eliminator | solve exponent `−0.045 ± 0.017`, merge `+0.043 ± 0.030` | engineering, negative | end to end `+0.051 ± 0.031` **worse** |
| 11.14 | border basis (demand-driven elimination) | `C₃` by at most `0.19 %` | accounting (ceiling) + engineering | `1.0009×` pooled; off by default |

**The verdict, restated with nothing left to try at `k = 3`.**  On
`E(F_{p³})` the best variant measured is the plain method at `528×` rho
(`2^{33.1}`).  Its linear algebra grows as `n^{0.68}`, faster than rho's
`n^{1/2}`, so it bottoms near `200×` rho around `2^{50}` and loses ground
after.  The double-large-prime variant is the only one with an exponent below
rho's — `n^{0.424 ± 0.016}` end to end over six seeds, against the theorem's
`4/9` — and it closes the measured `1,989×` at `n^{-1/18}`: past `2^{230}`,
extrapolated on those two exponents.  The four constant levers since §11.4
bought `5.5×` on `C₃`, and §11.12 showed before the last of them that no `C₃`
lever can move an exponent, because `C₃` is flat in `n`.

**What would count now is a change in the closing rate, not in a constant.**
Two exponents set it: the relation phase's, and the linear algebra's over
`|F| ∝ n^{1/k}` unknowns.  At `k = 3` the linear algebra is the obstacle —
`|F|² ∝ n^{2/3}` for any sparse solver — and the large-prime cure trades it
for a residual count that closes only at `n^{-1/18}`.  **The one structural
lever left is the extension degree.**  Derived, not measured: at `k = 4`,
`|F| ∝ n^{1/4}`, so the linear algebra is `|F|² ∝ n^{1/2}` — *the same
exponent as rho* — and the relation phase is `k!·|F|·C₄ ∝ n^{1/4}`.  The plain
method's `S / rho` then tends to a constant instead of growing without bound
as it does at `k = 3`.  That constant is `c_LA / c_rho`, the sparse solver's
cost against rho's, and it does **not** depend on the solve: it follows from
§11.7's measured Wiedemann counters and the cost of `F_{p⁴}` arithmetic, and
§11.16 derives it.  `C₄`, the cost of one `S₅` solve over `F_{p⁴}`, which
nothing in this module has measured, sets only *where* the relation phase
stops dominating.

**What does not count, and is not worth building here:**

- any lever on `C₃`, `C₄` or the linear-algebra constant at fixed `k`: flat in
  `n`, `engineering` by construction (§11.12);
- Mourrain's full border-basis iteration for the residuals whose staircase is
  not the generic one: `0–7` residuals per run of thousands, already solved
  correctly by the fallback oracle, and unable to move cost (§11.14);
- a larger merge cap, a different large-prime budget, or a different small
  base at `k = 3`: §11.11 measured the family and the method got worse.

### 11.16 Pre-registration: `k = 4`, derived before anything is built

**Written before any `F_{p⁴}` code exists.**  §11.15 named the extension degree
as the one lever that changes the closing rate.  Most of what decides the
`k = 4` question can be derived from counters this note has already measured,
so it is derived here first — as §11.12 did for the border basis — and the
build is registered only for the part the derivation cannot reach.
Everything below is derived or extrapolated, and says which.

**The setting.**  `E` over `F_{p⁴}`, base `F = {P : x(P) ∈ F_p}` with
`|F| ≈ p/2`, group order `n ≈ p⁴`.  A residual decomposes as a signed sum of
four base points; the Weil restriction of the symmetrised `S₅` gives four
equations over `F_p` in `e₁, …, e₄`, each of total degree `≤ 8`, with
`8⁴ = 2^{k(k−1)} = 4,096` solutions — against `4³ = 64` at `k = 3`, which §11.4
measured the solver attaining.

**Inputs, and where each comes from.**

| input | value | source |
|---|---:|---|
| residuals per relation | `k! = 24` | the count of §11.3; `≈ 6` at `k = 3` measured `0.15–0.18` against `1/6` |
| `F_p` multiplications per `F_{p⁴}` multiplication | `19` | schoolbook with `t⁴ = c`, the convention behind §11.1's `11` at `k = 3` |
| per `F_{p⁴}` inversion | `40` | norm to `F_{p²}`: `24` multiplications plus the `F_p` inversion, charged `16` as `Fp3::inv`'s count of `30` implies (it performs `14`) |
| per affine addition, `c_add` | **`97`** | one inversion and three multiplications, the decomposition that gives exactly `63` at `k = 3` |
| per multiplication mod `n` | `16` | schoolbook limbs, as §11.7 charges `9` for `n ≈ p³` |
| Wiedemann, per unknown² | `(5 + 3w) = 20` at row weight `w = k + 1 = 5` | §11.7 measured `17.0·N²` at `w = 4` at all four sizes; the split `2N` dot products, `3wN²` mat-vecs, `≈ 2N²` Berlekamp–Massey reproduces the `17` |
| unknowns after filtering, `φ` | `0.73–1.0` of `\|F\|` | `0.71–0.76` measured at `k = 3`; `1.0` if filtering finds nothing to drop |
| rho | `S ≈ 1.3` | the reference `AGENTS.md` fixes for this thread's unit |

**1. The asymptotic constant: `r∞ = c_LA / c_rho ≈ 0.34–0.63`.**  With
`N = φp/2` unknowns, the linear algebra costs `20·16·(φp/2)² = 80φ² p²` `F_p`
multiplications and rho `1.3 · p² · 97 = 126 p²`, both `∝ n^{1/2}`.  Their ratio
is `r∞ = 0.634 φ²`: **`0.34`** at `k = 3`'s filtering rate, `0.63` with none.
Below one either way, so the plain `k = 4` method, unlike `k = 3`'s, does have
a regime where it beats rho — by a factor between `1.6` and `3`, never more,
because the linear algebra it would converge to is itself `∝ n^{1/2}`.  This
number does not depend on the solve at all.

**2. Where the relation phase hands over.**  `24|F| = 12p` residuals at `C₄`
each put the relation phase at `12p·C₄ / 126p² = 0.095·C₄/p` of rho, so

```text
S / rho  =  0.095 · C₄ / p  +  r∞        crossover:  p* = 0.095 · C₄ / (1 − r∞)
```

**3. A floor on `C₄`, from this note's own solver.**  §11.4–11.6's design reads
the eigenvalues off the characteristic polynomial of a `D × D` multiplication
matrix.  That step alone measured `272,392 = 1.039·64³` per solve at `k = 3`
(`experiments/23_gaudry_lazy_elim_baseline.json`, `p = 271`, seed 1).  At
`D = 4,096` it is **`C₄ ≥ 7.1·10¹⁰`** before any Macaulay elimination — whose
matrix at the regularity degree `4·7 + 1 = 29` has `C(33, 4) = 40,920` columns
against `k = 3`'s `286` — or any normal form.  At `k = 3` the characteristic
polynomial was `31 %` of `C₃`.

**4. The crossover the floor implies — extrapolated on the exponents
`n^{1/4}` (relations) and `n^{1/2}` (linear algebra and rho).**

| `C₄` | `r∞ = 0.34` (`φ = 0.73`) | `r∞ = 0.63` (`φ = 1`) |
|---|---:|---:|
| floor, `7.1·10¹⁰` | **`n* ≈ 2^{133}`** | `2^{136.5}` |
| `3×` the floor | `2^{139}` | `2^{143}` |
| `10×` the floor | `2^{146}` | `2^{150}` |

At the smallest size this note runs, `p = 271` (`n ≈ 2^{32}`), the formula puts
`k = 4` at **`2.5·10⁷×` rho**, against `k = 3`'s `528×` at `2^{33}`.  Past
`n*` it would sit at `0.34–0.63×` rho; before it, the relation phase decides
everything.  Against `k = 3`'s double-large-prime crossover past `2^{230}`,
`k = 4` hands over about a hundred doublings sooner.  That is a statement
about exponents already derived, not a measurement.

**What a measurement can and cannot change.**  A measured `C₄` in this design
cannot fall below the floor, so it can only move `n*` **later** — by `4 log₂ m`
bits for a `C₄` that is `m×` the floor.  It cannot touch `r∞`.  What is not in
this design is not bounded by it: an `F₄` / sparse-FGLM solver, or the
symmetries Faugère, Gaudry, Huot and Renault use on curves with rational
torsion, lower `D` or the `D³`, and in this model a factor `f` off `C₄` moves
`n*` about `4 log₂ f` bits earlier.  Those are levers on a constant, and would
be registered as `engineering` if built.

**Two variants the derivation settles without a build.**

- **Joux–Vitse, decompositions into `k − 1 = 3` points.**  A residual is a
  three-point sum about `1/(6p)` of the time, so `≈ 3p²` residuals —
  `∝ n^{1/2}`, rho's exponent again — at `C′` each, where `C′` is the
  overdetermined `S₄` solve.  `S / rho` tends to `3C′/126 + r∞`: below one only
  if `C′ < 28` `F_p` multiplications (`15` if `φ = 1`).  `k = 3`'s two-point
  solve already costs `1,513`, so at any `C′` of that order this variant sits
  near **`36×` rho at every size**.  The solve cost enters its asymptote; it
  does not enter the full-decomposition variant's.
- **Double large primes.**  They trade the `p²` linear algebra for more
  residuals, `Õ(p^{3/2})` of them.  The extra residuals cost a factor
  `∝ p^{1/2}` of `C₄` each and save linear algebra `∝ p²`, so the variant beats
  the plain method only once `p^{1/2}` outgrows `C₄` up to constants — `p` of
  order `C₄²`, far past `n*`.  It is the asymptote (`n^{3/8}`, closing on rho
  as `n^{-1/8}`), not the crossover.

**The build this leaves.**  Only `C₄` — the one measured input the crossover
depends on, and the check that the `S₅` system behaves generically
(`D = 4,096`, regularity `29`).  Registered target: `C₄` per residual, for a
full solve in §11.4–11.6's design, at `p = 271` on at least three residuals,
with **every** output cross-checked against a meet-in-the-middle decomposition
oracle (zero mismatches), counted under this note's accounting.  Prediction:
`C₄ ≥ 7.1·10¹⁰`, so `n* ≥ 2^{133}`.  Inadmissible: a special curve family or
symmetry without registering it as a separate lever; changing the
multiplication charges above; quoting `n*` as anything but an extrapolation.
Abandoned, with whatever bound was reached recorded, if the degree-`29` solve
cannot run in this container (`15 GB`; a dense `36,824 × 40,920` matrix is
`6 GB` in `u32` and `12 GB` in `u64`) or takes more than a few hours a
residual.

**What it costs.**  None of the `k = 3` code carries over directly: `Fp3`,
`E3`, the symbolic `S₄` and the solver are all hard-wired to three.  The build
is `F_{p⁴}` arithmetic and a curve generator; `S₅` via
`Res_Y(S₄(x₁, x₂, x₃, Y), S₃(x₄, x₅, Y))`, symmetrised into `e₁, …, e₄`, which
is `≤ 495 · 9` coefficients; its Weil restriction; a four-variable Macaulay
solve at degree `29` with row selection, normal forms and a `4,096`-dimensional
eigen-solve; and a meet-in-the-middle oracle to check it.  That is days of
work, and runs of minutes to hours per residual.  Its payoff is one number
that can only confirm the floor or move `n*` later.

### 11.17 `C₄`, measured: `17×` the floor, and the crossover moves to `2^{150}`

**Runner:** `cargo run --release --example gaudry_quartic_c4 -- --p 269 --seed 1
--residuals 4 --json …`, and `-- --generic 2,3,4,5,6,7` for the scaling panel.
**Frozen:** `experiments/24_gaudry_quartic_c4.json` (and `.log`, with the phase
progress), `experiments/24_gaudry_quartic_generic.json` (and `.log`).
**Summary:** `python3 scripts/summarize_quartic_c4.py <c4> <generic>`
**Source:** `src/cryptanalysis/gaudry_quartic.rs`, `examples/gaudry_quartic_c4.rs`.
**Registered in advance:** §11.16, pushed (`4d4f4cb8`) before any `F_{p⁴}` code
existed (`220b8883`).

**One deviation from the registration.**  §11.16 named `p = 271`.  For
`p ≡ 3 (mod 4)` every binomial `t⁴ − c` is reducible, so `F_p[t]/(t⁴ − c)` is not
a field there; the run uses `p = 269`, the nearest prime `≡ 1 (mod 4)`.  `C₄`
depends on `p` only through root finding, `10⁻⁶` of it.

**What was built** (`gaudry_quartic`): `F_{p⁴}` under `Fp3`'s counting
convention; a curve with coefficients outside `F_{p²}`; `S₅` evaluated
numerically as `Res_Y(S₃(x₁, x₂, Y), S₄(x₃, x₄, x₅, Y))`; the symmetrised `S₅` by
interpolation — `495` `e`-monomials by `9` powers of `x_R`, solved once per curve
(`1.2·10⁹` multiplications, outside `C₄` as `k = 3`'s precomputation is) and
checked against fresh resultant evaluations; and the §11.4–11.6 solver in four
unknowns — Macaulay matrix at degree `29` with Macaulay's row selection, forward
elimination by leading column, normal forms memoised right to left, `M_{e₁}`,
and the characteristic polynomial, roots and eigenvectors through the same
counted routines as `k = 3`.  A meet-in-the-middle oracle over pairs of base
points checks every decomposition.  Tests: `S₅` vanishes on four-point sums and
is symmetric in its first four arguments; the solver recovers planted roots of
random quadrics and cubics in four unknowns at the Bézout dimension `d⁴`.

**1. The inputs §11.16 derived hold.**  `F_p` multiplications per `F_{p⁴}`
addition measure **`97.0`**, the derived figure.  The quotient has dimension
**`4,096`** on every residual at the Macaulay bound `29`: the `S₅` system is as
generic as §11.16 assumed.  The matrix is `41,780 × 40,920` (`3.4 GB` as `u16`),
rank `36,824`.

**2. `C₄`.**

| residual | `C₄` | elimination | normal forms | charpoly | eigenvectors | rational eigenvalues | solver = oracle | wall |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| 0, constructed | `1.227·10¹²` | `5.825·10¹¹` | `5.461·10¹¹` | `7.06·10¹⁰` | `2.77·10¹⁰` | 2 | `{41, 122, 180, 226}` = planted | `1,625 s` |
| 1 | `1.213·10¹²` | `5.825·10¹¹` | `5.460·10¹¹` | `7.06·10¹⁰` | `1.39·10¹⁰` | 1 | `∅ = ∅` | `1,460 s` |
| 2 | `1.213·10¹²` | `5.825·10¹¹` | `5.461·10¹¹` | `7.06·10¹⁰` | `1.39·10¹⁰` | 1 | `∅ = ∅` | `1,499 s` |
| 3 | `1.199·10¹²` | `5.825·10¹¹` | `5.460·10¹¹` | `7.06·10¹⁰` | `0` | 0 | `∅ = ∅` | `1,450 s` |

**`C₄ = 1.213·10¹²`, `17.0×` the registered floor.**  The `2.3 %` spread is the
eigenvector step alone, which scales with the number of rational eigenvalues;
elimination, normal forms and the characteristic polynomial agree to four
figures on every residual.  The cross-check agrees on all four, but three of
those agreements are on the empty set — a residual decomposes about one time in
`24` — so the constructed residual, whose planted decomposition both the solver
and the oracle return, is the check that carries weight.

Elimination is `48 %` of `C₄` and the normal forms `45 %`, each alone about `8×`
the floor.  The characteristic polynomial costs `1.027·D³`, a shade under the
`1.039` §11.16 carried over from `k = 3`: the generic runs show the coefficient
drifting down with `D` (`1.063` at `81`, `1.032` at `2,401`), so the floor's own
component was `1.2 %` high.  It is `6 %` of `C₄`.

**3. The same solver on generic systems predicts it.**  Four random equations of
degree `d` in four unknowns with a planted root, every root found, every
quotient at the Bézout dimension:

| `d` | columns | `D = d⁴` | `C` | elimination / normal forms | charpoly / `D³` | wall |
|---:|---:|---:|---:|---|---:|---:|
| 2 | 126 | 16 | `4.6·10⁴` | 44 % / 31 % | 1.010 | — |
| 3 | 715 | 81 | `6.2·10⁶` | 47 % / 41 % | 1.063 | — |
| 4 | 2,380 | 256 | `2.3·10⁸` | 46 % / 42 % | 1.054 | `1 s` |
| 5 | 5,985 | 625 | `3.7·10⁹` | 47 % / 44 % | 1.045 | `5 s` |
| 6 | 12,650 | 1,296 | `3.5·10¹⁰` | 47 % / 44 % | 1.038 | `58 s` |
| 7 | 23,751 | 2,401 | `2.4·10¹¹` | 47 % / 44 % | 1.032 | `419 s` |

`C ∝ D^{3.11}` over `d = 4 … 7` predicts `1.26·10¹²` at `D = 4,096`, **`4 %`**
above the `S₅` measurement.  Nothing about the summation polynomial makes this
system cheaper or dearer than a generic one of its degrees: `C₄` is what a
Macaulay solve of four octics in four unknowns costs, and `C_k` in this design
is `≈ D^{3.1}` with `D = 2^{k(k−1)}`.

**4. The crossover, extrapolated** on §11.16's formula
`p* = 0.095·C₄ / (1 − r∞)` and the exponents `n^{1/4}` and `n^{1/2}`:

| `C₄` | `r∞ = 0.34` | `r∞ = 0.63` |
|---|---:|---:|
| registered floor, `7.1·10¹⁰` | ~~`2^{133}`~~ | ~~`2^{136.5}`~~ |
| **measured, `1.21·10¹²`** | **`n* ≈ 2^{149.4}`** | **`2^{152.8}`** |

`17×` the floor moves `n*` `4 log₂ 17 ≈ 16` bits later — the only direction the
registration allowed.  At `p = 269` (`n ≈ 2^{32}`) the relation phase alone
would be **`4.3·10⁸×` rho**.  Against `k = 3`'s double-large-prime crossover
past `2^{230}`, the full-decomposition `k = 4` method still hands over about
**eighty doublings sooner**, and past `n*` it would sit at `0.34–0.63×` rho
rather than tend to zero.

**Class.**  A stage diagnostic: `C₄` prices one phase, and `n*` extrapolates the
method from it on derived exponents.  No end-to-end `S` exists for `k = 4`, and
none is claimed.  Nothing here changes an existing method's cost, so none of
§3's four classes applies.  §11.16's prediction — `C₄ ≥ 7.1·10¹⁰`,
`n* ≥ 2^{133}` — is met, and its abandonment condition was not reached: each
solve took under half an hour in `3.7 GB`.

**What this leaves.**  At `k = 4` the crossover is a question about one
constant.  Every factor `f` taken off `C₄` moves `n*` about `4 log₂ f` bits
earlier, and in this design `C₄ ≈ D^{3.1}`, so the levers are the ones that
shrink `D` or its exponent — `F₄` with sparse FGLM, or the Faugère–Gaudry–Huot–
Renault symmetries on curves with rational torsion.  Any of them would be
registered as `engineering`: none changes `r∞`, which caps a plain `k = 4`
method at `1.6–3×` better than rho however cheap the solve gets.

### 11.18 Pre-registration: measuring `r∞`, the one claim that puts `k = 4` below rho

**Written before the measurement exists.**  §11.16 derived that a plain `k = 4`
method tends to `r∞ = c_LA / c_rho = 0.34–0.63×` rho, and §11.17 carried it
into the scoreboard as "would eventually beat rho by `1.6–3×`".  It is the only
statement in §11 that puts anything below rho, and every input to it is
borrowed or derived: the Wiedemann constant from `k = 3` rescaled to row
weight `5`, the filtering fraction `φ` from `k = 3`, the decomposition rate from
a count, rho's `S` from the thread's convention.  None of them needs the `S₅`
solve, so all of them can be measured directly, at sizes where the method
runs end to end.

**What runs.**  Prime-order curves over `F_{p⁴}` for `p ∈ {269, 521, 769, 1033}`
(`n ≈ 2^{32.3}, 2^{36.1}, 2^{38.3}, 2^{40.0}`), two seeds each.  Relations come
from the meet-in-the-middle oracle of §11.17 — **not** the `S₅` solve: the
relation phase is not what is measured here, and its cost is reported beside
the table and kept out of the ratio.  Relations are filtered to a square core
and solved by the §11.7 Wiedemann, with every attempt's operations counted;
the logarithm is accepted only if `[d]G = Q`.  Rho runs `16` times per curve on
the same group, because a single rho run's `S` moves by `2×` between seeds.

**The statistic.**  On each curve, both measured in `F_p` multiplications:

```text
r  =  LA  /  rho  =  (la_ops · 16)  /  (S_rho · √n · c_add)
```

with multiplications mod `n` charged `16` and `c_add` the measured `97.0`, as
§11.16 did.  Both terms grow as `n^{1/2}` if §11.16 is right, so `r` should be
flat across the four sizes; its value is `r∞`.  Reported beside it, each input
§11.16 derived: the decomposition rate against `1/24`, `φ` against `0.73`, the
Wiedemann constant `la_ops / N²` against `20`, the row weight against `5`, rho's
`S` against `1.3`, and the linear algebra's exponent in `n` against `1/2`.

**Registered outcomes.**

| outcome | condition |
|---|---|
| **confirmed** | `r` flat within its seed spread and inside `0.34–0.63` |
| **corrected** | `r` flat and below `1`, but outside `0.34–0.63`: the cap becomes `1/r` |
| **withdrawn** | `r ≥ 1`, or `r` rising with `n`: a plain `k = 4` method never beats rho, and §11.16–11.17's "`1.6–3×`" comes off the scoreboard |

**Inadmissible:** changing the `16` charge or `c_add` after seeing the data;
dropping filtering or failed Wiedemann attempts from `la_ops`; excluding a
seed; using fewer than `16` rho runs for any curve; any run without a verified
logarithm.  An `S₅`-based relation phase would change nothing here, since
`r∞` does not depend on the solve.

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
