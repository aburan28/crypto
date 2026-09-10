# Finding points versus finding relations: residual walks over partial decompositions

**Module:** `src/cryptanalysis/residual_walk.rs`
**Bench:**  `cargo run --release --example residual_walk_bench -- --panel --json experiments/20_residual_walk_panel.json`
**Data:**   `experiments/20_residual_walk_panel.json`, `experiments/20_residual_walk_panel.log`
**Tables:** `python3 scripts/summarize_residual_walk_panel.py experiments/20_residual_walk_panel.json`

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
  literal `s_{t+1} = H(L(s_t))` construction.

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
```

The bench's single-instance mode also accepts `--strategies A,C1,R`,
`--dp BITS` (applied to the collision-preserving strategies only),
`--k`, `--budget OPS` and `--json FILE`.

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
