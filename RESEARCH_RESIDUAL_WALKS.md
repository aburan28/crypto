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

<!-- PANEL TABLES -->

## 5. Reading the numbers

**The count is the invariant.**  In every P1 row, A, B, C1 and C2 need
`samples/pred ≈ 1` — within the seed-to-seed spread, the same
`√(2n(B+1))` residuals whether the residual was drawn independently,
mutated locally, or produced by a collision-preserving walk.  Doubling
`B` four-fold (64 → 256) raises the count by about `2×`, as
`√(B+1)` says.  The walks do not find useful collisions faster than
independent sampling; they only find them *cheaper per residual*.

**Cheaper per residual is the whole gain, and it is bounded.**  B and C1
bring the per-residual cost from `≈ 3·log₂n` down to `≈ 2–6`
operations, which is why they finish `10–30×` faster than A at the same
count.  That gain saturates: one group operation per residual is the
floor, plain rho already sits on it, and rho needs `√(2(B+1))` times
fewer residuals.  Measured against rho's own total on the same curve
(the `×R` column, which includes rho's setup and replay), the best
hybrid is still tens of times slower at `B = 64` and worse as `B` grows.

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
however many relations it collects, and it solves `d` as soon as it
reaches that rank.  A collision-preserving walk with cheap incremental
updates can only ever generate relations inside the span of its update
vectors; to reach the full base it needs `r ≥ B` distinct updates, at
which point the coefficient vectors are dense and must be replayed (or
stored) — the `replay` column in P1.

**The filter trap (P4).**  Rejecting residuals outside `x < p/2^s`
multiplies the total work by `≈ 2^{s/2}` — the measured ratios track
`√(p/M)` — while the accepted-residual count stays on the
`√(2M(B+1))` birthday line.  Filtering for special-looking residuals
after the fact buys nothing, exactly as the model says.

**Distinguished points (P5).**  For C1, C2 and R, storing only
residuals with `dp` zero hash bits cuts the table by `≈ 2^{dp}` at a
cost of a few extra steps per walk; the relation count and the solve
are unchanged.  For A and B the module refuses the option, because a
residual collision between non-preserving states does not propagate
to a later distinguished point.

**Meet in the middle (P6).**  At `B ≈ (4n)^{1/4}` the 4-decomposition
search delivers about one weight-4 relation per target at `≈ B²/2`
operations each; `B + 1` such relations recover `d`, at a total that
grows like `n^{3/4}` rather than rho's `n^{1/2}`.  Doubling `B` finds
`≈ 16×` more matches per target at `4×` the cost per target, so the
per-relation cost falls, but the number of relations needed rises with
`B` and the total does not beat rho at any tested size.

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
