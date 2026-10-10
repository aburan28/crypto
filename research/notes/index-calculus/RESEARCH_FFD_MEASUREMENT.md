# First-Fall-Degree measurement on binary Semaev systems

**Module:** `src/cryptanalysis/ffd_harness.rs`
**Demo:**  `cargo run --release --example ffd_sweep_demo`
**Provenance:** Direct response to the "FFD is controversial / nobody
understands why some systems collapse" line in the ECDLP research-
direction map, and to Galbraith / Petit's 2015 calls for **empirical**
fall-degree numbers on Semaev systems.

## What the harness does

For each `n ∈ n_range`:

1. Pick random `b ∈ F_{2^n}` and `x_3 ∈ F_{2^n}`.
2. Form `S_3(X_1, X_2, x_3) = 0` over `F_{2^n}`.
3. **Weil-descend** into `n` quadratic equations over `F_2` in `2n`
   bit-variables (the bits of `X_1` and `X_2`).
4. Build the Macaulay matrix at degree `D ∈ {2, …, d_max}` over `F_2`
   (multilinear monomial basis, with field equations `xᵢ² = xᵢ`
   automatically enforced).
5. Compute the rank by bit-packed Gauss-Jordan over `F_2`.
6. Report `(rows, cols, rank, generic_prediction)` per degree, and
   the **first fall degree** = smallest `D` with `rank < rows` and
   `rank < cols`.

## Result on `n ∈ {3..7}` (seed `0x_FFD_DEAD`)

| `n` | vars | eqs | `D` | rows | cols | rank | gen rank | deficit vs. rows |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 6 | 3 | 2 | 3 | 22 | 3 | 3 | 0 |
| 3 | 6 | 3 | 3 | 21 | 42 | 20 | 21 | **1** |
| 3 | 6 | 3 | 4 | 66 | 57 | 54 | 57 | 12 |
| 4 | 8 | 4 | 2 | 4 | 37 | 4 | 4 | 0 |
| 4 | 8 | 4 | 3 | 36 | 93 | 35 | 36 | **1** |
| 4 | 8 | 4 | 4 | 148 | 163 | 127 | 148 | 21 |
| 5 |10 | 5 | 2 | 5 | 56 | 5 | 5 | 0 |
| 5 |10 | 5 | 3 | 55 |176 | 54 | 55 | **1** |
| 5 |10 | 5 | 4 |280 |386 |255 |280 | 25 |
| 6 |12 | 6 | 3 | 78 |299 | 77 | 78 | **1** |
| 6 |12 | 6 | 4 |474 |794 |441 |474 | 33 |
| 7 |14 | 7 | 3 |105 |470 |104 |105 | **1** |
| 7 |14 | 7 | 4 |742 |1471|700 |742 | 42 |

**Headline:** every system in the sweep has **first fall degree 3** —
exactly one syzygy appears at `D = 3` (rank deficit 1), then a large
collapse at `D = 4`.  This reproduces the Galbraith-Gebregiyorgis
empirical pattern (INDOCRYPT 2014, §4) for `S_3` over small `n`.

## Why it matters for the "FFD controversy"

Petit–Quisquater (ASIACRYPT 2012) assumed the solving degree tracks the
first fall degree, which for Semaev systems gives a subexponential
attack. Huang–Kosters–Yeo (CRYPTO 2015) introduced the last fall degree
and raised doubt on that assumption, and Kosters–Yeo showed the `S_3`
first fall degree is generically 2 for a structural reason (a trace
morphism), far below any solving degree. The FFD *looks* like a small
constant at small `n`, but extrapolation is unsafe. Proved bounds on
both degrees are ledgered in `RESEARCH_FALL_DEGREE_BOUNDS.md`.

This harness lets you **run the experiment yourself** in a few
seconds.  At `n = 7` the Macaulay matrix at `D = 4` already has 742
rows × 1471 columns; pushing `n = 8` is feasible, `n = 9` strains
laptop memory, and `n = 10` is currently out of reach.  That's
exactly the regime the controversy lives in.

## What to do with this output

- **For a researcher tracking FFD growth:** extend `d_max` and `n_max`,
  average over multiple `(b, x_3)` draws (the harness already
  supports it; flip the trial count in `run_sweep`).
- **For a researcher checking a Semaev variant:** swap in a different
  polynomial in `weil_descend_s3` — the rank-deficit measurement is
  agnostic to which quadratic system you feed it.
- **For a teacher illustrating Macaulay matrices:** the
  `(rows, cols, rank, generic_pred)` quartet is the cleanest possible
  illustration of Hilbert-series degeneration.

## Limits

- **`n ≤ 8` for interactive use.**  `n = 9` is ~10 minutes per `D`;
  `n = 10` exhausts dense-matrix memory.  A sparse F4-style solver
  would push to `n ≈ 12`; that's `src/cryptanalysis/groebner_f4.rs`.
- **Quadratic systems only.**  The bit-poly representation is hard-
  coded for degree 2.  Higher-degree systems (e.g. Semaev `S_4`) would
  require a more general monomial-set representation.
- **Single (b, x_3) per `n`.**  For paper-grade numbers, average over
  ~32 random draws and report the mean / max FFD.

## References

- M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil
  descent attacks on ECDLP*, Crypto 2015.
- S. Galbraith, S. Gebregiyorgis, *Summation polynomial algorithms for
  elliptic curves in characteristic two*, INDOCRYPT 2014.
- M. Kosters, S. L. Yeo, *Notes on summation polynomials*,
  arXiv:1503.08001, 2015 — the "backlash" the harness here is
  calibrated against.
- J.-C. Faugère, L. Perret, C. Petit, G. Renault, *Improving the
  complexity of index calculus algorithms in elliptic curves over
  binary fields*, EUROCRYPT 2012.
