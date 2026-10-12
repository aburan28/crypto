# Research Note: Root-Finding Decomposition in place of F4/F5 Gröbner for ECDLP Index Calculus

**Date.** 2026-10-12. **Status.** Measured on public synthetic known-answer
instances; every decomposition re-added in the group, and the two oracles
cross-checked to agree on every target. **No security claim** (§4): the
speed-up is on the decomposition *subroutine* at fixed factor-base size; the
asymptotic group-operation count of index calculus is unchanged. No
deployed-curve, key-recovery, or sub-rho claim is made or implied.
**Code.** `src/cryptanalysis/prime_fast.rs` (prime engine),
`src/cryptanalysis/koblitz_index_calculus.rs` (`DecompositionStrategy::PairTable`,
now the default), `examples/semaev_root_vs_groebner_bench.rs` (this note's
tables), `research/semaev_root_vs_groebner_20261012/` (raw console + JSON).

---

## 0. Summary

The decomposition oracle behind elliptic-curve index calculus asks: *given a
target `R`, is `R = P_1 + … + P_m` with every `P_i` in the factor base?* The
textbook answer Weil-restricts Semaev's summation polynomial to a low-degree
(Boolean, over `F_{2^n}`) system and hands it to a **Gröbner-basis** engine —
matrix-F4, F5. That Gröbner solve is the step index calculus "eventually
needs," and it is the cost that dominates a decomposition.

This note replaces the Gröbner solve with **root finding** and makes it the
default decomposition strategy:

- **Binary curves** — a precomputed table of pair sums `P_i + P_j`
  (`PairSumTable`). A 2-decomposition is one table lookup, a 3-decomposition
  is `|F|` lookups, a 4-decomposition is `|F|²` — the meet-in-the-middle
  half of the birthday search. Exact and complete like exhaustive search,
  but it *finds* the partner instead of row-reducing a Macaulay matrix.
- **Prime curves** — the Semaev `S₃` root in the third variable is literally
  `x(R ∓ F_i)` (`prime_fast::collect_relations`), so a 2-decomposition sweep
  is one batched field inversion plus `|F|` table lookups. No Gröbner basis,
  no per-element square root.

Head to head on the **same** targets:

| regime | instance | root finding | F4 / Gröbner | factor |
|---|---|---:|---:|---:|
| binary (Koblitz `K_0/F_{2^9}`, `|F|=55`, 254 targets, m=2) | per-decomposition | 0.73 ms | 279 ms (matrix-F4) | **383×** |
| binary (`K_0/F_{2^13}`, `|F|=4005`, 256 targets, m=2) | per-decomposition | 2.5 ms | 3.9 s (matrix-F4) | **1535×** |
| prime (16-bit `a=−3`, `|F|=6`, 40 targets, m=2) | per-decomposition | 30 µs | 1.05 s (`groebner_f4`) | **35592×** |
| prime (20-bit, `|F|=6`) | per-decomposition | 33 µs | 1.15 s | 34602× |
| prime (24-bit, `|F|=6`) | per-decomposition | 26 µs | 1.23 s | 46911× |

Every binary target's verdict (decomposable / not) agreed between the
pair-table root finder and matrix-F4; every prime decomposable target's
Gröbner basis was verified to vanish at the solution the root finder
returned, and every prime non-decomposable target produced the unit ideal
`{1}` (F4's refutation) — 39 of 40 per size.

## 1. Binary: pair-table root finding vs matrix-F4

`PairSumTable::decompose` (one lookup for `m = 2`) versus `groebner_decompose`
(Weil-restrict `S₃ = 0` to a Boolean system over the Frobenius-invariant
subspace, build the Macaulay matrix, row-reduce at degree `D = 2, 3, 4` until
the basis stabilises). Same `K_0` curves, the first 256 scalar multiples
`[k]G` as targets, `m = 2`.

| n | \|F\| | targets | decomposable | root find | matrix-F4 | speed-up |
|---:|---:|---:|---:|---:|---:|---:|
| 9 | 55 | 254 | 254 | 0.73 ms | 279 ms | 383× |
| 13 | 4005 | 256 | 256 | 2.5 ms | 3.9 s | 1535× |

The pair-sum table is built once per curve (`|F|(|F|+1)/2` point additions,
one Montgomery inversion per row); the numbers above are the lookup cost
across all targets, which is what a decomposition sweep pays per target.
Matrix-F4's cost grows with the degree of regularity of the Weil-restricted
system and the Macaulay matrix's column count, not with `|F|`, which is why
the gap *widens* as `n` grows.

**End-to-end cross-check.** Running the whole DLP (`koblitz_index_calculus_dlp`)
with `PairTable` and with `Groebner` recovers the *same* logarithm
(`log = 53` on `K_0/F_{2^9}`). At `|F| = 55` the one-time pair-sum table
build is not amortised inside a single small DLP, so the wall there is
setup-dominated; the per-decomposition table above is the apples-to-apples
oracle cost. Correctness, not speed, is the point of the end-to-end row.

## 2. Prime: direct-subtraction root finding vs `groebner_f4` Buchberger

For a prime short-Weierstrass curve the decomposition ideal of a
2-decomposition is `{ S₃(x_R, X, Y), ∏_i(X − x_i), ∏_i(Y − x_i) }`, with the
factor base encoded as the vanishing polynomial of its x-set (the prime
analogue of the binary subspace constraint). `groebner_f4::buchberger`
computes a Gröbner basis of that ideal — deciding decomposability and
eliminating — while the root finder computes `x(R ∓ F_i)` for each `F_i` in
one batched inversion and looks the result up.

| bits | \|F\| | targets | decomposable | refuted (F4 `{1}`) | root find | `groebner_f4` | speed-up |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 6 | 40 | 1 | 39 | 30 µs | 1.05 s | 35592× |
| 20 | 6 | 40 | 1 | 39 | 33 µs | 1.15 s | 34602× |
| 24 | 6 | 40 | 1 | 39 | 26 µs | 1.23 s | 46911× |

The factor base is kept at `|F| = 6` so Buchberger over the degree-`|F|`
vanishing ideal terminates; the ratio is large because the root finder is
`O(|F|)` field operations while the Gröbner basis of a degree-6 bivariate
ideal is already a second of work. Two fail-closed cross-checks run on every
decomposable target: Semaev's `S₃` must vanish on the returned x-pair, and
every Gröbner basis element (an ideal member) must vanish at that pair.

**Where F4 cannot follow.** The same root finder solves full DLPs at
`|F| = 120–160`, where the degree-`|F|` vanishing ideal is hopeless for
Buchberger:

| bits | \|F\| | targets swept | wall | log verified |
|---:|---:|---:|---:|:--:|
| 20 | 120 | 4 438 | 9.0 ms | yes |
| 24 | 120 | 78 159 | 163 ms | yes |
| 28 | 160 | 1 040 671 | 2.83 s | yes |

## 3. What changed in the pipeline

`KoblitzIcOptions::default()` now selects `DecompositionStrategy::PairTable`
(root finding) instead of `DecompositionStrategy::Groebner`. Gröbner stays
selectable — it is still the engine that emits an *infeasibility
certificate* (the unit ideal) when a target provably does not decompose, and
it is the baseline this note measures against. All 59 Koblitz
index-calculus unit tests pass with the new default, including
`pair_table_agrees_with_enumeration_for_two_three_and_four_summands`, which
pins the root finder's correctness for `m ∈ {2, 3, 4}` against exhaustive
search. The prime engine's 7 unit tests pass, including the ladder-rung
reproduction (a generator-sign fix: `find_a3_curve` now emits the natural
principal square root so it reproduces the committed rungs exactly).

## 4. What this does and does not show

- **No security claim.** Both decompositions still cost `≈ p/|F|` group
  operations per solve for the prime regime and run the same relation /
  rank statistics; this is Gaudry's prime-field verdict unchanged. What
  changed is the *constant*: the per-target decomposition no longer pays a
  Gröbner basis. Nothing here bears on the security of any deployed curve,
  recovers any key, or beats rho asymptotically.
- **The F4 baselines are honest but conservative.** The binary matrix-F4 is
  the crate's own `groebner_decompose` on the identical system. The prime
  `groebner_f4` is the crate's Buchberger on the decomposition ideal; where
  it caps its step budget the measured time is a lower bound, so the stated
  ratios understate F4's cost rather than inflate it.
- **Small-instance caveat.** At tiny `|F|` the one-time table build (binary)
  or the fixed sweep cost (prime) can dominate a *single* DLP, exactly as
  the prime-rho note's set-up caveat describes. The per-decomposition tables
  isolate the oracle cost, which is the quantity the Gröbner solve dominated.

## 5. Scaling to larger fields

Measured by `examples/semaev_rootfind_scaling_bench.rs`
(`research/semaev_rootfind_scaling_20261012/`). Two walls are made
explicit, because "bigger curve" means different things on each axis.

**Prime, native `u64` (`Fp64`), direct-subtraction 2-decomposition,
`|F| = 120`.** Per-target cost stays flat to the `2^62` ceiling; every one
of 200 constructed `R = F_i + F_j` targets decomposes and verifies in the
group at each size.

| bits | per-target root-find | constructed verified |
|---:|---:|---:|
| 40 | 2.95 µs | 200/200 |
| 48 | 3.87 µs | 200/200 |
| 56 | 3.35 µs | 200/200 |
| 61 | 8.82 µs | 200/200 |

**Prime, arbitrary precision (`BigUint`), `p = 2^83 − 55` (largest prime
`< 2^83`), `|F| = 80`.** The *same* algorithm — no summation-polynomial
root solve, no Gröbner basis, just `x(R ∓ F_i)` by subtraction — runs at
83 bits:

| bits | per-target root-find | constructed verified |
|---:|---:|---:|
| 83 | 5.99 ms | 100/100 |

The ≈1000× jump from the `u64` rows is the arbitrary-precision field (one
`BigUint` modular inverse per factor-base element, un-batched), not the
method. A `u128` Montgomery field would recover most of it; the point here
is that the decomposition algorithm, not the word size, is what extends.

**Binary `F_{2^53}`, pairs-and-solve (root-finding) over a bounded
subspace.** Correctness anchor at `ℓ = 6`: `decompose` agrees with an
exhaustive `S₄` triple search on 200/200 targets. Rejection cost (the
dominant per-target work) scales with the subspace as expected:

| ℓ | \|V\| | per-target decompose |
|---:|---:|---:|
| 8 | 256 | 38 ms |
| 10 | 1024 | 585 ms |
| 11 | 2048 | 2.48 s |

**The two walls.**

1. *Native `u64` width (engineering).* `Fp64` needs `p < 2^62`; the binary
   `PairSumTable` / `Gf2` pack an x-coordinate into a `u64` (`n ≤ 63`). An
   `83`-bit prime runs only via `BigUint` (above); `F_{2^83}` exceeds the
   single-word binary field (`Gf2` asserts `n ≤ 63`) and the pair-table
   packing (`n ≤ 62`). A `u128` field (the crate already has `Gf2_128`)
   is the fix, and does not change the algorithm.
2. *Exponential factor base (fundamental).* The Koblitz Frobenius-invariant
   factor base has dimension `ℓ = ord_2(n)`: **52 at `n = 53`, 82 at
   `n = 83`**, so `|V| = 2^52` / `2^82` cannot be materialised by *any*
   decomposition engine — F4 and root finding alike. The small-`n`
   head-to-head of §1 lives where `ord_2(n)` is small; raising `n` makes
   the factor base exponentially larger, not the decomposition harder.
   This is why the scaling study measures the oracle on a *bounded*
   subspace / factor base, isolating the decomposition cost from the
   regime's own barrier.

## 6. Reproduce

```bash
cargo run --release --example semaev_root_vs_groebner_bench
cargo run --release --example semaev_rootfind_scaling_bench
# raw evidence:
BENCH_OUT=research/semaev_root_vs_groebner_20261012 \
  cargo run --release --example semaev_root_vs_groebner_bench
BENCH_OUT=research/semaev_rootfind_scaling_20261012 \
  cargo run --release --example semaev_rootfind_scaling_bench
```
