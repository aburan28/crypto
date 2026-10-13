# Research Note: Root-Finding Decomposition in place of F4/F5 Gröbner for ECDLP Index Calculus

**Date.** 2026-10-12. **Status.** Measured on public synthetic known-answer
instances; every decomposition re-added in the group, and the two oracles
cross-checked to agree on every target. **No security claim** (§4): the
speed-up is on the decomposition *subroutine* at fixed factor-base size; the
asymptotic group-operation count of index calculus is unchanged. No
deployed-curve, key-recovery, or sub-rho claim is made or implied.
**Code.** `src/cryptanalysis/prime_fast.rs` (prime engine),
`src/cryptanalysis/koblitz_index_calculus.rs` (`DecompositionStrategy::PairTable`,
now the default), `examples/semaev_root_vs_groebner_bench.rs` (§1–2 tables),
`examples/semaev_rootfind_scaling_bench.rs` (§5),
`examples/semaev_rootfind_nist_bench.rs` (§6, NIST P-192…P-521);
raw console + JSON under `research/semaev_root_vs_groebner_20261012/`,
`research/semaev_rootfind_scaling_20261012/`,
`research/semaev_rootfind_nist_20261012/`.

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
  §7 measures the end-to-end effect directly and shows F4 winning the three
  smallest binary instances for exactly this reason.

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

## 6. NIST prime curves P-192, P-224, P-256, P-384, P-521

Measured by `examples/semaev_rootfind_nist_bench.rs`
(`research/semaev_rootfind_nist_20261012/`), on the five NIST P-curves with
their real FIPS 186-4 parameters (`a = −3` on all of them, which is exactly
the form the prime engine is built around). The question asked was: does
the root-finding decomposition *keep working* at 192–521-bit field sizes,
and does it still beat F4 there? Both, on every curve.

**Read this first.** Only the per-target decomposition *subroutine* is
measured, on a bounded synthetic factor base of small-`x` points. A fast
decomposition oracle says nothing about the ECDLP on these curves: the
relation-collection and linear-algebra stages of index calculus remain
exponentially out of reach at these sizes and the asymptotics are unchanged
(§4). This section is an engineering result about a subroutine, not a
security result about P-192…P-521.

**A. Head-to-head vs F4, same field, same ideal, same targets.**
`|F| = 6` (the largest at which a Gröbner basis over a 192–521-bit field is
still practical per target), 10 targets per curve: 5 constructed
`R = F_i + F_j`, 5 random on-curve points. Root finding is direct
subtraction with one batched inversion per target; F4 is the crate's
`groebner_f4::buchberger` on `{S₃(X, Y; x_R), ∏(X−x_i), ∏(Y−x_i)}`. Every
target passed the fail-closed cross-checks: a found pair re-adds to `R` in
the group, makes Semaev `S₃` vanish, and makes *every* Gröbner basis
element vanish; every refuted target produced the unit ideal `{1}` from F4.

| curve | bits | root-find / target | F4 / target | factor | decomposed / refuted |
|---|---:|---:|---:|---:|---:|
| P-192 | 192 | 120 µs | 414 ms | **3455×** | 5 / 5 |
| P-224 | 224 | 218 µs | 1.27 s | **5802×** | 5 / 5 |
| P-256 | 256 | 203 µs | 613 ms | **3021×** | 5 / 5 |
| P-384 | 384 | 337 µs | 1.09 s | **3235×** | 5 / 5 |
| P-521 | 521 | 747 µs | 2.54 s | **3395×** | 5 / 5 |

The factor is roughly flat in the field size because both sides pay the
same arbitrary-precision field arithmetic; what separates them is that the
root finder does `O(|F|)` field operations and the Gröbner solve does a
basis computation whose coefficient arithmetic is all in that field.

**B. Root-finder scaling, `|F| = 64`, generic (`BigUint`) field.** 200
random on-curve targets (the rejection path, which dominates in practice)
and 100 constructed `R = F_i ± F_j` controls per curve, every control
re-verified in the group and against `S₃`:

| curve | rejection / target | controls verified | control / target |
|---|---:|---:|---:|
| P-192 | 369 µs | 100/100 | 301 µs |
| P-224 | 428 µs | 100/100 | 268 µs |
| P-256 | 433 µs | 100/100 | 324 µs |
| P-384 | 583 µs | 100/100 | 394 µs |
| P-521 | 1.28 ms | 100/100 | 3.66 ms* |

\* run-to-run jitter on a shared desktop; the previous run measured
1.13 ms. The ratios in table A are stable to within ~2× across runs.

**C. Per-curve specialisation.** Each NIST prime has its own Solinas
reduction, so the natural "specialised code per curve" is a dedicated
field. The crate already ships one for P-256
(`ecc::p256_field::P256FieldElement`, constant-time Montgomery over
`U256`); the identical oracle on it, with the same batched inversion:

| P-256, `|F| = 64` | generic `BigUint` | `P256FieldElement` | factor |
|---|---:|---:|---:|
| rejection / target | 433 µs | 41 µs | **11×** |
| control / target | 324 µs | 28 µs | **12×** |

Both paths agree on every target (0 hits on the 200 random targets,
100/100 controls, same relations). P-192/224/384/521 run on the generic
field in this bench; their Solinas fields are the same follow-on
specialisation and would move them by a similar constant. Nothing about
the algorithm changes per curve — only the field it is instantiated over.

**Closing check.** Root finding beats F4 on all five NIST P-curves, with
every fail-closed check passing; the bench asserts this and exits non-zero
otherwise. The claim is bounded exactly as stated at the top of this
section.

## 7. End to end: complete DLP solves, baseline given every advantage

Measured by `examples/semaev_rootfind_e2e_bench.rs`
(`research/semaev_rootfind_e2e_20261012/`, see its `README.md` for which run
each row comes from). A complete index-calculus solve of a planted log —
factor base, (pair table), relation collection on a seeded target stream,
linear algebra, `[d]G = Q` — with the oracle swapped and everything else
fixed. **CPU time** (process + reaped children), not wall-clock. Every
solve verified. The baseline gets its best configuration: matrix-F4 at its
best `max_degree` (binary); Buchberger at its faster ordering and **msolve**
(an external state-of-the-art F4) charged only its own process CPU via a
paired trivial-ideal control (prime). Composed rows (sampled per-target ×
measured target count) are labelled.

**Binary Koblitz, single solve including the table build:**

| n | \|F\| | log₂ r | root finding | table | matrix-F4 | single solve |
|---:|---:|---:|---:|---:|---:|---|
| 5 | 21 | 4 | 4.3 ms | 1.4 ms | 0.3 ms | **F4 faster** (0.06×) |
| 9 | 55 | 7 | 3.9 ms | 1.5 ms | 1.7 ms | **F4 faster** (0.45×) |
| 13 | 4005 | 11 | 298 ms | 260 ms | 14.2 ms | **F4 faster** (0.05×) |
| 23 | 2025 | 21 | 90.8 ms | 89.2 ms | 11.1 s | **122×** |

Why F4 wins three of four: at `n ≤ 13` the subgroup is tiny, a solve needs a
few hundred decompositions, and matrix-F4 on those small descended systems
costs tens of µs per target — less than building `C(|F|+1, 2)` pair sums
that are then barely used. At `n = 23` F4 costs seconds and the 89 ms table
pays for itself in the first solve.

**Prime a = −3, all three oracles on the identical target stream:**

| bits | \|F\| | targets | root finding | Buchberger | msolve (F4) | vs msolve |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 8 | 6,544 | 23.1 ms | 133 s (5767×) | 2.31 s | **100×** |
| 20 | 12 | 75,683 | 416.1 ms | 6,614 s COMPOSED (15895×) | 43.77 s COMPOSED | **105×** |
| 24 | 12 | 1,501,376 | 12.30 s | 134,492 s COMPOSED (10937×) | 619 s COMPOSED | **50×** |

**Each method at its own best |F| (16-bit):** msolve's composed total is minimised at |F| = 16 (1.41 s); the root finder measured in full at |F| = 120 takes 6.9 ms (829 targets) → **205×** at each method's own best.

Why end-to-end is two orders of magnitude below per-target: the root
finder's solve is dominated by *generating* targets (a scalar multiplication
each), not deciding them — at 16 bits the oracle is 16% of its
total. Replacing a ms-scale Gröbner oracle with a µs one removes the oracle
from the critical path; what remains sets the ratio. The per-target tables
say how much faster the oracle got; this one says how much faster the solve
got; the second is the smaller number, as it must be.

Not measured end to end: `m = 3, 4`, the NIST curves, Buchberger beyond
16 bits (composed). At those sizes a complete solve is out of reach for every
method.

## 8. Reproduce

```bash
cargo run --release --example semaev_root_vs_groebner_bench
cargo run --release --example semaev_rootfind_scaling_bench
cargo run --release --example semaev_rootfind_nist_bench
# raw evidence:
BENCH_OUT=research/semaev_root_vs_groebner_20261012 \
  cargo run --release --example semaev_root_vs_groebner_bench
BENCH_OUT=research/semaev_rootfind_scaling_20261012 \
  cargo run --release --example semaev_rootfind_scaling_bench
F4_CONSTRUCTED=5 F4_RANDOM=5 BENCH_OUT=research/semaev_rootfind_nist_20261012 \
  cargo run --release --example semaev_rootfind_nist_bench
# end to end (needs msolve on PATH for the external-F4 column):
BENCH_OUT=research/semaev_rootfind_e2e_20261012 \
  cargo run --release --example semaev_rootfind_e2e_bench
E2E_SKIP_BINARY=1 E2E_FULL_BUDGET_S=300 BENCH_OUT=research/semaev_rootfind_e2e_20261012/run2 \
  cargo run --release --example semaev_rootfind_e2e_bench
```
