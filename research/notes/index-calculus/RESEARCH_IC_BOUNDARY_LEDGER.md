# The index-calculus boundary ledger: three regimes, one unit, every phase priced

**Harness:** `src/cryptanalysis/ic_boundary.rs` (pipelines, rho reference, ledger),
`src/cryptanalysis/ic_oracle_pricing.rs` (decomposition oracles per target),
`src/cryptanalysis/ic_corpus.rs` (Trimoska-style instance corpus)
**Command:** `cargo build --release --bin ic && ./target/release/ic boundary --oracles --out <file>`
**Frozen run:** `docs/ic/runs/ic-boundary-ledger-2026-09-21.json` (every number below is
read from it: the tables are rendered from its `ledger` and `oracle_pricing` fields, and the
report's own `markdown` field carries the same ladder; status `complete`, `all_verified`
true, 1,855 s on one Intel Xeon 2.10 GHz host, four cores, at commit `a696a93`)
**Corpus:** `docs/ic/corpus/` (`ic corpus`)
**Reporting rule:** `AGENTS.md` — boundary first, one table in one unit, ratio columns,
class labels, phases priced separately, falsification target declared in advance.

## Why this exists

The repository had every ingredient of an index-calculus attack on
elliptic curves — Semaev's summation polynomials `S₃…S₆` and their
symmetrised forms, Weil descent, a Boolean matrix-F4 with splitting, a
CDCL solver with native parity rows, XL, block Wiedemann, Frobenius-
invariant factor bases, pair tables, the Nagao encodings — and it had
measurements of most of them, in different units, on different curves,
at different sizes, with different things left out of the bill.  What
it did not have was the one thing `AGENTS.md` says a thread must state
before it starts: **a boundary, in a unit in which the boundary is a
constant, with every variant on the same axis.**  The prime-field
bench reported wall-clock exponents, the Koblitz ledger reported
charged wall ratios at four degrees, the binary work reported oracle
timings per target, and none of the three could be read against the
others.

This note fixes the axis.  Three regimes — a generic prime-field
curve, a generic (random) binary curve, and a Koblitz curve — run the
same planted discrete logarithm end to end, with the same relation
loop and the same elimination, and every phase of every variant is
counted in one unit against two boundaries that were written down
before anything was run.  The result is a table a later round can beat
or fail to beat by reading one column, and per-phase exponents that
say which phase would have to move.

## 1. The boundaries, stated before measuring

### 1.1 The unit

```
    S = total operations / √r
```

`r` is the prime order of the subgroup the logarithm lives in; an
*operation* is one affine group addition on the curve in question
(a doubling counts as one).  Every phase is inside the number: factor
base, target generation, decomposition oracle, pair table, linear
algebra, verification.  Work that is not a group addition — a modular
square root, an Artin–Schreier solve, one pair of the pairs-and-solve
loop, a hash probe, one multiply-subtract of the elimination, a
64-bit word XOR of a Macaulay reduction — is counted **exactly** in its
own native unit and converted at a factor **measured on the host at run
time** (nanoseconds per native unit over nanoseconds per addition on the
same curve representation).  Both the counts and the factors are in the
frozen report, so any reader can re-convert.  The counts are the
measurement; the factor is the only thing hardware touches.  Wall
time is carried as `s_wall` for practicality and is never the metric.

### 1.2 Boundary 1: the generic floor

A generic algorithm on a group of prime order `r` whose curve offers an
automorphism group of order `A` needs about `√(πr/2A)` operations
(Shoup's `Ω(√r)` with the `√A` the automorphisms buy).  In the unit:

| regime | automorphisms `A` | `S_floor = √(π/2A)` |
|:--|--:|--:|
| prime, generic curve | 2 (negation) | 0.886 |
| binary, random curve | 2 (negation) | 0.886 |
| Koblitz over `F_{2^n}` | `2n` (signed Frobenius) | `√(π/4n)`: 0.229 at `n = 15`, 0.138 at `n = 41` |

Nothing generic crosses it, and nothing in this note claims to.  Its
job is to be the denominator of the `vs floor` column so that the size
of the gap is a number rather than a feeling.

### 1.3 Boundary 1′: the counting ceiling on relations

A base of `F` signed points has at most `C(F+m−1, m)` `m`-sums, so a
uniform target of the group of order `#E` decomposes with probability at
most `p_ceiling = min(1, C(F+m−1, m)/#E)`, and the `K + 1` relations the
elimination needs cost at least `(K+1)/p_ceiling` targets.  The report
carries the measured relations-per-trial over `p_ceiling`
(`yield_over_ceiling`): a value near one says the base is as good as a
random set of its size and nothing structural is being exploited in the
yield; a value well above one would be a finding.

### 1.4 Boundary 2: the reference, measured

Pollard rho, on the same instance, in the same process, counted
exactly and verified as `[d]G = Q` on every run:

- prime and random-binary regimes: an r-adding walk with 16 jumps and
  distinguished points (van Oorschot–Wiener), no automorphism used, so
  its expected walk is `√(πr/2)` and its `A = 1`; the floor above is
  stated at `A = 2` because the negation map is a generic technique the
  reference does not implement, and the report carries the walk-only
  `s_walk` next to the total `s`;
- Koblitz regime: the repository's signed-Frobenius walk
  (`koblitz_signed_frobenius_rho_reference`), which quotients by all
  `2n` automorphisms, with its exact operation ledger.

Eight walks per target, three targets per instance.  A rho that does
not recover and verify the planted logarithm makes every ratio in its
row meaningless, so the report records `rho_verified_all` per instance
and the ledger's status is `complete` only when it is true everywhere.

### 1.5 What the three regimes are

| regime | curve | subgroup sizes | factor base | column map | oracles |
|:--|:--|:--|:--|:--|:--|
| `prime` | `y² = x³ + ax + b` over `F_p`, prime order (roster curves to 20 bits, generated prime-order curves at 22 and 24) | 2¹⁰ … 2²⁴ | the `2^{⌈bits/3⌉}` smallest abscissae, both signs | one per abscissa | Semaev `S₃` roots (`m = 2`), direct subtraction (`m = 2`), meet in the middle (`m = 2, 3`) |
| `char2` | random `y² + xy = x³ + ax² + b` over `F_{2^n}`, cofactor ≤ 8 | `n = 15, 18, 21, 24, 27` | low-order subspace `⟨1, z, …, z^{l−1}⟩`, `l = ⌈n/3⌉` | one per abscissa | `S₄` pairs-and-solve (`m = 3`), meet in the middle |
| `koblitz` | `K_a` over `F_{2^n}` | `n = 11 … 41` (the degrees with a usable subgroup) | Frobenius-invariant subspace of dimension near `⌈n/3⌉`, or an orbit union where none exists | one per signed Frobenius orbit (`λ^k` coefficients), and the same base with one column per abscissa as the control | meet in the middle, `S₄` pairs-and-solve on the invariant subspace |

All three share one relation loop (`R = [a]G + [b]Q`, oracle, sign lift,
one row) and one incremental Gauss–Jordan over `Z/rZ` that stops the
moment the target's column is pinned.  The only thing that differs
between a `char2` row and a `koblitz` row on the same curve is the
column map, which is exactly what the Frobenius is supposed to buy.

### 1.6 The falsification target, declared in advance

This round establishes baselines; the class of every later change is
judged against them by the test in `AGENTS.md` §3.  The conditions a
later round would have to meet to be called a result:

- **advance:** a variant whose fitted total exponent (`ops ∝ r^α`, four
  or more sizes, `R² ≥ 0.95`) is below the reference's, with every
  logarithm verified and every phase inside the count; or a
  `yield_over_ceiling` above `1.5` on a base that is not inside a
  proper subgroup, reproduced on a second seed;
- **crossover:** `S / S_rho < 1` on any instance with `r ≥ 2^{20}`, all
  phases charged, the reference verified on that instance;
- **abandon:** two consecutive rounds whose only movement is in the
  conversion factors or in the wall column.

Inadmissible: changing the unit, dropping a phase from the count,
choosing favourable seeds, counting dependent relations, or reading a
per-target oracle price as a pipeline cost without the trials it
implies.

A later round that claims a *speedup* over this baseline is also bound
by `AGENTS.md` §8: the frozen solver regression suite under
`research/index_calculus_baseline_20260914/regression/` run on
baseline and candidate, a matched full-DLP comparison on identical
instances with `speedup = baseline_total_operations /
candidate_total_operations` in this unit, and the evidence committed
with the claim.  The ledger here is the end-to-end side of that gate:
the counts a candidate has to lower, on the instances it has to lower
them on.

## 2. One table, one unit

Unit `S` throughout.  `vs rho` divides by the counted reference on the
same instance (setup and walk, eight walks per target); `vs floor`
divides by `√(π/2A)`.  `m` is the number of summands, `|F|` the signed
base size, `K` the number of unknowns (columns), `trials` the targets
drawn, and `yield/ceiling` the measured relations per trial over
`p_ceiling` of §1.3.  Means over three targets.  Every row's logarithm
was recovered and verified as `[d]G = Q`, and every rho run recovered
and verified its own; the class column is by the test of `AGENTS.md`
§3, argued row by row in §3.

### 2.1 The table: every variant on the largest instance of its regime

| regime, instance | variant | m | \|F\| | K | trials | yield/ceiling | S | vs rho | vs floor | correct | class |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| **prime**, `generated-24bit-10935329`, `r = 2^23.4`, `#E = r`, `A = 2` | generic floor `√(π/2A)` | | | | | | 0.886 | 0.23× | 1× | — | boundary |
| | Pollard rho, r-adding, counted (walk alone 3.54) | | | | | | 3.93 | 1× | 4.4× | ✓ | reference |
| | Semaev `S₃` roots | 2 | 512 | 256 | 9,784 | 1.05 | 4,349 | 1,107× | 4,907× | ✓ | baseline |
| | direct subtraction | 2 | 512 | 256 | 9,784 | 1.05 | 1,715 | 436× | 1,935× | ✓ | accounting |
| | meet in the middle | 2 | 512 | 256 | 9,784 | 1.05 | 232 | 59.1× | 262× | ✓ | engineering |
| | meet in the middle | 3 | 512 | 256 | 265 | 0.87 | **56.8** | **14.5×** | **64.1×** | ✓ | engineering |
| **binary**, `random-binary-n27-b845462`, `r = 2^24.4`, `#E = 6r`, `A = 2` | generic floor | | | | | | 0.886 | 0.38× | 1× | — | boundary |
| | Pollard rho, r-adding, counted (walk alone 1.98) | | | | | | 2.31 | 1× | 2.6× | ✓ | reference |
| | meet in the middle | 3 | 526 | 263 | 1,498 | 0.90 | **206** | **89.0×** | **232×** | ✓ | baseline |
| | `S₄` pairs-and-solve | 3 | 526 | 263 | 1,499 | 0.90 | 105,198 | 45,471× | 118,703× | ✓ | relabelling |
| **Koblitz**, `K_1 / GF(2^23)`, `r = 2^22.0`, `#E = 2r`, `A = 46` | generic floor | | | | | | 0.185 | 0.21× | 1× | — | boundary |
| | signed-Frobenius rho, counted (walk alone 0.19) | | | | | | 0.89 | 1× | 4.8× | ✓ | reference |
| | meet in the middle, one column per abscissa | 3 | 875 | 438 | 348 | 1.00 | 203 | 227× | 1,096× | ✓ | baseline (control) |
| | meet in the middle, signed-orbit columns | 3 | 875 | 20 | 18 | 1.00 | 188 | 210× | 1,017× | ✓ | advance, count |
| **Koblitz**, `K_0 / GF(2^41)`, `r = 2^39.0`, `#E = 4r`, `A = 82` | generic floor | | | | | | 0.138 | 0.69× | 1× | — | boundary |
| | signed-Frobenius rho, counted (walk alone 0.19) | | | | | | 0.20 | 1× | 1.45× | ✓ | reference |
| | meet in the middle, signed-orbit columns | 3 | 5,003 | 62 | 6,086 | 1.05 | **58.8** | **300×** | **425×** | ✓ | advance, count |
| **Koblitz**, `K_0 / GF(2^31)`, `r = 2^20.5`, `#E = 1492r`, `A = 62` | signed-Frobenius rho, counted | | | | | | 1.13 | 1× | 7.1× | ✓ | reference |
| | meet in the middle, one column per abscissa | 3 | 2,421 | 1,211 | 951 | 0.67 | 3,412 | 3,010× | 21,436× | ✓ | baseline (control) |
| | meet in the middle, signed-orbit columns | 3 | 2,421 | 41 | 49 | 0.74 | 2,489 | 2,195× | 15,635× | ✓ | advance, count |
| | `S₄` pairs-and-solve on the invariant subspace | 3 | 2,421 | 41 | 49 | 0.74 | 128,758 | 113,579× | 808,926× | ✓ | relabelling |

The abscissa-column control runs to `n = 31` (`koblitz_no_fold_max_degree`),
so the `n = 41` rung has no control row; `n = 23` is the largest
cofactor-2 rung on which both are measured, and `n = 31` the largest on
which the `S₄` oracle is.

### 2.2 The full ladder

Every instance, every variant, every phase (`FB` factor base and pair
table, `rel` targets and oracle, `LA` elimination, `verify`), in `S`.
`S_wall` is the same total priced by wall time instead of counts, a
practicality note only.

| regime | instance | log₂ r | variant | m | \|F\| | K | trials | yield/ceiling | S | S_wall | vs rho | vs floor | FB | rel | LA | verify | ok |
|:--|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| prime | bench-10bit | 9.7 | semaev_s3_roots_m2 | 2 | 32 | 16 | 31 | 0.75 | 43.3 | 58.0 | 1.91 | 48.9 | 0.64 | 42.2 | 0.14 | 0.349 | ✓ |
| prime | bench-10bit | 9.7 | direct_subtraction_m2 | 2 | 32 | 16 | 31 | 0.75 | 51.3 | 54.4 | 2.27 | 57.9 | 0.64 | 50.2 | 0.14 | 0.349 | ✓ |
| prime | bench-10bit | 9.7 | mitm_m2 | 2 | 32 | 16 | 31 | 0.75 | 45.9 | 45.4 | 2.03 | 51.8 | 19.0 | 26.4 | 0.14 | 0.349 | ✓ |
| prime | bench-10bit | 9.7 | mitm_m3 | 3 | 32 | 16 | 14 | 1.00 | 33.1 | 33.9 | 1.46 | 37.4 | 19.0 | 13.5 | 0.25 | 0.349 | ✓ |
| prime | bench-12bit | 11.9 | semaev_s3_roots_m2 | 2 | 32 | 16 | 95 | 1.03 | 75.0 | 88.1 | 4.04 | 84.7 | 0.34 | 74.4 | 0.05 | 0.214 | ✓ |
| prime | bench-12bit | 11.9 | direct_subtraction_m2 | 2 | 32 | 16 | 95 | 1.03 | 92.3 | 84.4 | 4.97 | 104 | 0.34 | 91.7 | 0.05 | 0.214 | ✓ |
| prime | bench-12bit | 11.9 | mitm_m2 | 2 | 32 | 16 | 95 | 1.03 | 55.7 | 55.1 | 3.00 | 62.8 | 8.81 | 46.6 | 0.05 | 0.214 | ✓ |
| prime | bench-12bit | 11.9 | mitm_m3 | 3 | 32 | 16 | 18 | 0.83 | 21.7 | 22.0 | 1.17 | 24.5 | 8.81 | 12.5 | 0.11 | 0.214 | ✓ |
| prime | bench-14bit | 14.0 | semaev_s3_roots_m2 | 2 | 64 | 32 | 161 | 1.25 | 102 | 116 | 11.5 | 116 | 0.50 | 102 | 0.05 | 0.141 | ✓ |
| prime | bench-14bit | 14.0 | direct_subtraction_m2 | 2 | 64 | 32 | 161 | 1.25 | 122 | 109 | 13.6 | 137 | 0.50 | 121 | 0.05 | 0.141 | ✓ |
| prime | bench-14bit | 14.0 | mitm_m2 | 2 | 64 | 32 | 161 | 1.25 | 64.2 | 58.4 | 7.20 | 72.4 | 16.8 | 47.1 | 0.05 | 0.141 | ✓ |
| prime | bench-14bit | 14.0 | mitm_m3 | 3 | 64 | 32 | 28 | 0.96 | 27.9 | 22.5 | 3.13 | 31.5 | 16.8 | 10.8 | 0.13 | 0.141 | ✓ |
| prime | bench-16bit | 16.0 | semaev_s3_roots_m2 | 2 | 128 | 64 | 360 | 0.92 | 594 | 654 | 93.3 | 670 | 1.66 | 592 | 0.04 | 0.081 | ✓ |
| prime | bench-16bit | 16.0 | direct_subtraction_m2 | 2 | 128 | 64 | 360 | 0.92 | 231 | 213 | 36.4 | 261 | 1.66 | 230 | 0.04 | 0.081 | ✓ |
| prime | bench-16bit | 16.0 | mitm_m2 | 2 | 128 | 64 | 360 | 0.92 | 94.9 | 86.8 | 14.9 | 107 | 34.0 | 60.7 | 0.04 | 0.081 | ✓ |
| prime | bench-16bit | 16.0 | mitm_m3 | 3 | 128 | 64 | 53 | 0.99 | 45.7 | 34.3 | 7.18 | 51.5 | 34.0 | 11.3 | 0.23 | 0.081 | ✓ |
| prime | bench-18bit | 18.0 | semaev_s3_roots_m2 | 2 | 128 | 64 | 674 | 1.13 | 192 | 220 | 45.6 | 216 | 0.24 | 191 | 0.01 | 0.041 | ✓ |
| prime | bench-18bit | 18.0 | direct_subtraction_m2 | 2 | 128 | 64 | 674 | 1.13 | 232 | 204 | 55.2 | 262 | 0.24 | 232 | 0.01 | 0.041 | ✓ |
| prime | bench-18bit | 18.0 | mitm_m2 | 2 | 128 | 64 | 674 | 1.13 | 81.1 | 72.7 | 19.3 | 91.6 | 16.4 | 64.7 | 0.01 | 0.041 | ✓ |
| prime | bench-18bit | 18.0 | mitm_m3 | 3 | 128 | 64 | 78 | 0.75 | 32.2 | 23.7 | 7.65 | 36.3 | 16.4 | 15.6 | 0.17 | 0.041 | ✓ |
| prime | bench-20bit | 20.0 | semaev_s3_roots_m2 | 2 | 256 | 128 | 1746 | 1.20 | 439 | 518 | 164 | 496 | 0.24 | 439 | 0.01 | 0.023 | ✓ |
| prime | bench-20bit | 20.0 | direct_subtraction_m2 | 2 | 256 | 128 | 1746 | 1.20 | 525 | 474 | 196 | 593 | 0.24 | 525 | 0.01 | 0.023 | ✓ |
| prime | bench-20bit | 20.0 | mitm_m2 | 2 | 256 | 128 | 1746 | 1.20 | 126 | 113 | 47.1 | 143 | 32.4 | 94.0 | 0.01 | 0.023 | ✓ |
| prime | bench-20bit | 20.0 | mitm_m3 | 3 | 256 | 128 | 118 | 0.93 | 45.0 | 30.3 | 16.8 | 50.8 | 32.4 | 12.3 | 0.25 | 0.023 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | semaev_s3_roots_m2 | 2 | 512 | 256 | 2793 | 1.03 | 692 | 838 | 232 | 781 | 0.29 | 692 | 0.01 | 0.016 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | direct_subtraction_m2 | 2 | 512 | 256 | 2793 | 1.03 | 870 | 780 | 292 | 981 | 0.29 | 869 | 0.01 | 0.016 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | mitm_m2 | 2 | 512 | 256 | 2793 | 1.03 | 165 | 136 | 55.3 | 186 | 72.7 | 92.1 | 0.01 | 0.016 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | mitm_m3 | 3 | 512 | 256 | 218 | 1.00 | 84.3 | 53.1 | 28.3 | 95.1 | 72.7 | 11.1 | 0.51 | 0.016 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | semaev_s3_roots_m2 | 2 | 512 | 256 | 9784 | 1.05 | 4,349 | 4,654 | 1,107 | 4,907 | 0.46 | 4,348 | 0.01 | 0.009 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | direct_subtraction_m2 | 2 | 512 | 256 | 9784 | 1.05 | 1,715 | 1,561 | 436 | 1,935 | 0.46 | 1,715 | 0.01 | 0.009 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | mitm_m2 | 2 | 512 | 256 | 9784 | 1.05 | 232 | 219 | 59.1 | 262 | 40.2 | 192 | 0.01 | 0.009 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | mitm_m3 | 3 | 512 | 256 | 265 | 0.87 | 56.8 | 38.1 | 14.5 | 64.1 | 40.2 | 16.1 | 0.51 | 0.009 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | mitm_m3 | 3 | 30 | 15 | 106 | 0.95 | 69.1 | 56.3 | 8.27 | 78.0 | 15.5 | 53.4 | 0.04 | 0.144 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | semaev_s4_pairs_and_solve_m3 | 3 | 30 | 15 | 117 | 0.89 | 1,294 | 1,283 | 155 | 1,460 | 11.9 | 1,282 | 0.05 | 0.144 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | mitm_m3 | 3 | 64 | 32 | 97 | 1.61 † | 80.0 | 60.6 | 10.6 | 90.2 | 30.2 | 49.6 | 0.05 | 0.101 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | semaev_s4_pairs_and_solve_m3 | 3 | 64 | 32 | 97 | 1.48 † | 2,652 | 2,643 | 351 | 2,992 | 18.7 | 2,633 | 0.04 | 0.101 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | mitm_m3 | 3 | 122 | 61 | 407 | 0.96 | 81.0 | 70.8 | 25.4 | 91.4 | 15.2 | 65.7 | 0.04 | 0.026 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | semaev_s4_pairs_and_solve_m3 | 3 | 122 | 61 | 408 | 0.96 | 8,439 | 8,338 | 2,652 | 9,522 | 7.89 | 8,431 | 0.04 | 0.026 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | mitm_m3 | 3 | 274 | 137 | 365 | 1.74 † | 102 | 89.6 | 40.5 | 115 | 36.7 | 65.0 | 0.15 | 0.020 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | semaev_s4_pairs_and_solve_m3 | 3 | 274 | 137 | 365 | 1.74 † | 17,425 | 17,346 | 6,926 | 19,661 | 10.6 | 17,414 | 0.15 | 0.020 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | mitm_m3 | 3 | 526 | 263 | 1498 | 0.90 | 206 | 194 | 89.0 | 232 | 37.3 | 168 | 0.20 | 0.007 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | semaev_s4_pairs_and_solve_m3 | 3 | 526 | 263 | 1499 | 0.90 | 105,198 | 104,992 | 45,471 | 118,703 | 8.04 | 105,190 | 0.20 | 0.007 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | mitm_m3_signed_orbit_columns | 3 | 45 | 3 | 5 | 0.82 | 39.0 | 38.6 | 2.24 | 146 | 33.0 | 5.58 | 0.02 | 0.328 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | mitm_m3_abscissa_columns_control | 3 | 45 | 23 | 20 | 0.90 | 54.4 | 57.1 | 3.13 | 203 | 33.0 | 20.7 | 0.27 | 0.328 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | mitm_m3_signed_orbit_columns | 3 | 79 | 4 | 4 | 1.00 | 74.0 | 70.4 | 5.51 | 301 | 70.8 | 2.83 | 0.01 | 0.290 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | mitm_m3_abscissa_columns_control | 3 | 79 | 40 | 29 | 1.00 | 92.1 | 89.8 | 6.85 | 375 | 70.8 | 20.7 | 0.28 | 0.290 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m3_signed_orbit_columns | 3 | 33 | 3 | 8 | 2.70 † | 36.5 | 37.0 | 1.91 | 160 | 20.6 | 15.6 | 0.01 | 0.316 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m3_abscissa_columns_control | 3 | 33 | 17 | 46 | 1.42 † | 106 | 103 | 5.56 | 465 | 20.6 | 85.4 | 0.16 | 0.316 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | semaev_s4_pairs_and_solve_m3_signed_orbit_columns | 3 | 33 | 3 | 8 | 2.70 † | 369 | 380 | 19.2 | 1,610 | 0.14 | 368 | 0.01 | 0.316 | ✓ |
| koblitz | K_1 / GF(2^17) | 16.0 | mitm_m2_signed_orbit_columns | 2 | 239 | 8 | 20 | 1.55 † | 116 | 111 | 32.2 | 538 | 112 | 3.44 | 0.00 | 0.082 | ✓ |
| koblitz | K_1 / GF(2^17) | 16.0 | mitm_m2_abscissa_columns_control | 2 | 239 | 120 | 231 | 1.56 † | 151 | 148 | 42.1 | 704 | 112 | 39.1 | 0.05 | 0.082 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | mitm_m3_signed_orbit_columns | 3 | 305 | 9 | 9 | 1.00 | 92.2 | 94.6 | 44.2 | 454 | 91.1 | 1.06 | 0.00 | 0.048 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | mitm_m3_abscissa_columns_control | 3 | 305 | 153 | 126 | 1.00 | 108 | 111 | 51.5 | 529 | 91.1 | 16.0 | 0.37 | 0.048 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | mitm_m3_signed_orbit_columns | 3 | 875 | 20 | 18 | 1.00 | 188 | 187 | 210 | 1,017 | 187 | 0.79 | 0.00 | 0.015 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | mitm_m3_abscissa_columns_control | 3 | 875 | 438 | 348 | 1.00 | 203 | 204 | 227 | 1,096 | 187 | 14.6 | 0.79 | 0.015 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | mitm_m3_signed_orbit_columns | 3 | 3771 | 66 | 55 | 1.00 | 34,554 | 38,510 | 8,226 | 209,968 | 34,518 | 35.4 | 0.18 | 0.107 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | mitm_m3_abscissa_columns_control | 3 | 3771 | 1886 | 216 | 1.00 | 34,650 | 38,640 | 8,249 | 210,548 | 34,518 | 131 | 0.15 | 0.107 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | mitm_m3_signed_orbit_columns | 3 | 2421 | 41 | 49 | 0.74 | 2,489 | 2,696 | 2,195 | 15,635 | 2,444 | 44.7 | 0.01 | 0.023 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | mitm_m3_abscissa_columns_control | 3 | 2421 | 1211 | 951 | 0.67 | 3,412 | 3,649 | 3,010 | 21,436 | 2,444 | 949 | 18.7 | 0.023 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | semaev_s4_pairs_and_solve_m3_signed_orbit_columns | 3 | 2421 | 41 | 49 | 0.74 | 128,758 | 128,755 | 113,579 | 808,926 | 0.20 | 128,757 | 0.01 | 0.023 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m3_signed_orbit_columns | 3 | 4663 | 64 | 493 | 0.98 | 857 | 1,006 | 2,052 | 5,882 | 716 | 141 | 0.00 | 0.003 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m3_signed_orbit_columns | 3 | 4681 | 61 | 921 | 2.05 † | 1,830 | 2,125 | 4,021 | 12,898 | 1,323 | 507 | 0.00 | 0.004 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m3_signed_orbit_columns | 3 | 5003 | 62 | 6086 | 1.05 | 58.8 | 65.9 | 300 | 425 | 16.9 | 41.9 | 0.00 | 0.000 | ✓ |

† the base lies inside a proper subgroup of `E`, where the ceiling of
§1.3 is loose by the subgroup's index; see §3.5, which is a correction
to the boundary and not a yield finding.

### 2.3 The instances and the reference

`#E/r` is the cofactor, `A` the automorphism count the floor and the
reference use, `rho S` the counted total (setup and walk), `rho walk S`
the walk alone, `steps/expected` the measured walk length over
`√(πr/2A)` for the reference's own `A` (`A = 1` for the r-adding walk,
`2n` for the signed-Frobenius walk).

| regime | instance | log₂ r | #E/r | A | S_floor | rho S | rho walk S | steps/expected | rho ok |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|
| prime | bench-10bit | 9.7 | 1 | 2 | 0.886 | 22.7 | 1.81 | 1.44 | ✓ |
| prime | bench-12bit | 11.9 | 1 | 2 | 0.886 | 18.6 | 4.23 | 3.38 | ✓ |
| prime | bench-14bit | 14.0 | 1 | 2 | 0.886 | 8.91 | 1.47 | 1.17 | ✓ |
| prime | bench-16bit | 16.0 | 1 | 2 | 0.886 | 6.36 | 2.35 | 1.87 | ✓ |
| prime | bench-18bit | 18.0 | 1 | 2 | 0.886 | 4.20 | 1.95 | 1.55 | ✓ |
| prime | bench-20bit | 20.0 | 1 | 2 | 0.886 | 2.68 | 1.59 | 1.27 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | 1 | 2 | 0.886 | 2.98 | 2.28 | 1.82 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | 1 | 2 | 0.886 | 3.93 | 3.54 | 2.83 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | 2 | 2 | 0.886 | 8.36 | 1.39 | 1.11 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | 8 | 2 | 0.886 | 7.56 | 1.38 | 1.10 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | 2 | 2 | 0.886 | 3.18 | 2.03 | 1.62 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | 8 | 2 | 0.886 | 2.52 | 1.80 | 1.44 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | 6 | 2 | 0.886 | 2.31 | 1.98 | 1.58 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | 2 | 22 | 0.267 | 17.4 | 0.28 | 1.05 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | 4 | 26 | 0.246 | 13.4 | 0.21 | 0.86 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | 44 | 30 | 0.229 | 19.2 | 0.27 | 1.19 | ✓ |
| koblitz | K_1 / GF(2^17) | 16.0 | 2 | 34 | 0.215 | 3.59 | 0.25 | 1.16 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | 2 | 38 | 0.203 | 2.09 | 0.21 | 1.04 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | 2 | 46 | 0.185 | 0.89 | 0.19 | 1.03 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | 12646 | 58 | 0.165 | 4.20 | 0.21 | 1.25 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | 1492 | 62 | 0.159 | 1.13 | 0.17 | 1.09 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | 596 | 74 | 0.146 | 0.42 | 0.15 | 1.02 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | 8012 | 78 | 0.142 | 0.46 | 0.13 | 0.89 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | 4 | 82 | 0.138 | 0.20 | 0.19 | 1.36 | ✓ |

The reference is not flat at these sizes and the table says why: the
walk alone sits at `1.4–4.2` in the prime and binary regimes against the
`1.25` its `A = 1` expectation gives, and at `0.13–0.28` on the Koblitz
curves against `√(π/4n)`, but the setup (sixteen jump points and the
`[a]G + [b]Q` start of each of eight walks, about 70 additions per
scalar at 24 bits) is a constant that `√r` does not amortise below
`2^20`.  So `rho S` falls from `22.7` to `3.9` over the prime ladder and
from `17.4` to `0.20` over the Koblitz one, and the `vs rho` column of
every variant *rises* with size for that reason alone.  The `vs floor`
column has no such drift.

## 3. Reading the table

### 3.1 The bottom line, against the targets of §1.6

No row meets either condition for a result.  No variant is below the
reference on any instance with `r ≥ 2^20` (the closest is the prime
meet-in-the-middle at `m = 3`, `14.5×` at `2^23.4`), and no fitted
total exponent is below the reference's (§4: every variant is at
`0.54` or above against a reference that fits at `0.25–0.30` here and
is `0.5` in the limit).  The `yield/ceiling` column never exceeds one
on a base that is not confined to a proper subgroup (§3.5).  What the
round establishes is the baseline table itself, with its class labels,
and the phase exponents that say what a later round would have to
move.

### 3.2 Prime field: the base costs `r^{2/3}`, the oracle is not the problem

Four variants, one relation loop, one elimination, identical trials
and identical relations for the three `m = 2` oracles (they see the
same targets and decide the same way, so their `trials` and
`yield/ceiling` columns agree to the digit).  Only the price of a
decision differs:

- **Semaev `S₃` roots (baseline).**  One modular square root per base
  abscissa per target.  Its `S` is set by the host's square-root
  factor, which is not a property of the algorithm: `63–171 ns` on six
  of the primes and `993–1,060 ns` on the 16- and 24-bit ones (Appendix
  A), which is why the row jumps to `4,349` at 24 bits.
- **Direct subtraction: accounting.**  The same enumeration paid in
  additions instead of roots.  Same trials, same yield, same
  elimination; the movement is entirely in the conversion factor, and
  it goes the other way at 10 bits (`51.3` against `43.3`).  Not a gain
  and not claimed as one.
- **Meet in the middle, `m = 2`: engineering.**  The `|F|²/2` pair sums
  are tabled once (`FB` rises from `0.46` to `40.2` at 24 bits) and each
  target costs one probe.  Trials and yield unchanged, so the ratio to
  the counting boundary is unmoved; `S` falls `7.4×` against direct
  subtraction at 24 bits.  With the oracle reduced to a probe, the
  relation phase is dominated by *making the target*: `192 × √r` over
  `9,784` trials is about `65` additions per trial, which is the two
  scalar multiplications of `R = [a]G + [b]Q`.
- **Meet in the middle, `m = 3`: engineering, and the best prime row.**
  Three summands raise the ceiling from `C(F+1,2)/r ≈ 0.012` to
  `C(F+2,3)/r ≈ 1`, so `265` trials replace `9,784`, at `|F|`
  subtractions and probes each.  `S = 56.8`, `14.5×` the reference,
  `64×` the floor, and `71%` of it is the pair table.  Over the ladder
  its `S` stays between `21.7` and `84.3`, oscillating with the base
  size `2^{⌈bits/3⌉}` that doubles every three bits, while the reference
  falls from `22.7` to `3.9`; that is the whole reason its `vs rho`
  climbs from `1.2×` at 12 bits to `14.5×` at 24.

The `yield/ceiling` column sits at `0.75–1.25` for `m = 2` and
`0.75–1.00` for `m = 3`: the smallest-abscissa base behaves like a
random set of its size, and nothing in the prime pipeline exploits
structure in the yield.

### 3.3 Random binary curve: same loop, and the algebraic oracle is a relabelling

- **Meet in the middle, `m = 3` (baseline).**  `S` from `69` to `206`
  over `n = 15…27`, `8×` to `89×` the reference, `78×` to `232×` the
  floor.  At `n = 27` the relation phase is `82%` of the total (`168` of
  `206`), at about `516` additions per trial, which is the `|F| = 526`
  subtractions the `m = 3` probe needs.  The pair table is `18%`.
- **`S₄` pairs-and-solve, `m = 3`: relabelling.**  Same base, same
  trials to within one, same yield, same elimination; the oracle
  replaces the table's `|F|²/2` entries by a per-target loop over
  `|F|²/2` pairs, each pair a quartic solve.  Memory falls to `O(1)` and
  `S` rises `19×` at `n = 15` and `510×` at `n = 27`; the count that
  fell is the table, and the headline number is not the count.  It is
  on the board because it is the only algebraic oracle that finishes a
  logarithm at these sizes, and its exponent (`1.06`, §4) is the price
  of having no table.

### 3.4 Koblitz curve: the fold is real, small, and not where the cost is

Two column maps on one base, so the trials and the oracle work are
identical by construction and only the number of unknowns differs.

- **One column per abscissa (control, baseline).**  The random-binary
  pipeline run on the Koblitz curve, against the signed-Frobenius
  reference that quotients by all `2n` automorphisms.  `227×` the
  reference at `n = 23`.
- **Signed-orbit columns: advance, count.**  `K` falls by the orbit
  length: `438 → 20` at `n = 23`, `1,211 → 41` at `n = 31`, and the
  trials fall with it (`348 → 18`, `951 → 49`) because the elimination
  pins the logarithm after `K + 1` relations.  That is the Gaudry /
  Galbraith–Menezes–Prasad fold reproduced with every phase counted,
  and it is a count advance: a generic algorithm has no orbits to fold.
  It is also small in `S` — `7%` at `n = 23`, `27%` at `n = 31`, `65%` at
  `n = 15` — because the phase the fold shrinks was already small.  At
  `n = 23` the factor base is `99.6%` of the total: `875` signed points
  need `438` Artin–Schreier solves and the pair table `|F|²/2 ≈
  383,000` additions, against `18` trials.  The fold cannot move that.
- **`S₄` pairs-and-solve on the invariant subspace: relabelling**, for
  the reason of §3.3: `128,758` against `2,489` at `n = 31`.
- **`n = 41`.**  The largest instance on the board, `r = 2^39`, cofactor
  4.  The base is capped at `6,000` signed points (`max_table_points`),
  which is dimension 12 rather than the `⌈n/3⌉ = 14` the rule would
  give, so `p_ceiling` is `0.0095` and `6,086` trials are needed: here
  the relation phase is `71%` and the table `29%`, and `S = 58.8` is
  `300×` the signed-Frobenius reference and `425×` the floor.  It is
  lower than the `n = 23` row's `188` because the base did not grow
  with `n`, not because anything improved; the same cap would make it
  rise again at `n = 47`.
- **The large-cofactor rungs** (`n = 29, 31, 37, 39`, cofactors `596` to
  `12,646`) have a subgroup far smaller than the field, so a base sized
  by `n` is enormous relative to `√r` and `S` runs to `34,554`.  They
  are honest rows (every logarithm verified), they are what a Koblitz
  curve with a large cofactor costs, and they are why the fold's
  exponent fit in §4 has `R² = 0.71`.

### 3.5 An accounting correction: the counting ceiling is loose on confined bases

Six rows carry `yield/ceiling` above `1.4` (marked † in §2.2), and
none of them is a yield finding.  The ceiling of §1.3 spreads the
`C(F+m−1, m)` sums over all of `#E`; if every `m`-sum lies in a
proper subgroup `H`, they are spread over `|H|` and the true ceiling is
larger by the index `[E : H]`.  That is the case on every marked row:

- `random-binary-n18` and `-n24` (`a = 0`, even `n`): every element of
  the low-order subspace `⟨1, z, …, z^{l−1}⟩` has trace zero (`Tr(1) = n
  mod 2 = 0` and `Tr(z^i) = 0` for `i < l` on both moduli), so every base
  point satisfies `Tr(x) = Tr(a)` and lies in `2E`, an index-2
  subgroup; all `32` and all `137` base abscissae, checked.  Corrected,
  the rows read `0.80` and `0.87`.
- `K_1 / GF(2^17)` (`a = 1`, `Tr(a) = 1`): a Frobenius-invariant subspace
  that does not contain `F₂` lies in the trace-zero hyperplane, so the
  base points lie *outside* `2E`, no sum of three of them can reach the
  odd-order subgroup (the census found `m = 3` yields nothing, which is
  why this rung runs at `m = 2`), and every sum of two lies in `2E`.
  Corrected by the same index 2: `0.78`.
- `K_0 / GF(2^15)` and `K_0 / GF(2^39)` (`a = 0`, cofactors `44` and
  `8,012`): the base is again in `2E`, which accounts for a factor 2;
  the rest is the finer class structure of a large cofactor and, at
  `n = 15`, the variance of two or three relations per run.

This is an **accounting** correction to the boundary, not a change to
any algorithm or any count, and it moves no `S`.  The next round should
tighten `decomposition_probability_ceiling` to divide by the order of
the subgroup generated by the base's `m`-sums (the trace condition
gives the index-2 part in closed form; the cofactor classes the
pipeline already computes for admissibility give the rest); until then
the † rows are to be read with the corrected values above.

### 3.6 What the ledger replaces

`docs/ic/BOUNDARY_TARGETS.md` carried, for the Koblitz regime, an
"online charged crossover" at `n = 41` (`0.29×` rho on the online phase,
`5.0×` amortised, `111×` with the full build), and whole-process
*wall* wins at `n = 53` on one host.  Those rows stay: they are
wall-clock measurements of a different stack (a materialised
point-defined base with an exact support table).  What this ledger adds
next to them is the one number that survives hardware: with every phase
counted, on this stack, `n = 41` is `300×` the signed-Frobenius
reference.  The prime and binary regimes had no operation-counted rows
at all; they have the ladders above now, with their `relation_yield`,
`rank` and `end_to_end_dlp` stages filled from the same run.

## 4. Every phase priced: the exponents

Least-squares fits of `log(operations)` against `log r` (and against
`log #E`) over the ladder, per variant and per phase, from the report's
`fits` array.  `trials` is the target count, not a cost.  Rho's `0.5`
is the mark; its *measured* exponent here is lower because its setup
is a constant at these sizes (§2.3), and every variant is above both.

| regime | variant | phase | α (r) | R² | α (#E) | R² | sizes |
|:--|:--|:--|--:|--:|--:|--:|--:|
| prime | rho_reference | total | 0.274 | 0.903 | 0.274 | 0.903 | 8 |
| prime | semaev_s3_roots_m2 | total | 0.905 | 0.959 | 0.905 | 0.959 | 8 |
| prime | semaev_s3_roots_m2 | relations | 0.907 | 0.959 | 0.907 | 0.959 | 8 |
| prime | semaev_s3_roots_m2 | factor_base | 0.435 | 0.852 | 0.435 | 0.852 | 8 |
| prime | semaev_s3_roots_m2 | linear_algebra | 0.203 | 0.740 | 0.203 | 0.740 | 8 |
| prime | semaev_s3_roots_m2 | trials | 0.568 | 0.989 | 0.568 | 0.989 | 8 |
| prime | direct_subtraction_m2 | total | 0.850 | 0.995 | 0.850 | 0.995 | 8 |
| prime | direct_subtraction_m2 | relations | 0.852 | 0.995 | 0.852 | 0.995 | 8 |
| prime | direct_subtraction_m2 | factor_base | 0.435 | 0.852 | 0.435 | 0.852 | 8 |
| prime | direct_subtraction_m2 | linear_algebra | 0.203 | 0.740 | 0.203 | 0.740 | 8 |
| prime | direct_subtraction_m2 | trials | 0.568 | 0.989 | 0.568 | 0.989 | 8 |
| prime | mitm_m2 | total | 0.662 | 0.996 | 0.662 | 0.996 | 8 |
| prime | mitm_m2 | relations | 0.670 | 0.994 | 0.670 | 0.994 | 8 |
| prime | mitm_m2 | factor_base | 0.647 | 0.961 | 0.647 | 0.961 | 8 |
| prime | mitm_m2 | linear_algebra | 0.203 | 0.740 | 0.203 | 0.740 | 8 |
| prime | mitm_m2 | trials | 0.568 | 0.989 | 0.568 | 0.989 | 8 |
| prime | mitm_m3 | total | **0.601** | 0.982 | 0.601 | 0.982 | 8 |
| prime | mitm_m3 | relations | 0.511 | 0.993 | 0.511 | 0.993 | 8 |
| prime | mitm_m3 | factor_base | 0.647 | 0.961 | 0.647 | 0.961 | 8 |
| prime | mitm_m3 | linear_algebra | 0.623 | 0.964 | 0.623 | 0.964 | 8 |
| prime | mitm_m3 | trials | 0.329 | 0.988 | 0.329 | 0.988 | 8 |
| char2 | rho_reference | total | 0.302 | 0.973 | 0.271 | 0.936 | 5 |
| char2 | mitm_m3 | total | **0.624** | 0.986 | 0.564 | 0.963 | 5 |
| char2 | mitm_m3 | relations | 0.640 | 0.983 | 0.571 | 0.935 | 5 |
| char2 | mitm_m3 | factor_base | 0.578 | 0.953 | 0.541 | 0.998 | 5 |
| char2 | mitm_m3 | linear_algebra | 0.706 | 0.960 | 0.651 | 0.976 | 5 |
| char2 | mitm_m3 | trials | 0.366 | 0.957 | 0.318 | 0.864 | 5 |
| char2 | semaev_s4_pairs_and_solve_m3 | total | 1.058 | 0.989 | 0.960 | 0.974 | 5 |
| char2 | semaev_s4_pairs_and_solve_m3 | relations | 1.060 | 0.989 | 0.962 | 0.974 | 5 |
| char2 | semaev_s4_pairs_and_solve_m3 | factor_base | 0.412 | 0.966 | 0.383 | 0.998 | 5 |
| char2 | semaev_s4_pairs_and_solve_m3 | linear_algebra | 0.704 | 0.962 | 0.645 | 0.966 | 5 |
| char2 | semaev_s4_pairs_and_solve_m3 | trials | 0.357 | 0.946 | 0.309 | 0.847 | 5 |
| koblitz | rho_reference | total | 0.254 | 0.917 | 0.171 | 0.627 | 11 |
| koblitz | mitm_m3_signed_orbit_columns | total | **0.544** | 0.709 | 0.520 | 0.933 | 10 |
| koblitz | mitm_m3_signed_orbit_columns | relations | 0.646 | 0.830 | 0.571 | 0.937 | 10 |
| koblitz | mitm_m3_signed_orbit_columns | factor_base | 0.508 | 0.643 | 0.499 | 0.894 | 10 |
| koblitz | mitm_m3_signed_orbit_columns | linear_algebra | 0.214 | 0.523 | 0.219 | 0.787 | 10 |
| koblitz | mitm_m3_signed_orbit_columns | trials | 0.357 | 0.853 | 0.307 | 0.907 | 10 |
| koblitz | mitm_m3_abscissa_columns_control | total | 0.753 | 0.593 | 0.632 | 0.967 | 7 |
| koblitz | mitm_m3_abscissa_columns_control | relations | 0.585 | 0.657 | 0.444 | 0.873 | 7 |
| koblitz | mitm_m3_abscissa_columns_control | factor_base | 0.845 | 0.621 | 0.676 | 0.917 | 7 |
| koblitz | mitm_m3_abscissa_columns_control | linear_algebra | 0.802 | 0.833 | 0.430 | 0.553 | 7 |
| koblitz | mitm_m3_abscissa_columns_control | trials | 0.356 | 0.811 | 0.245 | 0.885 | 7 |

### 4.1 Which phase dominates, and where

| regime, variant | phase that dominates | its share, smallest → largest instance | what it is |
|:--|:--|:--|:--|
| prime, `S₃` roots and direct subtraction | relations | `97% → 99.99%` | `|F| ∝ r^{1/3}` roots or subtractions per target, `r^{0.57}` targets: `r^{0.9}` |
| prime, meet in the middle `m = 2` | relations | `57% → 83%` | `r^{0.57}` targets at `65` additions each, the two scalar multiplications that make a target |
| prime, meet in the middle `m = 3` | pair table (factor base) | `57% → 71%` (relations at 12 bits) | `|F|²/2 ∝ r^{2/3}` additions, once |
| binary, meet in the middle `m = 3` | relations | `77% → 82%` | `|F|` subtractions per target, `r^{0.37}` targets |
| binary, `S₄` pairs-and-solve | relations | `99.1% → 99.99%` | `|F|²/2` quartic solves per target |
| Koblitz, signed-orbit columns | factor base, until the cap | `85–99.6%` on `n = 11…31` (56% at `n = 15`), `72%` at `n = 39`, `29%` at `n = 41` | Artin–Schreier lifts and the pair table; at `n = 41` the capped base pushes `71%` into relations |
| Koblitz, abscissa columns (control) | factor base | `61% → 92%` (`n = 11…23`); linear algebra rises to `18.7` at `n = 31` | as above, plus a `K ≈ |F|/2` elimination |

Linear algebra is under `1%` of every row except the Koblitz control at
`n = 31` (`0.5%`) and is not the phase to price further at these sizes;
its fitted `0.62–0.80` on the dense elimination says where it would go.

### 4.2 The Koblitz fit, restricted (derived, not in the report)

The report fits every rung.  Re-fitting the frozen rows by the same
least squares over only the rungs with cofactor `≤ 4` (`n = 11, 13, 19,
23, 41` for the `m = 3` fold; `n = 17` runs at `m = 2` and is excluded;
the control stops at `n = 31`) gives a cleaner picture of the two
column maps on comparable instances.  These are derived from the
report's per-instance counts and are not in its `fits` array:

| variant | phase | α (r) | R² | rungs |
|:--|:--|--:|--:|--:|
| signed-orbit columns | total | 0.509 | 0.981 | 5 |
| signed-orbit columns | relations | 0.615 | 0.938 | 5 |
| signed-orbit columns | factor_base | 0.452 | 0.950 | 5 |
| signed-orbit columns | linear_algebra | 0.175 | 0.975 | 5 |
| abscissa columns (control) | total | 0.624 | 0.992 | 4 |
| abscissa columns (control) | relations | 0.455 | 1.000 | 4 |
| abscissa columns (control) | factor_base | 0.664 | 0.988 | 4 |
| abscissa columns (control) | linear_algebra | 0.617 | 0.994 | 4 |
| signed-Frobenius rho | total | 0.280 | 0.948 | 6 |

The fold's `0.51` against the control's `0.62` is the count advance of
§3.4 read as an exponent: the elimination falls from `0.62` to `0.18`
and the trials with it.  The factor-base exponent `0.45` is depressed
by the `n = 41` cap (a base of dimension `⌈n/3⌉` would put the pair
table at `r^{2/3}`), so the `0.51` is optimistic and the honest
statement is: at or above one half, with the table as the term that
sets it.

### 4.3 Extrapolation, marked as such

Against rho's limiting one half, every total exponent on the board is
higher, so every gap widens: the prime `m = 3` row, `14.5×` at `2^23.4`
with `α = 0.60`, extrapolates to about `46×` at `2^40` and `240×` at
`2^64`; the binary row, `89×` at `2^24.4` with `0.62`, to about `340×`
at `2^40`.  These rest on the fitted exponents above, over eight and
five sizes, and on the reference's one half rather than its measured
`0.27–0.30`, which would make them worse.  There is no size at which
the fits cross.

## 5. The decomposition oracles, per target

The pipelines above use the oracles that finish a logarithm at these
sizes.  The algebraic ones — matrix-F4 with splitting over `F₂` and
CDCL with native parity rows — do not, so they are priced per target on
the same Semaev systems (`ic boundary --oracles`, the report's
`oracle_pricing` field): a Koblitz curve `K_a / F_{2^n}`, its
Frobenius-invariant subspace base (or the orbit-union base where the
subspace is the subfield), eight targets per cell drawn as `[a]G +
[b]Q`, the Weil-descended symmetrised `S₃` (`m = 2`) or chained `S₄`
(`m = 3`) system in `m·l + (m−2)·n` unknowns, and every oracle asked
the same question of the same target.  `disagreements` counts targets
on which any two oracles that returned an answer differed; it is zero
in every cell (`all_agree` is true), and every witness sums to its
target in the group.

### 5.1 The systems

`FFD` is the first fall degree over the eight draws (the smallest
Macaulay degree with a non-trivial syzygy); the Macaulay columns give
rows × columns and the rank at degrees 2 and 3.

| n | dim | m | \|F\| | K | unknowns | eqs | deg | eq/var | FFD | Macaulay D=2 rows×cols (rank) | D=3 | hit rate | disagreements |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|:--|--:|--:|
| 9 | 6 | 2 | 55 | 4 | 12 | 9 | 2 | 0.75 | 2 | 9×49 (8) | 117×259 (104) | 1.00 | 0 |
| 9 | 6 | 3 | 55 | 4 | 27 | 18 | 3 | 0.67 | 3 | 9×70 (9) | 261×1672 (260) | 1.00 | 0 |
| 11 | 10 | 2 | 991 | 46 | 20 | 11 | 2 | 0.55 | 2 | 11×121 (10) | 231×1111 (210) | 1.00 | 0 |
| 13 | 12 | 2 | 4005 | 155 | 24 | 13 | 2 | 0.54 | 2 | 13×169 (12) | 325×1885 (300) | 1.00 | 0 |
| 13 | 12 | 3 | 4005 | 155 | 49 | 26 | 3 | 0.53 | not profiled (49 > 48-unknown cap) | — | — | 1.00 | 0 |
| 15 | 5 | 2 | 33 | 3 | 10 | 15 | 2 | 1.50 | 2–3 | 15×36 (15) | 165×156 (136) | 0.00 | 0 |
| 15 | 5 | 3 | 33 | 3 | 30 | 30 | 3 | 1.00 | 3 | 15×96 (15) | 480×2246 (479) | 0.38 | 0 |
| 17 | 8 | 2 | 239 | 8 | 16 | 17 | 2 | 1.06 | 2 | 17×81 (16) | 289×585 (272) | 0.25 | 0 |
| 23 | 11 | 2 | 2071 | 46 | 22 | 23 | 2 | 1.05 | 2 | 23×144 (22) | 529×1464 (506) | 0.50 | 0 |

The first fall degree is `2` for every `m = 2` system and `3` for
every `m = 3` one, which is the chained-Semaev ladder of
`research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md` (no fall later than 3) reproduced
on these draws; it does not move with `n` between 9 and 23.  The hit
rate is 1 while the base's sums cover the group (`n ≤ 13`, where
`C(F+m−1, m) ≫ #E`) and falls to `0–0.5` once they do not.

### 5.2 What one target costs

`native` is the oracle's own exact count (median over the targets it
found, and over those it refuted); `GAE/target` converts it at the
host factor and averages over all eight; `projected S` is
`GAE/target × trials_floor / √r` with `trials_floor = (K+1)/p_ceiling`,
an **extrapolation** of what the relation phase alone would cost if
this oracle drove the pipeline, with no table and no factor base in it.
`vs floor` divides that by `√(π/4n)`.

| n | m | oracle | unit | found / refuted / no answer | native, found | native, refuted | ms, found | ms, refuted | GAE/target | GAE/refutation | projected S | vs floor |
|--:|--:|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 9 | 2 | enumerate | group additions | 8/0/0 | 3.00 | — | 0.001 | — | 5.75 | — | 2.55 | 8.6 |
| 9 | 2 | meet in the middle | pair-table probes | 8/0/0 | 1.00 | — | 0.000 | — | 1.02 | — | 0.45 | 1.5 |
| 9 | 2 | matrix-F4, splitting | word XORs | 8/0/0 | 46,150 | — | 1.885 | — | 1.70e2 | — | 75.4 | 255 |
| 9 | 2 | CDCL, native XOR | conflicts | 8/0/0 | 60.0 | — | 0.586 | — | 4.12e3 | — | 1.83e3 | 6.19e3 |
| 9 | 3 | enumerate | group additions | 8/0/0 | 3.50 | — | 0.002 | — | 5.00 | — | 2.22 | 7.5 |
| 9 | 3 | meet in the middle | pair-table probes | 8/0/0 | 1.00 | — | 0.001 | — | 2.02 | — | 0.90 | 3.0 |
| 9 | 3 | `S₄` pairs-and-solve | pairs | 8/0/0 | 64.0 | — | 0.032 | — | 2.68e2 | — | 119 | 402 |
| 9 | 3 | matrix-F4, splitting | word XORs | 8/0/0 | 2,041,484 | — | 27.2 | — | 7.53e3 | — | 3.34e3 | 1.13e4 |
| 9 | 3 | CDCL, native XOR | conflicts | 8/0/0 | 409 | — | 9.167 | — | 7.25e4 | — | 3.22e4 | 1.09e5 |
| 11 | 2 | enumerate | group additions | 8/0/0 | 1.00 | — | 0.001 | — | 1.00 | — | 1.49 | 5.6 |
| 11 | 2 | meet in the middle | pair-table probes | 8/0/0 | 1.00 | — | 0.001 | — | 1.02 | — | 1.52 | 5.7 |
| 11 | 2 | matrix-F4, splitting | word XORs | 8/0/0 | 677,942 | — | 10.3 | — | 1.84e3 | — | 2.75e3 | 1.03e4 |
| 11 | 2 | CDCL, native XOR | conflicts | 8/0/0 | 284 | — | 2.678 | — | 1.52e4 | — | 2.27e4 | 8.48e4 |
| 13 | 2 | enumerate | group additions | 8/0/0 | 1.00 | — | 0.002 | — | 1.00 | — | 3.49 | 14.2 |
| 13 | 2 | meet in the middle | pair-table probes | 8/0/0 | 1.00 | — | 0.002 | — | 1.02 | — | 3.55 | 14.4 |
| 13 | 2 | matrix-F4, splitting | word XORs | 8/0/0 | 2,731,272 | — | 26.8 | — | 6.39e3 | — | 2.23e4 | 9.07e4 |
| 13 | 2 | CDCL, native XOR | conflicts | 8/0/0 | 222 | — | 2.858 | — | 2.46e4 | — | 8.56e4 | 3.48e5 |
| 13 | 3 | enumerate | group additions | 8/0/0 | 2.00 | — | 0.001 | — | 2.00 | — | 6.97 | 28.4 |
| 13 | 3 | meet in the middle | pair-table probes | 8/0/0 | 1.00 | — | 0.002 | — | 2.02 | — | 7.03 | 28.6 |
| 13 | 3 | `S₄` pairs-and-solve | pairs | 8/0/0 | 4,096 | — | 0.362 | — | 2.41e4 | — | 8.41e4 | 3.42e5 |
| 15 | 2 | enumerate | group additions | 0/8/0 | — | 33.0 | — | 0.007 | 33.0 | 33.0 | 154 | 674 |
| 15 | 2 | meet in the middle | pair-table probes | 0/8/0 | — | 1.00 | — | 0.000 | 0.011 | 0.011 | 0.051 | 0.22 |
| 15 | 2 | matrix-F4, splitting | word XORs | 0/8/0 | — | 17,925 | — | 0.607 | 32.6 | 32.6 | 152 | 666 |
| 15 | 2 | CDCL, native XOR | conflicts | 0/8/0 | — | 5.00 | — | 0.387 | 1.54e3 | 1.54e3 | 7.18e3 | 3.14e4 |
| 15 | 3 | enumerate | group additions | 3/5/0 | 229 | 594 | 0.055 | 0.133 | 4.67e2 | 5.94e2 | 182 | 794 |
| 15 | 3 | meet in the middle | pair-table probes | 3/5/0 | 8.00 | 33.0 | 0.003 | 0.010 | 24.7 | 33.4 | 9.60 | 42.0 |
| 15 | 3 | `S₄` pairs-and-solve | pairs | 3/5/0 | 203 | 528 | 0.152 | 0.369 | 1.12e3 | 1.44e3 | 435 | 1.90e3 |
| 15 | 3 | matrix-F4, splitting | word XORs | 3/5/0 | 171,194,633 | 1,698,176,490 | 859 | 8,255 | 2.37e6 | 3.60e6 | 9.22e5 | 4.03e6 |
| 15 | 3 | CDCL, native XOR | conflicts | 0/0/8 | — | — | — | — | 3.97e7 (at the 200,000-conflict budget) | — | 1.55e7 | 6.76e7 |
| 17 | 2 | enumerate | group additions | 2/6/0 | 15.5 | 239 | 0.006 | 0.062 | 1.83e2 | 2.39e2 | 25.7 | 120 |
| 17 | 2 | meet in the middle | pair-table probes | 2/6/0 | 1.00 | 1.00 | 0.002 | 0.002 | 0.26 | 0.011 | 0.037 | 0.17 |
| 17 | 2 | matrix-F4, splitting | word XORs | 2/6/0 | 3,042,354 | 4,979,528 | 36.6 | 60.0 | 8.36e3 | 9.26e3 | 1.18e3 | 5.47e3 |
| 17 | 2 | CDCL, native XOR | conflicts | 2/6/0 | 944 | 2,052 | 9.871 | 22.3 | 7.21e4 | 8.38e4 | 1.01e4 | 4.72e4 |
| 23 | 2 | enumerate | group additions | 4/4/0 | 468 | 2,071 | 0.179 | 0.766 | 1.31e3 | 2.07e3 | 60.1 | 325 |
| 23 | 2 | meet in the middle | pair-table probes | 4/4/0 | 1.00 | 1.00 | 0.003 | 0.002 | 0.51 | 0.011 | 0.023 | 0.13 |
| 23 | 2 | matrix-F4, splitting | word XORs | 4/4/0 | 68,297,378 | 226,738,826 | 406 | 1,352 | 2.00e5 | 3.14e5 | 9.16e3 | 4.96e4 |
| 23 | 2 | CDCL, native XOR | conflicts | 4/4/0 | 29,247 | 188,443 | 561 | 9,529 | 1.29e7 | 2.45e7 | 5.91e5 | 3.20e6 |

Engine detail, from the same cells: matrix-F4 at `n = 23, m = 2` made
`2,588` reductions over `1,296` splits with `1,283` infeasible branches
and never needed a Macaulay degree above the cap; at `n = 15, m = 3` it
made `8,937` reductions over `3,631` splits.  The CDCL solver's parity
rows contributed `176` implied rows at `n = 23` and it returned no
spurious model in any cell (every model was checked in the group).

### 5.3 Reading the oracle ledger

- **The order of the oracles is the order of the ledger, at every
  size.**  Per target at `n = 23, m = 2`: one probe (`0.5` additions),
  enumeration (`1,310`), matrix-F4 (`2.0 × 10⁵`), CDCL (`1.3 × 10⁷`).  At
  `n = 15, m = 3`, `30` unknowns: probes (`25`), enumeration (`467`),
  `S₄` pairs-and-solve (`1,118`), matrix-F4 (`2.4 × 10⁶`), CDCL with no
  answer on any of eight targets at the budget.  Five orders of
  magnitude separate the tabled probe from the best algebraic solver on
  the same system.
- **Refuting costs more than finding**, for the solvers that must
  exhaust: matrix-F4 `3.3×` more at `n = 23` (`2.3 × 10⁸` word XORs
  against `6.8 × 10⁷`) and `10×` at `n = 15, m = 3`; CDCL `6.4×` at
  `n = 23` (`188,000` conflicts against `29,000`).  A relation search
  refutes most of its targets (the hit rate is `0.5` at `n = 23`, and
  `p_ceiling` at a useful base size is well under one), so the
  refutation price is the one a pipeline pays.
- **The projected relation phase for a tabled probe is under the
  floor** (`0.023` against `0.127` at `n = 23`), and that is the caveat
  of §1.6 made numeric: the probe is cheap because the `|F|²/2` pair
  table already did the work, and on the pipeline row for the same
  degree (§2.1, `K_1 / GF(2^23)`) the table is `99.6%` of `S = 188`.  A
  per-target price is not a pipeline cost.
- **Driving the pipeline with matrix-F4 instead of the table** would put
  the relation phase alone at a projected `S ≈ 9.2 × 10³` at `n = 23`
  (`5 × 10⁴` times the floor), `49×` the *whole* measured cost of the
  tabled fold on that degree, and CDCL at `5.9 × 10⁵`.  The F4
  refutation cost grows `46×` from `n = 17` to `n = 23` at `m = 2`,
  two points and so not a fit; it is enough to say that the algebraic
  route does not reach the table at any size on this ladder.
- **Where the algebraic oracles stand as a target for the next round:**
  the number to beat at `n = 23, m = 2` is `2.0 × 10⁵` additions per
  target (`6.8 × 10⁷` word XORs found, `2.3 × 10⁸` refuted), with
  `FFD = 2`, `22` unknowns and `23` equations, agreeing with the
  complete oracle on every target; at `n = 15, m = 3`, `2.4 × 10⁶` with
  `FFD = 3` and `30` unknowns.  The table's price on the same systems is
  the boundary they would have to cross to matter.

### 5.4 Re-run on the merged tree, with the lift-checked tally

A review finding on the first CI round was correct: the `S₄` tally
behind §5.2 counted a disagreement only when the oracle refuted a target
the exhaustive search could decompose, so a spurious witness would not
have counted.  The pricing now lifts every `S₄` witness over the signed
base with the pipeline's own sign-lift, charges the lift to the row,
records lift failures, and tests agreement in both directions.  The
thirteen cells were re-run under that tally on the merged tree (commit
`7b2e1b5`, `docs/ic/runs/ic-oracle-pricing-lifted-2026-09-21.json`, the
same seed and targets as the frozen run):

- every verdict, hit rate, first fall degree and SAT conflict count is
  identical to §5.1 and §5.2, `disagreements` is zero in every cell, and
  no `S₄` witness failed to lift (`lift_failures = 0` on all three
  `m = 3` cells); the `S₄` prices move only with the host factor
  (`284` against `268` additions per target at `n = 9`);
- matrix-F4 does `1.4×` to `2.7×` fewer word XORs than the frozen run,
  because the merge brought in main's F4 work since the freeze: the
  sparse structured elimination for binary Macaulay matrices (#405) and
  the native-elimination speed-up (5fa6138f).  Same verdicts, same
  splits, fewer XORs per reduction: **engineering** by the §3 test, in
  main's solver rather than in anything this note built, and the ratio
  to the tabled probe is unmoved at five orders of magnitude.

| n | m | F4 word XORs, found (frozen → re-run) | refuted (frozen → re-run) | GAE per target | projected relation-phase S | total XORs, ratio |
|--:|--:|:--|:--|:--|:--|--:|
| 9 | 2 | 46,150 → 32,114 | — | 170 → 167 | 75.4 → 74.2 | 1.44× |
| 9 | 3 | 2,041,484 → 1,057,972 | — | 7.53e3 → 3.80e3 | 3.34e3 → 1.69e3 | 1.98× |
| 11 | 2 | 677,942 → 450,018 | — | 1.84e3 → 1.25e3 | 2.75e3 → 1.87e3 | 1.52× |
| 13 | 2 | 2,731,272 → 1,795,996 | — | 6.39e3 → 4.17e3 | 2.23e4 → 1.45e4 | 1.52× |
| 15 | 2 | — | 17,925 → 11,329 | 32.6 → 20.8 | 152 → 97.1 | 1.57× |
| 15 | 3 | 171,194,633 → 64,818,554 | 1,698,176,490 → 628,830,074 | 2.37e6 → 8.73e5 | 9.22e5 → 3.40e5 | 2.70× |
| 17 | 2 | 3,042,354 → 1,647,510 | 4,979,528 → 2,683,981 | 8.36e3 → 4.47e3 | 1.18e3 → 628 | 1.86× |
| 23 | 2 | 68,297,378 → 32,153,918 | 226,738,826 → 105,986,387 | 2.00e5 → 8.84e4 | 9.16e3 → 4.06e3 | 2.14× |

The GAE column moves less than the XOR column on the small cells because
a re-run is a new calibration and the host factors of the two runs
differ; the native counts are the measurement.  The numbers to beat in
§5.3 therefore read, on the merged tree: `8.8 × 10⁴` additions per
target at `n = 23, m = 2` (`3.2 × 10⁷` word XORs found, `1.06 × 10⁸`
refuted) and `8.7 × 10⁵` at `n = 15, m = 3`; the projected F4-driven
relation phase at `n = 23` is `S ≈ 4.1 × 10³`, `22×` the whole measured
cost of the tabled fold.  §5.2 stays as the measurement at the frozen
commit; this section is the current one, and the scoreboard's F4 rows
carry both.

## 6. The algorithms on the shelf, and where each one is priced

The user's question was also an inventory question: which solvers and
which polynomial systems exist here, and what does each cost.  In one
place:

| what | where | what it does | priced in this note as |
|:--|:--|:--|:--|
| Semaev `S₃` (prime and binary) | `ec_index_calculus::semaev_s3_in_x3`, `binary_semaev::binary_semaev_s3`, `ic_boundary::semaev_s3_quadratic` | the summation polynomial as a quadratic in one abscissa; two-point decomposition by root finding | prime `semaev_s3_roots_m2` |
| symmetrised `S₄` (binary, any `b`) | `semaev_decomp::SubspaceOracle` (this round: general `b`, general subspace basis), `binary_semaev_s4::weil_descend_s4` | quartic in the third abscissa; pairs-and-solve with the subspace polynomial `L_V`; Weil descent to a Boolean system | `char2` and `koblitz` `semaev_s4_pairs_and_solve_m3`; oracle cells |
| `S₅`, `S₆` (prime) | `semaev_higher` | Semaev's recurrence and a Sylvester resultant; degrees `8`, `16` in the last variable | not priced here: no pipeline uses them at these sizes |
| symmetrised summation polynomials (Faugère–Gaudry–Huot–Renault) | `symmetrized_semaev`, `koblitz_symmetrised` | elementary-symmetric change of variables; `S₃` monomials `17 → ≤ 10` | oracle cells use the symmetrised binary `S₄` |
| matrix-F4 over `F₂` with splitting | `koblitz_groebner::matrix_f4_f2_counted`, `solve_boolean_system` | Macaulay matrices at degrees `2..3`, bit-packed elimination, propagation, splitting on stall; exact 64-bit word XOR count (process-wide counter added this round) | oracle cells `matrix_f4_splitting` |
| Buchberger over `F₂` | `pq_groebner_f2::groebner_basis_f2` | the reference engine the F4 is tested against | not priced (same answers, far slower) |
| F4 over `F_p` | `f4_fp::f4`, `f4_fp::solve` | degree-bounded F4 with dense modular elimination | not on a pipeline here; used by the `E(F_{p³})` Gaudry work (`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md`) |
| XL, crossbred | `pq_xl`, `crossbred` | linearisation at a degree bound; guess-and-linearise | not priced here |
| CDCL SAT with native parity rows | `sat::Solver`, `semaev_sat` | Gauss–Jordan on XOR rows at the propagation fixpoint, Macaulay rows as implied clauses | oracle cells `cdcl_sat_native_xor` |
| meet in the middle (pair table) | `ic_boundary::PairTable` (this note), `koblitz_index_calculus::PairSumTable` (the production one, with compact and folded tiers) | all pair sums keyed by abscissa; one probe (`m = 2`) or `|F|` subtractions and probes (`m = 3`) | every regime |
| enumeration | `ic_oracle_pricing::enumerate_counted`, `koblitz_index_calculus::enumerate_decompose` | `|F|^{m−1}` group additions | oracle cells `enumerate` |
| dense Gauss–Jordan over `Z/rZ` | `ic_boundary::IncrementalGauss` | incremental reduced echelon form, multiply-subtracts counted | every regime (`linear_algebra` phase) |
| relation filtering + block Wiedemann | `koblitz_sparse_la` | singleton and clique removal, structured elimination, Krylov sequence with a matrix Berlekamp–Massey | not needed at these matrix sizes; measured at `n^{0.68}` on `E(F_{p³})` in `research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` |
| Nagao / Riemann–Roch encodings | `research/nagao_relations/`, `ecc2k130/codegen/nagao*.py` | incidence and norm encodings without a summation polynomial; function-first support test `H | L_V` | not priced: a correctness panel only, no full-attack `S` (their README says so) |
| first fall degree | `koblitz_groebner::first_fall_degree`, `ffd_harness` | smallest `D` whose Macaulay matrix has a non-trivial syzygy | every oracle cell |

## 7. The corpus for external solvers

`ic corpus` writes the family of instances the Trimoska–Ionica–Dequen
generator produces — the Weil-descended symmetrised `S₄` of
`y² + xy = x³ + x² + 1` over `F_{2^n}` with the three abscissae confined
to `⟨1, z, …, z^{l−1}⟩` — in four formats sharing one variable numbering:
DIMACS with native `x` parity lines (the CryptoMiniSat / WDSat
convention), plain CNF, one GF(2) polynomial per line, and a Magma
script computing the Gröbner basis; plus an `INFO` file with the
target, the label, the planted witness and its SAT assignment.

Two facts about it worth knowing.  Our encoder, written from the
algebra rather than from upstream's code, produces upstream's instance
size exactly at `n = 19, l = 6`: **767 variables, 2 364 clauses, 52
parity rows** against their `p cnf 767 2416` with 52 `x`-lines
(`docs/ic/corpus/n19l6.json`).  And the labels are certified: a
satisfiable instance plants a sum of three factor-base points and is
confirmed by a model of this crate's solver, an unsatisfiable one is
refuted by the complete pairs-and-solve search — upstream's
`n19l6-19-U` is satisfiable despite its name (`research/notes/ecc2k130/RESEARCH_TRIMOSKA_BENCHMARKS.md`),
which is the failure mode the certification exists for.

## 8. What does not count, and what this is not

- Not a claim about a deployed curve.  The largest instance is a
  40-bit Koblitz subgroup; the point is a concrete boundary to iterate
  against, measured rather than extrapolated.
- Not a wall-clock benchmark.  The wall columns are practicality
  notes; the metric is the count, and the conversion factors are in
  the report so the count can be re-priced on other hardware.
- Not best-of-breed composition.  Every row is one process on one
  instance with every phase inside it; the oracle cells' projected
  relation phases are extrapolations and are labelled as such.
- Not a speedup claim.  Nothing here is faster than anything; the round
  establishes the baseline table and the phase exponents, and by
  `AGENTS.md` §8 no gain against them can be claimed without the
  frozen regression suite and a matched full-DLP comparison alongside.
- The Koblitz `S_floor` uses all `2n` automorphisms; the reference
  walk uses them too.  The prime and random-binary references use none,
  so their `vs rho` column is generous to index calculus by the `√2`
  the negation map would buy; the `vs floor` column is not.

## 9. Reproducing

```bash
cargo build --release --bin ic
./target/release/ic boundary --quick                 # a few seconds, every regime
./target/release/ic boundary --oracles --out out.json # the frozen run, about an hour
./target/release/ic corpus --degree 19 --dimension 6 --sat 5 --unsat 5 --dir corpus/n19l6
cargo test --release --lib cryptanalysis::ic_       # oracle agreement, planted logs, fits
cargo test --release --test ic_framework            # the CLI end to end
```

The harness is deterministic in its seeds: the same seed gives the same
curves, targets, trials and counts on any host; only the nanosecond
factors and the wall columns move.

The tables in §2, §4 and §5 and the rows of the scoreboard panel are
rendered from the frozen report by `docs/ic/tools/boundary_ledger_tables.py`
and `docs/ic/tools/boundary_scoreboard_rows.py`; the ledger twins were
updated from it by `docs/ic/tools/boundary_ledger_update.py`.  Re-run
them on a new report rather than editing numbers by hand.

## 10. Round 2: an engineering ledger, and one artifact caught

**Frozen run:** `docs/ic/runs/ic-boundary-ledger-round2-2026-09-21.json`
(`ic boundary --out …`, status `complete`, `all_verified` true, 2,026 s
on the same host as §2, at commit `4d0e9f9b`; 24 instances, 498 rows.
The oracle cells of §5 are unchanged and stay frozen in their own file).
**Holdout:** `…-round2-holdout-2026-09-21.json`, a fresh seed
(`1213743172`) and two targets on the largest instances of each regime.
**Diagnostic:** `…-round2-unguarded-diagnostic-2026-09-21.json`
(`--unguarded-targets`), the measurement of §10.2.
**Comparison:** `ic-boundary-round2-comparison-2026-09-21.json`, written
by `docs/ic/tools/boundary_round_compare.py`, which pairs every Round-2
row with the rung it was built on and is the `AGENTS.md` §8 evidence
(§10.7).

### 10.1 What the first round said to move, and what was moved

§4.1 named the phase that dominates each row: the pair table on the
prime `m = 3` row (`|F|²/2` additions, once), the *making* of the targets
on the prime `m = 2` row (two scalar multiplications per trial, about
`65` additions each at 24 bits), the `|F|` subtractions per target and
the table behind them on the binary row, and the table on every Koblitz
row until the `n = 41` cap.  §3.5 owed an accounting correction: the
counting ceiling divides by all of `#E` where a base confined to a
subgroup can only reach that subgroup.  Round 2 moves those and nothing
else — same curves, same subgroups, same seeds, same relation loop and
elimination — so every first-round row stays on the table as the
*before* mark of the rung built on it.

| rung (suffix) | what changed | what did not | class by the `AGENTS.md` §3 test |
|:--|:--|:--|:--|
| `_negfold` | the pair table is built once per pair up to negation: `{P, Q}` and `{−P, −Q}` sum to `±(P+Q)`, whose abscissa was already the key, so the first point above each abscissa is added to every signed point above a later one and to itself — `|F|²/4` additions instead of `|F|²/2`, the same key set, the same probe | trials, relations, yield, the elimination, every count outside the table | **engineering**: `S` falls by the table's share, the ratio to the counting boundary is flat |
| `_frobfold` (Koblitz) | the table is built once per pair up to `⟨σ, −1⟩`, keyed by the normal-basis canonical form of the sum's abscissa (`FrobeniusCanon`, the repository's rotation key): `σ^t(P) + σ^t(Q) = σ^t(P + Q)`, so all `2n` images of a pair share one entry — `|F|²/(4n)` entries and additions.  A probe canonicalises the target's abscissa once; the entry's shift and the target's give the rotation `t`; the summands `σ^t(P), σ^t(Q)` are read from a per-point orbit index; one addition confirms `±R`.  Canonicalisations are counted in their own unit and converted at a measured `ns_per_canon`, and every folded probe is cross-checked (`frobfold_mismatches`, zero on every row of every run) | as above.  The `2n` is the automorphism group the floor already credits to a generic algorithm, and the reference walk uses it too | **engineering** |
| `_walk` | targets come from a 16-jump r-adding walk with the coefficients tracked modulo `r`, one addition per target, instead of `R = [a]G + [b]Q` with two scalar multiplications.  Fresh jumps at every restart and one guard for the whole run, for the reason in §10.2 | the oracle, the yield per target against the ceiling, the elimination, the verification | **engineering**: this is the reference's own target generator |
| `mitm_m2_…` | two summands wherever the cofactor classes admit them and the exact floor fits the budget — the first round's census of 64 targets could not see a ceiling of `10⁻⁵`, which the exact formula computes instead of sampling | the base, the table, the loop | **engineering**: the ceiling moves by the counting formula of §1.3, the ratio to it does not |
| `_balanced` (Koblitz, `n ≥ 37`) | a base sized so the folded table balances a two-summand walk's trials, `|F| ≈ 1.2·(4#E)^{1/3}` signed points within `2^23` folded entries, where that is at least half again the first-round base | the pipeline | **engineering**: the base size was a parameter the full table's budget had set |
| exact ceiling (every row) | `min(1, N₀/r)` with `N₀` the `m`-multisets of base points whose cofactor classes `[r]P` cancel, reported next to the uniform `C(F+m−1, m)/#E`, with `yield_over_ceiling_exact` and `trials_floor_exact` | every count | **accounting**: the §3.5 correction, no `S` moves |

Nothing in that list does what a generic algorithm cannot: the two folds
use automorphisms the floor's `A` already grants, the walk is the
reference's own device, and the base size and summand count are
parameters of the same counting bound.  So the round can only be
engineering, and the table reads accordingly — `S` falls, the ratio to
the floor falls with it by the same factor, and the columns that would
mark an advance (`yield/ceiling` against the exact ceiling, and the
fitted exponents) do not move.

### 10.2 An artifact caught before it was reported: a repeated target is a collision, not a relation

**The mechanism.**  A relation loop that decomposes the **same group
element twice** hands the elimination two rows with the same factor-base
part and different `(a, b)`.  Those two rows pin the logarithm by
themselves: `a₁ + b₁d = a₂ + b₂d`.  That is a generic collision resolved
through the factor base, not a relation search, and it arrives after
about `√r` targets whatever the oracle costs — so a run that takes it is
measuring rho with extra steps and reporting the result as index
calculus.

**The rate is about one run in four, and it does not depend on the
instance.**  A two-summand search draws `T ≈ K/p` targets, of which
`T²/2r` pairs collide, and a fraction `p` of those are decomposable, so
the expected number of repeated rows is `K²/(2rp)`; with
`p ≈ C(F+1,2)/#E ≈ F²/2r` and `K = F/2` that is `1/4`, independently of
`F`, `K` and `r`.  Measured on the unguarded diagnostic
(`--unguarded-targets`, 198 rows):

| target source | `m` | runs pinned by a repeated target |
|:--|--:|--:|
| `[a]G + [b]Q` per trial | 2 | **11 of 45 = 24.4%** |
| `[a]G + [b]Q` per trial | 3 | 0 of 75 = 0% |
| r-adding walk, one jump table for every segment | 2 | **32 of 39 = 82.1%** |
| r-adding walk, one jump table for every segment | 3 | 2 of 39 = 5.1% |

The `24.4%` is the predicted `1/4`.  The three-summand rows are immune
at these sizes because they need far fewer targets than `√r`.  The walk
is worse than the random draw because segments sharing one step function
**merge** exactly as rho's walks do, so the second segment runs into the
first segment's path deliberately rather than by birthday.

**The fix, and what it cost.**  A repeat carries no relation the matrix
does not already hold — its factor-base part is one the elimination has
— so skipping it costs the relation search nothing.  Both target sources
are now guarded by one hash insert per target, counted as
`target_guard_probes`, priced as a lookup, and reported with
`repeated_targets_skipped`; the walk additionally draws sixteen fresh
jumps at every restart so no two segments share a step function.  The
loop counts `repeated_column_rows` and `pinned_by_repeated_row`, which
are zero on every row of all three guarded runs.  Two tests hold the
ends: `an_unguarded_draw_pins_the_logarithm_by_a_repeated_target` and
`a_shared_jump_table_lets_walk_segments_merge_and_pin_by_a_repeated_row`.
Against the unguarded diagnostic on the same rows, the guard costs
`1.10×` in `S` on the prime 24-bit two-summand walk, `2.72×` at
`n = 41` and `2.90×` on the binary `n = 27` row — that is the size of
what the artifact had been hiding.

**It is a correction to the first round too.**  Round 1 drew targets
without the guard, so 22 of its 63 rows moved when the guarded ladder
reran them, by `0.88×` to `1.59×` in `S` — the largest being
`bench-18bit`, where all three `m = 2` variants ran 674 trials and now
run 1,068 (`S` `81.1 → 119` for meet in the middle).  §2.2's `m = 2`
rows are to be read with that correction, and the comparison file
carries it row by row.

**What is *not* an artifact.**  A cycle among two-summand relations —
rows over *distinct* targets whose factor-base parts sum to zero — is
ordinary linear algebra and the classical way a two-summand index
calculus closes.  It is what the corrected rows do: at `n = 41` the
elimination finishes with 28 rows over 62 columns.  The line is whether
the group element repeats (a collision) or the column vectors combine (a
relation).

### 10.3 The table: every rung on the largest instance of its regime

| regime, instance | variant | m | \|F\| | K | trials | yield/ceiling (exact) | S | was | vs rho | vs floor | ok | class |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| **prime**, `generated-24bit-10935329`, `r = 2^23.4`, `#E = 1r`, `A = 2` | generic floor `√(π/2A)` | | | | | | 0.886 | | 0.23× | 1× | — | boundary |
| | Pollard rho, r-adding, counted (walk alone 3.54) | | | | | | 3.93 | | 1× | 4.43× | ✓ | reference |
| | Semaev `S₃` roots | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 5,129 |  | 1,305× | 5,787× | ✓ | baseline |
| | direct subtraction | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 2,037 |  | 518× | 2,298× | ✓ | accounting |
| | meet in the middle | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 268 |  | 68.3× | 303× | ✓ | engineering |
| | meet in the middle | 3 | 512 | 256 | 265 | 0.87 (0.87) | **56.8** |  | **14.5×** | 64.1× | ✓ | engineering |
| | + negation-folded table | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 248 | 268 | 63.2× | 280× | ✓ | engineering |
| | + walk targets | 2 | 512 | 256 | 7,974 | 0.93 (0.93) | **24.2** | 248 | **6.16×** | 27.3× | ✓ | engineering |
| | + negation-folded table | 3 | 512 | 256 | 265 | 0.87 (0.87) | 37.0 | 56.8 | 9.41× | 41.7× | ✓ | engineering |
| | + walk targets | 3 | 512 | 256 | 259 | 0.89 (0.89) | 31.7 | 37.0 | 8.06× | 35.7× | ✓ | engineering |
| **binary**, `random-binary-n27-b845462`, `r = 2^24.4`, `#E = 6r`, `A = 2` | generic floor `√(π/2A)` | | | | | | 0.886 | | 0.38× | 1× | — | boundary |
| | Pollard rho, r-adding, counted (walk alone 1.98) | | | | | | 2.31 | | 1× | 2.61× | ✓ | reference |
| | meet in the middle | 3 | 526 | 263 | 1,498 | 0.90 (0.90) | **206** |  | **89.1×** | 232× | ✓ | engineering |
| | `S₄` pairs-and-solve | 3 | 526 | 263 | 1,499 | 0.90 (0.90) | 108,654 |  | 46,965× | 122,603× | ✓ | relabelling |
| | + negation-folded table | 3 | 526 | 263 | 1,498 | 0.90 (0.90) | 191 | 206 | 82.7× | 216× | ✓ | engineering |
| | + walk targets | 3 | 526 | 263 | 1,482 | 0.90 (0.90) | 169 | 191 | 73.2× | 191× | ✓ | engineering |
| | + negation-folded table | 2 | 526 | 263 | 112,779 | 1.06 (1.05) | 1,645 | 206 | 711× | 1,856× | ✓ | engineering |
| | + walk targets | 2 | 526 | 263 | 75,384 | 1.14 (1.13) | **71.2** | 1,645 | **30.8×** | 80.4× | ✓ | engineering |
| **Koblitz**, `K_0 / GF(2^41)`, `r = 2^39.0`, `#E = 4r`, `A = 82` | generic floor `√(π/2A)` | | | | | | 0.138 | | 0.71× | 1× | — | boundary |
| | signed-Frobenius rho, counted (walk alone 0.19) | | | | | | 0.20 | | 1× | 1.41× | ✓ | reference |
| | meet in the middle, signed-orbit columns | 3 | 5,003 | 62 | 6,086 | 1.05 (1.04) | **58.8** |  | **300×** | 425× | ✓ | advance, count |
| | + negation-folded table | 3 | 5,003 | 62 | 6,086 | 1.05 (1.04) | 50.4 | 58.8 | 257× | 364× | ✓ | engineering |
| | + Frobenius-folded table | 3 | 5,003 | 62 | 6,086 | 1.05 (1.04) | 45.2 | 58.8 | 231× | 327× | ✓ | engineering |
| | + walk targets | 3 | 5,003 | 62 | 5,645 | 1.07 (1.07) | 41.1 | 45.2 | 210× | 297× | ✓ | engineering |
| | two summands, folded table, walk targets | 2 | 5,003 | 62 | 5,346,681 | 0.99 (0.96) | 8.15 | 58.8 | 41.6× | 58.9× | ✓ | engineering |
| | balanced base, three summands | 3 | 20,501 | 251 | 468 | 0.77 (0.77) | 12.3 | 41.1 | 62.7× | 88.7× | ✓ | engineering |
| | balanced base, two summands | 2 | 20,501 | 251 | 901,943 | 1.12 (1.12) | **5.12** | 8.15 | **26.2×** | 37.0× | ✓ | engineering |

The full Round-2 ladder, every instance and rung, with the phase split:

<details><summary>every Round-2 row on every instance</summary>

| regime | instance | log₂ r | variant | m | \|F\| | K | table | targets | trials | y/c | exact | S | was | vs rho | vs floor | FB | rel | LA | ok |
|:--|:--|--:|:--|--:|--:|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| prime | bench-10bit | 9.7 | mitm_m2_negfold | 2 | 32 | 16 | negation | random | 31 | 0.76 | 0.76 | 36.7 | 45.7 | 1.62 | 41.5 | 9.93 | 26.4 | 0.09 | ✓ |
| prime | bench-10bit | 9.7 | mitm_m2_negfold_walk | 2 | 32 | 16 | negation | walk | 28 | 0.66 | 0.66 | 25.9 | 36.7 | 1.14 | 29.3 | 9.93 | 15.6 | 0.07 | ✓ |
| prime | bench-10bit | 9.7 | mitm_m3_negfold | 3 | 32 | 16 | negation | random | 14 | 1.00 | 1.00 | 23.9 | 32.8 | 1.06 | 27.0 | 9.93 | 13.5 | 0.16 | ✓ |
| prime | bench-10bit | 9.7 | mitm_m3_negfold_walk | 3 | 32 | 16 | negation | walk | 13 | 1.00 | 1.00 | 26.3 | 23.9 | 1.16 | 29.7 | 9.93 | 15.9 | 0.14 | ✓ |
| prime | bench-12bit | 11.9 | mitm_m2_negfold | 2 | 32 | 16 | negation | random | 95 | 1.03 | 1.03 | 51.6 | 55.7 | 2.78 | 58.2 | 4.70 | 46.6 | 0.05 | ✓ |
| prime | bench-12bit | 11.9 | mitm_m2_negfold_walk | 2 | 32 | 16 | negation | walk | 93 | 1.11 | 1.11 | 20.8 | 51.6 | 1.12 | 23.4 | 4.70 | 15.8 | 0.06 | ✓ |
| prime | bench-12bit | 11.9 | mitm_m3_negfold | 3 | 32 | 16 | negation | random | 18 | 0.83 | 0.83 | 17.6 | 21.7 | 0.95 | 19.8 | 4.70 | 12.6 | 0.10 | ✓ |
| prime | bench-12bit | 11.9 | mitm_m3_negfold_walk | 3 | 32 | 16 | negation | walk | 17 | 0.84 | 0.84 | 17.2 | 17.6 | 0.92 | 19.4 | 4.70 | 12.1 | 0.10 | ✓ |
| prime | bench-14bit | 14.0 | mitm_m2_negfold | 2 | 64 | 32 | negation | random | 161 | 1.25 | 1.25 | 56.1 | 64.1 | 6.29 | 63.3 | 8.74 | 47.2 | 0.05 | ✓ |
| prime | bench-14bit | 14.0 | mitm_m2_negfold_walk | 2 | 64 | 32 | negation | walk | 154 | 1.20 | 1.20 | 18.7 | 56.1 | 2.09 | 21.1 | 8.74 | 9.73 | 0.05 | ✓ |
| prime | bench-14bit | 14.0 | mitm_m3_negfold | 3 | 64 | 32 | negation | random | 28 | 0.96 | 0.96 | 19.8 | 27.8 | 2.22 | 22.3 | 8.74 | 10.8 | 0.13 | ✓ |
| prime | bench-14bit | 14.0 | mitm_m3_negfold_walk | 3 | 64 | 32 | negation | walk | 29 | 0.94 | 0.94 | 17.2 | 19.8 | 1.93 | 19.4 | 8.74 | 8.15 | 0.16 | ✓ |
| prime | bench-16bit | 16.0 | mitm_m2_negfold | 2 | 128 | 64 | negation | random | 360 | 0.92 | 0.92 | 78.9 | 95.0 | 12.4 | 89.0 | 18.0 | 60.7 | 0.04 | ✓ |
| prime | bench-16bit | 16.0 | mitm_m2_negfold_walk | 2 | 128 | 64 | negation | walk | 314 | 1.13 | 1.13 | 22.5 | 78.9 | 3.53 | 25.4 | 18.0 | 4.32 | 0.04 | ✓ |
| prime | bench-16bit | 16.0 | mitm_m3_negfold | 3 | 128 | 64 | negation | random | 53 | 0.99 | 0.99 | 29.7 | 45.8 | 4.67 | 33.5 | 18.0 | 11.3 | 0.25 | ✓ |
| prime | bench-16bit | 16.0 | mitm_m3_negfold_walk | 3 | 128 | 64 | negation | walk | 53 | 0.99 | 0.99 | 23.9 | 29.7 | 3.75 | 26.9 | 18.0 | 5.54 | 0.22 | ✓ |
| prime | bench-18bit | 18.0 | mitm_m2_negfold | 2 | 128 | 64 | negation | random | 1,068 | 1.08 | 1.08 | 111 | 119 | 26.4 | 125 | 8.37 | 102 | 0.01 | ✓ |
| prime | bench-18bit | 18.0 | mitm_m2_negfold_walk | 2 | 128 | 64 | negation | walk | 1,581 | 0.89 | 0.89 | 16.7 | 111 | 3.97 | 18.8 | 8.37 | 8.24 | 0.02 | ✓ |
| prime | bench-18bit | 18.0 | mitm_m3_negfold | 3 | 128 | 64 | negation | random | 78 | 0.75 | 0.75 | 24.2 | 32.2 | 5.76 | 27.3 | 8.37 | 15.6 | 0.17 | ✓ |
| prime | bench-18bit | 18.0 | mitm_m3_negfold_walk | 3 | 128 | 64 | negation | walk | 80 | 0.72 | 0.72 | 18.8 | 24.2 | 4.48 | 21.3 | 8.37 | 10.3 | 0.12 | ✓ |
| prime | bench-20bit | 20.0 | mitm_m2_negfold | 2 | 256 | 128 | negation | random | 1,746 | 1.20 | 1.20 | 110 | 126 | 41.1 | 125 | 16.4 | 94.0 | 0.01 | ✓ |
| prime | bench-20bit | 20.0 | mitm_m2_negfold_walk | 2 | 256 | 128 | negation | walk | 2,258 | 0.95 | 0.95 | 22.4 | 110 | 8.34 | 25.3 | 16.4 | 5.99 | 0.01 | ✓ |
| prime | bench-20bit | 20.0 | mitm_m3_negfold | 3 | 256 | 128 | negation | random | 118 | 0.93 | 0.93 | 29.0 | 45.0 | 10.8 | 32.7 | 16.4 | 12.3 | 0.24 | ✓ |
| prime | bench-20bit | 20.0 | mitm_m3_negfold_walk | 3 | 256 | 128 | negation | walk | 120 | 0.94 | 0.94 | 23.4 | 29.0 | 8.72 | 26.4 | 16.4 | 6.71 | 0.32 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | mitm_m2_negfold | 2 | 512 | 256 | negation | random | 3,118 | 1.00 | 1.00 | 139 | 176 | 46.8 | 157 | 36.5 | 103 | 0.02 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | mitm_m2_negfold_walk | 2 | 512 | 256 | negation | walk | 3,242 | 1.06 | 1.06 | 39.8 | 139 | 13.4 | 44.9 | 36.5 | 3.23 | 0.02 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | mitm_m3_negfold | 3 | 512 | 256 | negation | random | 218 | 1.00 | 1.00 | 48.1 | 84.3 | 16.2 | 54.3 | 36.5 | 11.1 | 0.51 | ✓ |
| prime | generated-22bit-3290411 | 21.7 | mitm_m3_negfold_walk | 3 | 512 | 256 | negation | walk | 210 | 1.00 | 1.00 | 41.2 | 48.1 | 13.8 | 46.5 | 36.5 | 4.16 | 0.49 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | mitm_m2_negfold | 2 | 512 | 256 | negation | random | 11,612 | 0.99 | 0.99 | 248 | 268 | 63.2 | 280 | 20.4 | 228 | 0.01 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | mitm_m2_negfold_walk | 2 | 512 | 256 | negation | walk | 7,974 | 0.93 | 0.93 | 24.2 | 248 | 6.16 | 27.3 | 20.4 | 3.84 | 0.00 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | mitm_m3_negfold | 3 | 512 | 256 | negation | random | 265 | 0.87 | 0.87 | 37.0 | 56.8 | 9.41 | 41.7 | 20.4 | 16.1 | 0.51 | ✓ |
| prime | generated-24bit-10935329 | 23.4 | mitm_m3_negfold_walk | 3 | 512 | 256 | negation | walk | 259 | 0.89 | 0.89 | 31.7 | 37.0 | 8.06 | 35.7 | 20.4 | 10.8 | 0.53 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | mitm_m3_negfold | 3 | 30 | 15 | negation | random | 106 | 0.95 | 0.96 | 67.2 | 69.0 | 8.04 | 75.9 | 13.8 | 53.3 | 0.06 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | mitm_m3_negfold_walk | 3 | 30 | 15 | negation | walk | 140 | 0.71 | 0.71 | 52.4 | 67.2 | 6.27 | 59.1 | 13.8 | 38.5 | 0.04 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | mitm_m2_negfold | 2 | 30 | 15 | negation | random | 631 | 1.05 | 1.02 | 196 | 69.0 | 23.5 | 221 | 13.8 | 182 | 0.01 | ✓ |
| char2 | random-binary-n15-b524b | 14.0 | mitm_m2_negfold_walk | 2 | 30 | 15 | negation | walk | 528 | 1.25 | 1.20 | 63.9 | 196 | 7.65 | 72.1 | 13.8 | 50.0 | 0.01 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | mitm_m3_negfold | 3 | 64 | 32 | negation | random | 97 | 1.61 | 0.84 | 74.3 | 80.0 | 9.83 | 83.9 | 24.6 | 49.6 | 0.05 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | mitm_m3_negfold_walk | 3 | 64 | 32 | negation | walk | 107 | 1.72 | 0.90 | 59.0 | 74.3 | 7.80 | 66.5 | 24.6 | 34.2 | 0.07 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | mitm_m2_negfold | 2 | 64 | 32 | negation | random | 908 | 2.28 | 1.03 | 226 | 80.0 | 29.9 | 255 | 24.6 | 201 | 0.01 | ✓ |
| char2 | random-binary-n18-b6507 | 15.0 | mitm_m2_negfold_walk | 2 | 64 | 32 | negation | walk | 1,206 | 1.75 | 0.79 | 108 | 226 | 14.3 | 122 | 24.6 | 83.5 | 0.01 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | mitm_m3_negfold | 3 | 122 | 61 | negation | random | 407 | 0.96 | 0.96 | 77.3 | 81.0 | 24.3 | 87.3 | 11.6 | 65.7 | 0.04 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | mitm_m3_negfold_walk | 3 | 122 | 61 | negation | walk | 468 | 0.84 | 0.84 | 64.6 | 77.3 | 20.3 | 72.8 | 11.6 | 52.9 | 0.04 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | mitm_m2_negfold | 2 | 122 | 61 | negation | random | 10,149 | 0.99 | 0.98 | 557 | 81.0 | 175 | 629 | 11.6 | 546 | 0.00 | ✓ |
| char2 | random-binary-n21-b1b6f3b | 20.0 | mitm_m2_negfold_walk | 2 | 122 | 61 | negation | walk | 10,919 | 0.89 | 0.88 | 78.5 | 557 | 24.7 | 88.5 | 11.6 | 66.9 | 0.00 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | mitm_m3_negfold | 3 | 274 | 137 | negation | random | 365 | 1.74 | 0.87 | 88.9 | 102 | 35.3 | 100 | 23.7 | 65.0 | 0.15 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | mitm_m3_negfold_walk | 3 | 274 | 137 | negation | walk | 386 | 1.63 | 0.81 | 79.4 | 88.9 | 31.5 | 89.6 | 23.7 | 55.5 | 0.13 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | mitm_m2_negfold | 2 | 274 | 137 | negation | random | 16,514 | 2.05 | 1.02 | 686 | 102 | 273 | 774 | 23.7 | 662 | 0.00 | ✓ |
| char2 | random-binary-n24-b5fc9da | 21.0 | mitm_m2_negfold_walk | 2 | 274 | 137 | negation | walk | 14,315 | 1.87 | 0.92 | 62.1 | 686 | 24.7 | 70.1 | 23.7 | 38.4 | 0.00 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | mitm_m3_negfold | 3 | 526 | 263 | negation | random | 1,498 | 0.90 | 0.90 | 191 | 206 | 82.7 | 216 | 22.7 | 168 | 0.20 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | mitm_m3_negfold_walk | 3 | 526 | 263 | negation | walk | 1,482 | 0.90 | 0.90 | 169 | 191 | 73.2 | 191 | 22.7 | 146 | 0.19 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | mitm_m2_negfold | 2 | 526 | 263 | negation | random | 112,779 | 1.06 | 1.05 | 1,645 | 206 | 711 | 1,856 | 22.7 | 1,622 | 0.00 | ✓ |
| char2 | random-binary-n27-b845462 | 24.4 | mitm_m2_negfold_walk | 2 | 526 | 263 | negation | walk | 75,384 | 1.14 | 1.13 | 71.2 | 1,645 | 30.8 | 80.4 | 22.7 | 48.5 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | mitm_m3_signed_orbit_columns_negfold | 3 | 45 | 3 | negation | random | 5 | 0.82 | 0.82 | 23.6 | 39.0 | 1.36 | 88.4 | 17.7 | 5.59 | 0.02 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | mitm_m3_signed_orbit_columns_frobfold | 3 | 45 | 3 | frobenius | random | 5 | 0.82 | 0.82 | 11.3 | 39.0 | 0.65 | 42.2 | 5.24 | 5.69 | 0.02 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 45 | 3 | frobenius | walk | 4 | 1.00 | 1.00 | 19.8 | 11.3 | 1.14 | 74.1 | 5.24 | 14.2 | 0.03 | ✓ |
| koblitz | K_1 / GF(2^11) | 10.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 45 | 3 | frobenius | walk | 4 | 1.00 | 0.98 | 19.2 | 39.0 | 1.11 | 71.9 | 5.24 | 13.7 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | mitm_m3_signed_orbit_columns_negfold | 3 | 79 | 4 | negation | random | 4 | 1.00 | 1.00 | 40.0 | 74.0 | 2.98 | 163 | 36.8 | 2.83 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | mitm_m3_signed_orbit_columns_frobfold | 3 | 79 | 4 | frobenius | random | 4 | 1.00 | 1.00 | 10.6 | 74.0 | 0.79 | 43.3 | 7.49 | 2.84 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 79 | 4 | frobenius | walk | 4 | 1.00 | 1.00 | 18.8 | 10.6 | 1.40 | 76.5 | 7.49 | 11.0 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^13) | 11.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 79 | 4 | frobenius | walk | 6 | 1.65 | 0.84 | 18.5 | 74.0 | 1.38 | 75.4 | 7.49 | 10.8 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m3_signed_orbit_columns_negfold | 3 | 33 | 3 | negation | random | 8 | 2.70 | 0.89 | 27.2 | 36.5 | 1.42 | 119 | 11.3 | 15.6 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m3_signed_orbit_columns_frobfold | 3 | 33 | 3 | frobenius | random | 8 | 2.70 | 0.89 | 21.5 | 36.5 | 1.12 | 94.2 | 5.10 | 16.1 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 33 | 3 | frobenius | walk | 10 | 2.12 | 0.70 | 31.5 | 21.5 | 1.64 | 138 | 5.10 | 26.1 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 33 | 3 | frobenius | walk | 45 | 11.23 | 3.05 | 65.2 | 36.5 | 3.41 | 285 | 5.10 | 59.8 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m3_signed_orbit_columns_frobfold_walk_balanced | 3 | 91 | 4 | frobenius | walk | 4 | 1.00 | 1.00 | 31.2 | 31.5 | 1.63 | 136 | 14.3 | 16.5 | 0.03 | ✓ |
| koblitz | K_0 / GF(2^15) | 9.6 | mitm_m2_signed_orbit_columns_frobfold_walk_balanced | 2 | 91 | 4 | frobenius | walk | 6 | 3.76 | 0.61 | 29.8 | 65.2 | 1.55 | 130 | 14.3 | 15.1 | 0.01 | ✓ |
| koblitz | K_1 / GF(2^17) | 16.0 | mitm_m2_signed_orbit_columns_negfold | 2 | 239 | 8 | negation | random | 20 | 1.55 | 0.77 | 60.3 | 116 | 16.8 | 281 | 56.8 | 3.44 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^17) | 16.0 | mitm_m2_signed_orbit_columns_frobfold | 2 | 239 | 8 | frobenius | random | 22 | 1.47 | 0.73 | 9.96 | 116 | 2.77 | 46.4 | 6.16 | 3.72 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^17) | 16.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 239 | 8 | frobenius | walk | 15 | 2.24 | 1.12 | 9.18 | 9.96 | 2.56 | 42.7 | 6.16 | 2.94 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | mitm_m3_signed_orbit_columns_negfold | 3 | 305 | 9 | negation | random | 9 | 1.00 | 1.00 | 47.2 | 92.2 | 22.6 | 232 | 46.0 | 1.06 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | mitm_m3_signed_orbit_columns_frobfold | 3 | 305 | 9 | frobenius | random | 9 | 1.00 | 1.00 | 5.28 | 92.2 | 2.53 | 26.0 | 4.16 | 1.07 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 305 | 9 | frobenius | walk | 8 | 1.00 | 1.00 | 6.20 | 5.28 | 2.97 | 30.5 | 4.16 | 1.98 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^19) | 18.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 305 | 9 | frobenius | walk | 67 | 1.46 | 1.37 | 5.99 | 92.2 | 2.87 | 29.4 | 4.16 | 1.78 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | mitm_m3_signed_orbit_columns_negfold | 3 | 875 | 20 | negation | random | 18 | 1.00 | 1.00 | 94.7 | 188 | 106 | 512 | 93.9 | 0.79 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | mitm_m3_signed_orbit_columns_frobfold | 3 | 875 | 20 | frobenius | random | 18 | 1.00 | 1.00 | 6.35 | 188 | 7.11 | 34.4 | 5.53 | 0.80 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 875 | 20 | frobenius | walk | 19 | 1.00 | 1.00 | 6.33 | 6.35 | 7.09 | 34.3 | 5.53 | 0.78 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^23) | 22.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 875 | 20 | frobenius | walk | 281 | 1.00 | 0.98 | 6.21 | 188 | 6.96 | 33.6 | 5.53 | 0.67 | 0.00 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | mitm_m3_signed_orbit_columns_negfold | 3 | 3771 | 66 | negation | random | 59 | 0.99 | 0.99 | 17,311 | 34,555 | 4,121 | 105,189 | 17,274 | 36.8 | 0.22 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | mitm_m3_signed_orbit_columns_frobfold | 3 | 3771 | 66 | frobenius | random | 59 | 0.99 | 0.99 | 728 | 34,555 | 173 | 4,426 | 690 | 38.5 | 0.22 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 3771 | 66 | frobenius | walk | 53 | 1.00 | 1.00 | 713 | 728 | 170 | 4,335 | 690 | 23.6 | 0.15 | ✓ |
| koblitz | K_1 / GF(2^29) | 15.4 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 3771 | 66 | frobenius | walk | 350 | 1.38 | 0.30 | 703 | 34,555 | 167 | 4,271 | 690 | 13.2 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | mitm_m3_signed_orbit_columns_negfold | 3 | 2421 | 41 | negation | random | 49 | 0.74 | 0.74 | 1,268 | 2,489 | 1,119 | 7,968 | 1,224 | 44.7 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | mitm_m3_signed_orbit_columns_frobfold | 3 | 2421 | 41 | frobenius | random | 49 | 0.74 | 0.74 | 98.2 | 2,489 | 86.6 | 617 | 50.3 | 47.9 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 2421 | 41 | frobenius | walk | 51 | 0.70 | 0.70 | 100 | 98.2 | 88.4 | 630 | 50.3 | 49.9 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^31) | 20.5 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 2421 | 41 | frobenius | walk | 4,151 | 0.93 | 0.62 | 59.2 | 2,489 | 52.2 | 372 | 50.3 | 8.86 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m3_signed_orbit_columns_negfold | 3 | 4663 | 64 | negation | random | 493 | 0.98 | 0.98 | 499 | 857 | 1,195 | 3,426 | 358 | 141 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m3_signed_orbit_columns_frobfold | 3 | 4663 | 64 | frobenius | random | 493 | 0.98 | 0.98 | 161 | 857 | 387 | 1,108 | 11.3 | 150 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 4663 | 64 | frobenius | walk | 496 | 0.96 | 0.96 | 161 | 161 | 385 | 1,103 | 11.3 | 150 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 4663 | 64 | frobenius | walk | 91,721 | 1.07 | 0.92 | 19.4 | 857 | 46.4 | 133 | 11.3 | 8.13 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m3_signed_orbit_columns_frobfold_walk_balanced | 3 | 9177 | 125 | frobenius | walk | 175 | 0.70 | 0.70 | 94.1 | 161 | 225 | 646 | 41.8 | 52.3 | 0.01 | ✓ |
| koblitz | K_0 / GF(2^37) | 27.8 | mitm_m2_signed_orbit_columns_frobfold_walk_balanced | 2 | 9177 | 125 | frobenius | walk | 72,038 | 1.40 | 1.32 | 48.2 | 19.4 | 115 | 331 | 41.8 | 6.38 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m3_signed_orbit_columns_negfold | 3 | 4681 | 61 | negation | random | 921 | 2.05 | 1.02 | 1,170 | 1,831 | 2,569 | 8,242 | 662 | 508 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m3_signed_orbit_columns_frobfold | 3 | 4681 | 61 | frobenius | random | 921 | 2.05 | 1.02 | 561 | 1,831 | 1,233 | 3,955 | 19.8 | 541 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 4681 | 61 | frobenius | walk | 950 | 1.98 | 0.99 | 569 | 561 | 1,251 | 4,013 | 19.8 | 550 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 4681 | 61 | frobenius | walk | 224,721 | 4.84 | 0.71 | 110 | 1,831 | 241 | 774 | 19.8 | 90.1 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m3_signed_orbit_columns_frobfold_walk_balanced | 3 | 21529 | 277 | frobenius | walk | 275 | 0.94 | 0.94 | 522 | 569 | 1,146 | 3,678 | 391 | 131 | 0.11 | ✓ |
| koblitz | K_0 / GF(2^39) | 26.0 | mitm_m2_signed_orbit_columns_frobfold_walk_balanced | 2 | 21529 | 277 | frobenius | walk | 38,806 | 1.93 | 0.90 | 397 | 110 | 873 | 2,800 | 391 | 6.62 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m3_signed_orbit_columns_negfold | 3 | 5003 | 62 | negation | random | 6,086 | 1.05 | 1.04 | 50.4 | 58.8 | 257 | 364 | 8.45 | 41.9 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m3_signed_orbit_columns_frobfold | 3 | 5003 | 62 | frobenius | random | 6,086 | 1.05 | 1.04 | 45.2 | 58.8 | 231 | 327 | 0.24 | 45.0 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m3_signed_orbit_columns_frobfold_walk | 3 | 5003 | 62 | frobenius | walk | 5,645 | 1.07 | 1.07 | 41.1 | 45.2 | 210 | 297 | 0.24 | 40.9 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m2_signed_orbit_columns_frobfold_walk | 2 | 5003 | 62 | frobenius | walk | 5,346,681 | 0.99 | 0.96 | 8.15 | 58.8 | 41.6 | 58.9 | 0.24 | 7.91 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m3_signed_orbit_columns_frobfold_walk_balanced | 3 | 20501 | 251 | frobenius | walk | 468 | 0.77 | 0.77 | 12.3 | 41.1 | 62.7 | 88.7 | 3.80 | 8.48 | 0.00 | ✓ |
| koblitz | K_0 / GF(2^41) | 39.0 | mitm_m2_signed_orbit_columns_frobfold_walk_balanced | 2 | 20501 | 251 | frobenius | walk | 901,943 | 1.12 | 1.12 | 5.12 | 8.15 | 26.2 | 37.0 | 3.80 | 1.32 | 0.00 | ✓ |

</details>

### 10.4 Reading the rungs

**The folds do what their counts say, and no more.**  The negation fold
halves the table's additions and leaves every other count identical:
prime `m = 3` at 24 bits, `S` `56.8 → 37.0` with the same 265 trials, the
same 0.87 yield and the same elimination; Koblitz `m = 3` at `n = 41`,
`58.8 → 50.4` with the same 6,086 trials.  The Frobenius fold shrinks
the Koblitz table by the orbit length: `6,249,999 → 152,439` entries at
`n = 41`, which takes the factor-base phase from `16.9` to `0.24` in
`S`, and it pays `2.7` back in canonicalisations (one per probe, at
`50.8 ns` against an addition's `682.6 ns`, so `0.074` of an addition
each).  Net `50.4 → 45.2`.  On the large-cofactor rungs, where the base
was almost the whole cost, the same fold is worth far more: `n = 29`
goes `34,554 → 728`, `n = 31` `2,489 → 97.8`, `n = 37` `857 → 162`,
`n = 39` `1,831 → 562`.

**The walk is worth what the target generation cost.**  At `m = 3` a
target costs two scalar multiplications against `|F|` subtractions and
probes, so replacing the former by one addition moves little: prime
`37.0 → 31.7`, binary `191 → 169`, Koblitz `45.2 → 41.1`.  At `m = 2`
the target *was* the cost — one probe per target against `65` additions
to make it — and the walk is decisive: prime `248 → 24.2`, binary
`1,645 → 71.2`.  The relation phase of the prime `m = 2` row falls from
`228` to `3.84` in `S`, and the scalar multiplications from `23,224` to
`136`.

**Two summands beat three once targets are cheap.**  The exact ceiling
of §10.5 makes `m = 2` admissible where the first round's 64-target
census saw nothing: at `n = 41` the two-summand ceiling is `5.9 × 10⁻⁶`,
so the row needs `5.3` million walked targets — but each is one
addition, one canonicalisation and one probe, and the whole row costs
`S = 8.15` against the three-summand `41.1`.  On the balanced base it is
`5.12`, the best Koblitz row on the board at `26.2×` the
signed-Frobenius reference.  Note what pays for it: the elimination
finishes with 28 rows over 62 columns on the narrow base, through cycles
in the factor-base graph rather than through full rank — the classical
way a two-summand index calculus closes, and (per §10.2) verified not to
be a collision.

**The balanced base trades table for trials, and wins at `n ≥ 37`.**
Sizing the base to `1.2·(4#E)^{1/3}` takes `|F|` from `5,003` to
`20,501` at `n = 41` and `K` from `62` to `251`.  The folded table grows
to `2.6` million entries (`S` `0.24 → 3.80`) and the two-summand trials
fall `5.3` million `→ 902` thousand, for a net `8.15 → 5.12`.  At
`m = 3` the same base is `41.1 → 12.3`.

**Per instance, best row to best row**, the round is worth `1.3×` to
`2.9×` on the prime and binary ladders and `3×` to `50×` on the Koblitz
one, where the folds bite hardest: the ratios of the best first-round row
to the best Round-2 row run `1.26–2.35` (prime), `1.25–2.89` (binary) and
`3.4–49` (Koblitz).  Nothing crosses the reference.  The closest rows are
the prime `m = 2` walk at `6.2×` rho and the Koblitz balanced `m = 2` at
`26.2×`; at 12 bits the prime `m = 3` walk sits at `0.9×`, which is
rho's unamortised setup at small `r` and not a crossover (§2.3).

### 10.5 The exact ceiling: the † rows corrected

§3.5 marked six rows with a `yield/ceiling` above `1.4` and argued, by
hand, that the ceiling was loose by the index of the subgroup the base's
sums land in.  The exact ceiling computes that index instead of arguing
it: it counts the `m`-multisets of base points whose cofactor classes
cancel and divides by `r`, so a base confined to a proper subgroup is
measured against what it can actually reach.  Every row now carries
both.

| instance | rows | cofactor | classes | `yield/ceiling` | exact | index |
|:--|:--|--:|--:|--:|--:|--:|
| `random-binary-n18-b6507` | `m = 3` and `m = 2`, all rungs | 8 | 4 | 1.48–2.28 | 0.77–1.03 | 1.92–2.21 |
| `random-binary-n24-b5fc9da` | all rungs | 8 | 4 | 1.63–2.05 | 0.81–1.02 | 2.00–2.02 |
| `K_0 / GF(2^15)` | all rungs | 44 | 13 | 1.37–2.70 | 0.45–0.89 | 3.04 |
| `K_1 / GF(2^17)` | all `m = 2` rungs | 2 | 1 | 1.47–2.24 | 0.73–1.12 | 2.00 |
| `K_1 / GF(2^29)` | `m = 2` walk | 12,646 | 3,133 | 1.38 | 0.30 | 4.59 |

The hand-argued factor 2 of §3.5 is confirmed where the base is in `2E`
(`1.92–2.21` measured, against the exact 2 the trace argument gives, the
remainder being the multiset boundary terms), and the cases §3.5 could
only gesture at now have numbers: the index is `3.04` at `n = 15`, where
a cofactor of 44 leaves 13 classes, and `4.59` at `n = 29`, where 12,646
leaves 3,133.  No count moved and no `S` moved; this is the accounting
correction §3.5 asked for, and the `†` rows of §2.2 are to be read at
the exact column from here on.

One row still exceeds its exact ceiling: `K_0 / GF(2^15)`'s two-summand
folded walk, at `3.05`.  That instance has `r = 2^9.6` with a cofactor of
44 and a base of 33 signed points, so the run finds two or three
relations in a handful of trials and the ratio is the variance of a
three-run mean at single digits, not a structural yield.  It is flagged
here rather than smoothed away.

### 10.6 The exponents, refitted

Fitted total exponents of every Round-2 variant next to the rung it was
built on, from the run's `fits` array.  The folds move nothing — they
divide a phase by a constant — and the walk moves the `m = 2` rows,
because it replaces a per-target cost that grew with the scalar size by
one that does not.

| regime | variant | α (r) | R² | sizes | α of the rung it was built on |
|:--|:--|--:|--:|--:|--:|
| prime | rho reference | 0.274 | 0.903 | 8 | — |
| prime | semaev_s3_roots_m2 | 0.933 | 0.965 | 8 | — |
| prime | direct_subtraction_m2 | 0.871 | 0.997 | 8 | — |
| prime | mitm_m2 | 0.678 | 0.998 | 8 | — |
| prime | mitm_m3 | 0.601 | 0.982 | 8 | — |
| prime | mitm_m2_negfold | 0.681 | 0.997 | 8 | 0.678 |
| prime | mitm_m2_negfold_walk | 0.525 | 0.980 | 8 | 0.681 |
| prime | mitm_m3_negfold | 0.578 | 0.989 | 8 | 0.601 |
| prime | mitm_m3_negfold_walk | 0.554 | 0.982 | 8 | 0.578 |
| char2 | rho reference | 0.302 | 0.973 | 5 | — |
| char2 | mitm_m3 | 0.624 | 0.986 | 5 | — |
| char2 | semaev_s4_pairs_and_solve_m3 | 1.055 | 0.988 | 5 | — |
| char2 | mitm_m3_negfold | 0.618 | 0.985 | 5 | 0.624 |
| char2 | mitm_m3_negfold_walk | 0.637 | 0.987 | 5 | 0.618 |
| char2 | mitm_m2_negfold | 0.788 | 0.999 | 5 | 0.624 |
| char2 | mitm_m2_negfold_walk | 0.476 | 0.979 | 5 | 0.788 |
| koblitz | rho reference | 0.254 | 0.917 | 11 | — |
| koblitz | mitm_m3_signed_orbit_columns | 0.544 | 0.709 | 10 | — |
| koblitz | mitm_m3_signed_orbit_columns_negfold | 0.556 | 0.734 | 10 | 0.544 |
| koblitz | mitm_m3_signed_orbit_columns_frobfold | 0.588 | 0.831 | 10 | 0.544 |
| koblitz | mitm_m3_signed_orbit_columns_frobfold_walk | 0.558 | 0.827 | 10 | 0.588 |
| koblitz | mitm_m2_signed_orbit_columns_frobfold_walk | 0.446 | 0.789 | 11 | 0.544 |
| koblitz | mitm_m3_signed_orbit_columns_frobfold_walk_balanced | 0.471 | 0.861 | 4 | 0.558 |
| koblitz | mitm_m2_signed_orbit_columns_frobfold_walk_balanced | 0.430 | 0.821 | 4 | 0.446 |

### 10.7 Against the targets of §1.6, and the `AGENTS.md` §8 gate

**Neither condition of §1.6 is met, and the round does not claim them.**
No variant is below the reference on any instance with `r ≥ 2^20`: the
closest are the prime `m = 2` walk at `6.2×` rho (`2^23.4`) and the
Koblitz balanced `m = 2` at `26.2×` (`2^39.0`).  No fitted total exponent
is below the reference's: rho fits at `0.25–0.30` here and the best
Round-2 rows at `0.43–0.53` (§10.6).  No `yield/ceiling` exceeds `1.5`
against the exact ceiling on a base outside a proper subgroup (§10.5).
Every logarithm on every row of all three runs was recovered and
verified, every rho run recovered and verified its own, and no row was
pinned by a repeated target.

**The §8 comparison.**  103 candidate rows paired with the rung each was
built on, plus 34 holdout rows on a fresh seed, in
`ic-boundary-round2-comparison-2026-09-21.json`.  The speedups
(`baseline_total_operations / candidate_total_operations`, same curves,
subgroups, factor bases, targets and seeds, every phase charged) run
from `0.13` to `49.2` with a median of `1.41`; the holdout runs `0.07`
to `67.1`, median `1.31`.  Fourteen of the 103 are below one and are
reported as such: they are the rungs where a lever costs more than it
saves — the two-summand rows before the walk is applied to them
(`mitm_m2_negfold` against a three-summand baseline), and the balanced
base at the sizes where its larger table is not yet repaid.  A rung that
loses is a measurement, not a failure, and it stays on the table.

The comparison file also records what the target guard did to the
*first* round's rows: 41 of 63 reproduce the frozen Round-1 counts
exactly, and 22 moved, by `0.88×` to `1.59×` in `S`.  Those 22 are the
rows §10.2 corrects, and the correction is to Round 1's numbers, not to
this round's.

The frozen WDSat regression suite of `AGENTS.md` §8
(`research/index_calculus_baseline_20260914/regression/`) does not apply
to this round and was not run, for the reason its own README gives: it
measures one SAT-solver stage on sixty fixed Weil-descended inputs and
"is not an adapter" for any other pipeline.  None of the six levers
touches a solver stage — the pair table, the target generator, the
counting bound and the base size are the phases *around* the oracle, and
the oracle that finishes a logarithm here is a hash probe.  The matched
suite §8 asks for instead is the comparison file: every candidate row
paired with the rung it was built on, on identical curves, subgroups,
factor bases, targets and seeds, with every phase charged exclusively,
`speedup = baseline_total_operations / candidate_total_operations` per
repeat and as the mean, a fresh-seed holdout on the largest instances,
and a check that the first-round rows of the new run reproduce the
frozen counts of §2 exactly.  No runtime claim is made, so no paired
wall-clock interval is owed; the wall columns stay what they were in §2,
a practicality note.

### 10.8 What does not count, and what this is not

- **Not a crossover.**  The best row of the round is `6.2×` the counted
  reference on the prime ladder and `26.2×` the signed-Frobenius
  reference on the Koblitz one.  The `0.9×` at 12 bits is rho's setup
  being unamortised below `2^20` (§2.3), which the `vs floor` column
  does not flatter.
- **Not an advance.**  Every rung is engineering by the §3 test.  The
  two folds use automorphisms the floor's `A` already grants a generic
  algorithm, the walk is the reference's own target generator, and the
  summand count and base size are parameters of the same counting bound.
  The exact ceiling is accounting.  `yield/ceiling` against the exact
  ceiling stays at or below one on every base that is not confined to a
  proper subgroup.
- **Not an exponent result.**  The two-summand walk rows fit at `0.43`
  to `0.53`, which brackets rho's limiting one half, but over four to
  eight rungs whose base size steps by `2^{1/3}` and with the factor
  base and relation phases pulling opposite ways (§10.6).  §1.6 asks for
  an exponent below *the reference's*; the reference fits at `0.25–0.30`
  on these ladders.  The balanced rows in particular are four rungs at
  cofactors 44, 596, 8,012 and 4, reported with their `R²` and not
  offered as a law.
- **Not a claim about a deployed curve**, and not a wall-clock
  benchmark.  §8 stands unchanged.
- **The collision correction cuts both ways.**  §10.2 removes a
  mechanism that had been making some rows *look* cheaper, including 22
  of the first round's.  The Round-2 gains in §10.4 are measured after
  that removal, on rows where it applies equally to the before mark and
  the after mark.

### 10.9 Reproducing

```bash
cargo build --release --bin ic
./target/release/ic boundary --out round2.json
./target/release/ic boundary --seed 1213743172 --repeats 2 --prime-bits 22,24 \
    --char2-degrees 24,27 --koblitz-degrees 37,39,41 --s4-max-degree 24 --out holdout.json
./target/release/ic boundary --unguarded-targets --prime-bits 20,22,24 \
    --char2-degrees 21,24,27 --s4-max-degree 20 --koblitz-degrees 23,37,39,41 --out unguarded.json
python3 docs/ic/tools/boundary_round_compare.py round2.json \
    --baseline docs/ic/runs/ic-boundary-ledger-2026-09-21.json --holdout holdout.json --out comparison.json
python3 docs/ic/tools/boundary_ledger_tables.py round2.json
python3 docs/ic/tools/boundary_scoreboard_rows.py round2.json
cargo test --release --lib cryptanalysis::ic_
```

## 11. Round 3: the restart was the cost, and the family has a ceiling

**Frozen runs.**  `docs/ic/runs/ic-boundary-ledger-round3-2026-09-21.json`
(the ladder), `…-round3-holdout-2026-09-21.json` (a fresh seed), and a
**matched headline pair** at eight repeats:
`…-round2-headline-2026-09-21.json`, run on a binary built from a git
worktree at the Round-2 commit `acd192dc`, against
`…-round3-headline-2026-09-21.json` on this round's.  The intermediate
that the round discarded is frozen too — `…-round3-rotation-*.json`,
§11.1 — because it is the measurement that justified discarding it.
**Comparison:** `ic-boundary-round3-comparison-2026-09-21.json`.

### 11.1 What Round 2 left on the table

§10.4 reported the walk as an unqualified win: one addition per target
instead of two scalar multiplications.  That was true of the *steps* and
false of the *restarts*.  The guard of §10.2 forces a restart whenever a
segment reaches a group element the run has already decomposed, and
Round 2 answered each restart by drawing sixteen fresh jumps — sixteen
`[a]G + [b]Q`, two scalar multiplications each, about `1,170` additions
on a 24-bit curve — so that no two segments could share a step function
and merge.  On a rung that restarts often that is the dominant cost of
the relation phase: at `n = 27` the binary two-summand row restarted 132
times, drew 2,117 jumps, and spent about two thirds of its whole
relation phase on them.

Nothing about the guard requires a fresh jump *table*.  Two segments
merge when they share a step function, and a step function is a pair
(jump set, index map); Round 3 keeps one jump set for the run and gives
each segment its own index map.  Sixteen offsets `[c]G + [d]Q` are
pooled with the jumps at setup, and a restart adds one pooled offset to
the current point — one addition, the coefficients add modulo `r` — and
changes the map.  The restart therefore costs one group operation
instead of thirty-two scalar multiplications, and `walk_jumps` is 32 for
a whole run however often it restarts, which is what the test
`a_walk_restart_costs_one_addition_and_still_does_not_merge` asserts.

**A hypothesis about the index map, tested and not supported.**  The map
was at first a *rotation*: segment `j` selected `jump[(h(P) + j) mod 16]`.
That gives consecutive segments different step functions, but only
sixteen distinct ones in total, so after enough restarts two segments
share one and can merge as §10.2's shared jump table let them.  Some
Koblitz rows did move the wrong way under it — at `n = 37` the ladder's
guard-forced restarts went from 18 to 129 and the row got `1.35×` worse,
and on the holdout its walk steps per relation went `14,359 → 21,517` —
so the round replaced the rotation with a **permutation**: Fisher–Yates
over the sixteen indices at each restart, `16!` index maps instead of
`16`, two colliding segments diverging again within a step or two
because two random permutations agree on about one index in sixteen.
Shuffling sixteen entries costs no group operation, so the restart stays
one addition.

**Then the hypothesis was checked against the counter that measures it
directly, and it failed.**  A merge shows up as a run of targets the
guard has already seen, and `repeated_targets_skipped` counts exactly
those.  Per thousand trials, on the holdout's two-summand walk rows:

| instance | Round 2, fresh jumps | rotation | permutation |
|:--|--:|--:|--:|
| `generated-22bit` | 0.47 | 0.45 | 0.51 |
| `generated-24bit` | 0.76 | 0.45 | 0.46 |
| `random-binary-n24` | 2.92 | 4.28 | 3.61 |
| `random-binary-n27` | 1.21 | 0.76 | 0.66 |
| `K_0 / GF(2^37)` | 0.19 | 0.21 | 0.55 |
| `K_0 / GF(2^39)` | 0.38 | 0.99 | 1.31 |
| `K_0 / GF(2^41)` | 0.00 | 0.01 | 0.01 |

There is no systematic difference between sixteen step functions per
run, `16!` of them, and a fresh sixteen at every restart.  The rows that
moved the wrong way under the rotation moved within their own spread
(§11.4), not because their segments merged.

So the permutation **is not a measured improvement and is not reported
as one.**  It is retained because it removes a mechanism the ledger has
already been bitten by once, at a cost of zero group operations, and
because it is the version that shipped; the rotation runs stay frozen
beside it as what the comparison actually showed.  `repeated_column_rows`
and `pinned_by_repeated_row` are zero on every row of every run under
both.

The round's second lever gives the prime and binary regimes the
`_balanced` row the Koblitz regime got in §10.1, sized by the law of
§11.2 rather than by the `⌈bits/3⌉` rule.  The third is the law itself.

### 11.2 The family shape law, derived before measuring

Every two-summand row on this ledger pays exactly two things, and both
are functions of the one parameter it has, the base size `F` in signed
points.

- **The pair table**, one entry per pair of base points up to the fold:
  `F²/4` additions up to negation, `F²/(4n)` when the Frobenius folds it
  too.  Write `t` for that **table fold**, `1` or `n`.
- **The relations.**  A two-summand row is an *edge* on the `K = F/2k`
  columns, with `k` the **column fold** — `1` when a column is an
  abscissa, `n` when it is a signed Frobenius orbit.  A random graph
  first carries a cycle at about half as many edges as vertices, and it
  is that cycle which closes the elimination, so the row needs about
  `K/2` relations.  Each costs `1/p` targets with
  `p = C(F+1,2)/#E ≈ F²/2#E` the counting ceiling of §1.3, at `c` group
  operations per target — `c = 1` for a walk step.

So the family's whole cost:

```text
    ops(F) = F²/(4t) + c·#E/(2kF),
    least at  F* = (c·#E·t/k)^{1/3},  where  ops = 0.75·(c·#E·t/k)^{2/3} / t
```

and, in the unit of §1.1 with `c = 1`,

```text
    S_family = 0.75·(#E·t/k)^{2/3} / (t·√r).
```

**The two folds are independent**, which is the correction this round
had to make to its own first formula: a Koblitz row can fold its columns
without folding its table, and the ladder runs rows of exactly that
shape, so a single fold in both terms is the wrong model for them.  On
the prime and binary ladders both folds are one and nothing changes.

**The law is derived for two summands**, and the `vs family` column on
an `m = 3` row therefore reads "how this row compares with the best
*two-summand* member of the family on the same instance" — a common mark
across the ledger, not that row's own optimum.

**And it is a model, not a bound.**  The `K/2` is the first-cycle
threshold for a simple graph; the factor-base graph is a *multi*graph,
where two relations over the same pair of columns are a cycle of length
two and appear after about `√K` relations.  The two thresholds differ by
`√K/2`, so the relation term is an over-estimate somewhere in that
range, and a row can finish under `S_family` without doing anything a
counting argument forbids.  The Koblitz `n = 37` two-summand walk does,
at `0.79`, with its measured yield *below* the exact ceiling (`0.84`) —
which rules out the alternative explanation, that it beat the counting
bound.

### 11.3 The family cannot cross rho, and the ladder shows it converging

The unit was chosen so that rho is a constant: a counted r-adding walk
costs `√(πr/2A)` operations, so `S_rho = √(π/2A)` — `1.25` at `A = 1`,
`0.886` at `A = 2` — at every size.  §11.2's law says the best member of
the pair-table two-summand family costs `0.75·r^{1/6}` on a prime-order
curve, which is **unbounded**.  A constant and an unbounded increasing
function cross at most once, and these cross at `r ≈ 22`.

So no choice of base size, table fold or target generator inside this
family gives a crossover at any size that matters, and the gap widens as
`r^{1/6}`.  Against rho's asymptotic `1.25`:

| `r` | `S_family` | vs rho |
|:--|--:|--:|
| `2^24` | 12.0 | 9.6× |
| `2^48` | 192 | 153× |
| `2^64` | 1,219 | 973× |
| `2^128` | 1.98×10⁶ | 1.58×10⁶× |
| `2^160` | 7.99×10⁷ | 6.4×10⁷× |

Those are **extrapolations from the law**, not measurements, and they
use rho's asymptote rather than its measured `S`, which is the less
flattering choice: at `2^{24}` rho measures `3.93` because its setup is
not amortised, so the *measured* ratio there is `3.6×` and not `9.6×`.

This is a stronger statement than the earlier rounds could make.  §10.8
said "not a crossover" about the rows it had measured.  This says the
family has none to find.

**The ladder is the evidence that the law describes these rows**, and it
gives it by converging.  The balanced prime row takes the base size the
law prescribes at every rung:

| instance | `log₂ r` | `\|F\|` | `F*` | `S` | `S_family` | `S / S_family` |
|:--|--:|--:|--:|--:|--:|--:|
| bench-10bit | 9.7 | 10 | 9 | 32.90 | 2.30 | 14.33 |
| bench-12bit | 11.9 | 16 | 16 | 21.51 | 2.97 | 7.23 |
| bench-14bit | 14.0 | 26 | 25 | 14.28 | 3.77 | 3.78 |
| bench-16bit | 16.0 | 40 | 40 | 12.36 | 4.76 | 2.60 |
| bench-18bit | 18.0 | 64 | 64 | 10.72 | 6.00 | 1.79 |
| bench-20bit | 20.0 | 102 | 102 | 11.97 | 7.56 | 1.58 |
| generated-22bit | 21.7 | 148 | 149 | 10.47 | 9.15 | 1.15 |
| generated-24bit | 23.4 | 222 | 222 | 14.03 | 11.17 | 1.26 |

The same check on the Koblitz ladder, where the Frobenius fold divides
the family's cost by `n`: at `n = 41` the law puts the family's floor at
`S_family = 4.17`, which is `30.1×` the asymptotic signed-Frobenius rho
`√(π/4n) = 0.138`, and the best measured row sits at `31.5×`.  The
family is within half a tenth of its own optimum, and its own optimum is
thirty times the reference.

**That convergence is why the fitted exponent must not be
extrapolated.**  The balanced prime row fits `ops ∝ r^{0.407}`
(`R²` 0.97, 8 sizes) while the law predicts `r^{2/3}` — `S = 0.75 r^{1/6}`
means `ops = 0.75 r^{2/3}`.  The fit is flatter because the row starts
fourteen times above the law and is still falling toward it.  A reader
who extrapolated `0.407` would predict a crossover with rho that the law
says does not exist; that is the §5 mistake in a new costume, a number
measured inside a transient and read as an asymptote, and the fit is
reported here only beside the prediction it disagrees with.

### 11.4 The resolution of the measurement, measured

A row of one of these runs is a mean over repeats that differ only in
the target drawn.  A two-summand search finishes when the factor-base
graph first carries a cycle, and when that happens belongs to the draw,
so the same row, on the same instance, on the same binary, costs
several times more on one repeat than on another.  That spread is the
resolution of every comparison in this note and until this round it was
nowhere in it.  `docs/ic/tools/boundary_repeat_spread.py` reports it as
`max/min` of the per-repeat operation count:

| run | rows | `m` | median `max/min` | 75th | 90th | max |
|:--|--:|--:|--:|--:|--:|--:|
| Round-2 ladder, 3 repeats | 69 | 2 | 1.37 | 2.24 | 3.07 | 7.52 |
| Round-2 ladder, 3 repeats | 97 | 3 | 1.08 | 1.23 | 1.72 | 25.10 |
| headline pair, 8 repeats | 24 | 2 | 1.76 | 2.36 | 4.35 | 5.17 |
| headline pair, 8 repeats | 32 | 3 | 1.28 | 1.43 | 1.56 | 1.69 |

Two consequences, and both change how this round reports itself.

**Three repeats cannot resolve this round's gains.**  Round 2's levers
moved rows by factors of two to fifty and three repeats were enough.
Round 3's move them by tens of percent.  So the headline is run at
**eight repeats on two binaries** — the previous round's, built from a
worktree at its commit, and this one's — on the same instances with the
same seed, run concurrently so they saw the same machine.

**A mean of ratios is the wrong average.**  The paired ratios on the
Koblitz `n = 41` two-summand row were `0.425, 0.479, 0.662, 0.790,
1.110, 1.473, 2.208, 4.277`; their arithmetic mean is `1.43`, which
reads a spread as a gain.  `boundary_round_compare.py` now reports the
**geometric** mean with a two-sided 95% `t`-interval on the logs, and a
column that says whether the interval clears one.  Every speedup in
§11.6 is quoted that way, and the ones whose interval straddles one are
reported as not moved — which is a result, not a gap in the evidence.

## Appendix A. The conversion factors, as measured

Nanoseconds per native unit on the run's host, per instance, from the
report's `calibration` records.  A count divided by the host's
`ns_per_add` for the same curve is what every non-addition unit was
converted at; the square-root column is the one that moves the prime
`S₃` row (§3.2).

| regime | instance | add | double | sqrt | AS solve | S₄ pair | lookup | row op | word xor | legendre | inversion | frobenius |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| prime | bench-10bit | 140 | 165 | 63.5 | — | — | 2.63 | 3.57 | 0.50 | 49.5 | 95.5 | — |
| prime | bench-12bit | 153 | 177 | 77.0 | — | — | 2.64 | 3.57 | 0.50 | 66.5 | 109 | — |
| prime | bench-14bit | 166 | 190 | 123 | — | — | 2.62 | 3.59 | 0.50 | 84.0 | 123 | — |
| prime | bench-16bit | 177 | 203 | 993 | — | — | 2.71 | 3.43 | 0.49 | 87.9 | 132 | — |
| prime | bench-18bit | 195 | 218 | 155 | — | — | 2.94 | 3.51 | 0.49 | 107 | 144 | — |
| prime | bench-20bit | 203 | 227 | 171 | — | — | 2.64 | 3.51 | 0.49 | 118 | 159 | — |
| prime | generated-22bit-3290411 | 215 | 240 | 168 | — | — | 3.10 | 3.64 | 0.50 | 129 | 169 | — |
| prime | generated-24bit-10935329 | 224 | 248 | 1,060 | — | — | 3.00 | 3.51 | 0.49 | 131 | 180 | — |
| char2 | random-binary-n15-b524b | 238 | 238 | — | 13.1 | 676 | 2.70 | 3.55 | 0.49 | — | — | 14.2 |
| char2 | random-binary-n18-b6507 | 301 | 299 | — | 9.92 | 848 | 2.65 | 3.55 | 0.50 | — | — | 15.3 |
| char2 | random-binary-n21-b1b6f3b | 344 | 345 | — | 22.5 | 960 | 2.82 | 3.57 | 0.50 | — | — | 15.0 |
| char2 | random-binary-n24-b5fc9da | 387 | 387 | — | 13.0 | 1,050 | 2.59 | 3.51 | 0.49 | — | — | 15.0 |
| char2 | random-binary-n27-b845462 | 459 | 460 | — | 32.8 | 1,292 | 3.19 | 3.59 | 0.50 | — | — | 16.5 |
| koblitz | K_1 / GF(2^11) | 183 | 182 | — | 9.11 | — | 2.74 | 3.54 | 0.50 | — | — | 14.3 |
| koblitz | K_0 / GF(2^13) | 210 | 211 | — | 10.7 | — | 2.56 | 3.54 | 0.49 | — | — | 14.2 |
| koblitz | K_0 / GF(2^15) | 239 | 239 | — | 13.7 | 642 | 2.61 | 3.57 | 0.49 | — | — | 14.4 |
| koblitz | K_1 / GF(2^17) | 266 | 266 | — | 15.4 | — | 2.69 | 3.55 | 0.49 | — | — | 14.2 |
| koblitz | K_1 / GF(2^19) | 313 | 313 | — | 17.9 | — | 2.61 | 3.53 | 0.49 | — | — | 15.0 |
| koblitz | K_1 / GF(2^23) | 373 | 372 | — | 25.4 | — | 3.32 | 3.52 | 0.49 | — | — | 15.1 |
| koblitz | K_1 / GF(2^29) | 487 | 486 | — | 38.1 | — | 5.02 | 3.52 | 0.50 | — | — | 17.0 |
| koblitz | K_0 / GF(2^31) | 519 | 521 | — | 42.4 | 1,495 | 3.85 | 3.56 | 0.49 | — | — | 16.6 |
| koblitz | K_0 / GF(2^37) | 641 | 642 | — | 67.8 | — | 3.91 | 3.52 | 0.49 | — | — | 17.3 |
| koblitz | K_0 / GF(2^39) | 685 | 682 | — | 69.5 | — | 3.86 | 3.59 | 0.50 | — | — | 17.5 |
| koblitz | K_0 / GF(2^41) | 721 | 717 | — | 72.8 | — | 4.19 | 3.61 | 0.50 | — | — | 17.5 |

The word-XOR factor (`0.49–0.50 ns`) against an addition (`183–721 ns`)
is what turns matrix-F4's `2.3 × 10⁸` XORs into `3.1 × 10⁵` additions
in §5; the factors for the binary regimes rise with `n` because the
field arithmetic is a portable word-by-word implementation, and a
carry-less-multiply build would move every binary row down by the same
factor without changing any count or any ratio.
