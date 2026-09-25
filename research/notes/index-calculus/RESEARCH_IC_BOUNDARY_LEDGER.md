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

**Three summands obey the same law, with a worse constant.**  The
`vs family` column on an `m = 3` row reads "how this row compares with
the best *two-summand* member on the same instance", which is a common
mark and not that row's own optimum.  But the `m = 3` family's own
optimum is worth deriving, because it says whether the summand count is
a way out.  It is not.  An `m = 3` row builds the same pair table,
`F²/(4t)`; a target costs `F` probes rather than one, since the oracle
subtracts each base point and looks the remainder up; and a target
decomposes with `p₃ = C(F+2,3)/#E ≈ F³/6#E`.  So its relation term is

```text
    F · ρ·K · 6#E/F³  =  3ρ·#E/(kF)
```

with `ρ` the relations needed per column — `1/2` for the first cycle of
a graph, about `0.82` for the 2-core of a 3-uniform hypergraph.  That is
the *same* `1/F` shape as the two-summand term with a coefficient three
to five times larger, so the `m = 3` optimum sits at the same
`F ∝ (#E·t/k)^{1/3}` and costs roughly two to three times more.  The
exponent does not move.

**So nothing inside this family changes the exponent** — not the base
size, not the two folds, not the summand count.  They move the constant
in front of `r^{1/6}` and nothing else, which is what §11.3 turns into a
statement about crossovers.

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
They are also generous to index calculus in the other direction, for the
reason §2.3 gives: the prime reference walk takes `A = 1`, where a
negation-aware one would take `A = 2` and sit at `0.886`, making every
ratio in the table a further `√2` larger.

**The Koblitz case, which is the one the campaign cares about.**  There
both folds are `n`, so `S_family = 0.75·(cof·r)^{2/3}/(n√r)` against a
signed-Frobenius reference at `A = 2n`, `S_rho = √(π/4n)`.  The fold
buys a factor `n`, and the reference gets `√n` of it back, so the ratio
grows as `r^{1/6}/√(log r)` — slower than the prime case, and still
unbounded:

| curve | `r` | `S_family` | vs signed-Frobenius rho |
|:--|:--|--:|--:|
| `GF(2^41)`, cof 4 | `2^39` | 4.17 | 30× |
| `GF(2^83)`, cof 2 | `2^82` | 187 | 1.9×10³ |
| `GF(2^163)`, cof 2 | `2^162` | 9.8×10⁵ | 1.4×10⁷ |
| `GF(2^233)`, cof 2 | `2^232` | 2.2×10⁹ | 3.8×10¹⁰ |
| `GF(2^283)`, cof 2 | `2^282` | 5.9×10¹¹ | 1.1×10¹³ |

The first row is measured to `31.5×` against the law's `30.1×`; the rest
are extrapolation.  At `sect163k1`'s size the whole pair-table family,
folded by the full `⟨σ, −1⟩` and sized at its own optimum, is about
`10^7` times the cost of the automorphism-aware walk it is competing
with.

This is a stronger statement than the earlier rounds could make.  §10.8
said "not a crossover" about the rows it had measured.  This says the
family has none to find.

**The conclusion does not rest on the `K/2`.**  That heuristic is the
shakiest part of §11.2, and the runs show it is only roughly right: the
measured relations used to close, over columns, run from `1.13` at
`K = 5` down to `0.08` at `K = 277`, falling as the base grows.  So
suppose instead that a row needs `K^{1−a}` relations for some `a ≥ 0`
— `a = 0` is the `K/2` heuristic, larger `a` is the faster closing the
data hints at.  Redoing the minimisation gives
`F* ∝ (#E·t)^{1/(3+a)}` and, on a prime-order curve,

| relations needed | `F*` | `S_family` |
|:--|:--|:--|
| `K` (`a = 0`) | `#E^{1/3}` | `r^{1/6}` |
| `K^{3/4}` | `#E^{1/3.25}` | `r^{0.115}` |
| `K^{1/2}` | `#E^{1/3.5}` | `r^{1/14}` |
| `K^{1/4}` | `#E^{1/3.75}` | `r^{0.033}` |

`S_family ∝ r^{2/(3+a) − 1/2}`, whose exponent is positive for every
`a < 1`.  A row would have to close on a number of relations growing
slower than *any* power of `K` to escape, and a bounded number of
two-summand relations cannot pin a logarithm over `K` unknowns.  The
crossover with a flat reference therefore does not exist for any
plausible relation count; only the rate at which the gap opens depends
on `a`.

**What the law does not cover, said plainly.**  It is a statement about
the *pair-table* family, and it rests on that family's two costs: a
table of every pair, `F²/4t`, and `Θ(K)` relations from targets that
decompose with the counting probability.  An index calculus with an
**algebraic** decomposition oracle — Semaev's summation polynomials
solved by F₄ or by SAT, which is what the published ECDLP attacks use —
builds no such table, so the `F²` term is not its term and this law says
nothing about it.  §5 prices those oracles on the same systems: at the
sizes that fit here they cost `10^4` to `10^7` group-addition
equivalents per target against the probe's `0.5` to `25`, which is why
no row in this ledger uses one to finish a logarithm.  Whether their
asymptotics are better is the question this ladder cannot reach, and
naming it is the honest end of the round rather than a claim about it.

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

### 11.5 The table: every rung on the largest instance of its regime

Means over three targets, from the frozen ladder.  `vs family` is §11.2's
model at the row's own folds, and is a model, not a bound.

| regime, instance | variant | m | \|F\| | K | trials | y/c (exact) | S | was | vs rho | vs floor | vs family | ok | class |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| **prime**, `generated-24bit-10935329`, `r = 2^23.4`, `#E = 1r`, `A = 2` | generic floor `√(π/2A)` | | | | | | 0.886 | | 0.23× | 1× | | — | boundary |
| | family optimum, `k = 1` | | 222 | | | | 11.2 | | 2.84× | 12.6× | 1× | — | model |
| | Pollard rho, r-adding, counted | | | | | | 3.93 | | 1× | 4.43× | | ✓ | reference |
| | Semaev `S₃` roots | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 5,123 | | 1,304× | 5,781× | 459× | ✓ | baseline |
| | direct subtraction | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 2,037 | | 518× | 2,298× | 182× | ✓ | accounting |
| | meet in the middle | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 268 | | 68.3× | 303× | 24.0× | ✓ | engineering |
| | meet in the middle | 3 | 512 | 256 | 265 | 0.87 (0.87) | 56.8 | | 14.5× | 64.1× | 5.09× | ✓ | engineering |
| | + negation-folded table | 2 | 512 | 256 | 11,612 | 0.99 (0.99) | 248 | 268 | 63.2× | 280× | 22.2× | ✓ | engineering |
| | + walk targets | 2 | 512 | 256 | 9,873 | 0.99 (0.99) | 24.1 | 248 | 6.14× | 27.2× | 2.16× | ✓ | engineering |
| | + negation-folded table | 3 | 512 | 256 | 265 | 0.87 (0.87) | 37.0 | 56.8 | 9.42× | 41.8× | 3.31× | ✓ | engineering |
| | + walk targets | 3 | 512 | 256 | 262 | 0.89 (0.89) | 31.4 | 37.0 | 7.99× | 35.4× | 2.81× | ✓ | engineering |
| | **+ base at the family optimum** | 2 | 222 | 111 | 31,023 | 1.02 (1.02) | **14.3** | 24.1 | **3.63×** | 16.1× | 1.28× | ✓ | engineering |
| **binary**, `random-binary-n27-b845462`, `r = 2^24.4`, `#E = 6r`, `A = 2` | generic floor | | | | | | 0.886 | | 0.38× | 1× | | — | boundary |
| | family optimum, `k = 1` | | 512 | | | | 41.6 | | 18.0× | 46.9× | 1× | — | model |
| | Pollard rho, r-adding, counted | | | | | | 2.31 | | 1× | 2.61× | | ✓ | reference |
| | meet in the middle | 3 | 526 | 263 | 1,498 | 0.90 (0.90) | 206 | | 89.1× | 233× | 4.96× | ✓ | engineering |
| | `S₄` pairs-and-solve | 3 | 526 | 263 | 1,499 | 0.90 (0.90) | 119,093 | | 51,477× | 134,383× | 2,865× | ✓ | relabelling |
| | + negation-folded table | 3 | 526 | 263 | 1,498 | 0.90 (0.90) | 191 | 206 | 82.7× | 216× | 4.61× | ✓ | engineering |
| | + walk targets | 3 | 526 | 263 | 1,462 | 0.93 (0.93) | 167 | 191 | 72.1× | 188× | 4.01× | ✓ | engineering |
| | + negation-folded table | 2 | 526 | 263 | 112,779 | 1.06 (1.05) | 1,645 | 206 | 711× | 1,856× | 39.6× | ✓ | engineering |
| | + walk targets | 2 | 526 | 263 | 130,973 | 0.99 (0.98) | 51.3 | 1,645 | 22.2× | 57.9× | 1.23× | ✓ | engineering |
| | **+ base at the family optimum** | 2 | 254 | 127 | 217,845 | 0.99 (0.96) | **50.7** | 51.3 | **21.9×** | 57.2× | 1.22× | ✓ | engineering |
| **Koblitz**, `K_0 / GF(2^41)`, `r = 2^39.0`, `#E = 4r`, `A = 82` | generic floor | | | | | | 0.138 | | 0.71× | 1× | | — | boundary |
| | family optimum, `k = 41` | | 13,021 | | | | 4.18 | | 21.4× | 30.2× | 1× | — | model |
| | signed-Frobenius rho, counted | | | | | | 0.20 | | 1× | 1.41× | | ✓ | reference |
| | meet in the middle, signed-orbit columns | 3 | 5,003 | 62 | 6,086 | 1.05 (1.04) | 58.8 | | 301× | 425× | 4.05× | ✓ | advance, count |
| | + negation-folded table | 3 | 5,003 | 62 | 6,086 | 1.05 (1.04) | 50.4 | 58.8 | 257× | 364× | 3.47× | ✓ | engineering |
| | + Frobenius-folded table | 3 | 5,003 | 62 | 6,086 | 1.05 (1.04) | 45.1 | 58.8 | 231× | 326× | 10.7× | ✓ | engineering |
| | + walk targets | 3 | 5,003 | 62 | 5,110 | 1.14 (1.13) | 37.2 | 45.1 | 190× | 269× | 8.81× | ✓ | engineering |
| | two summands, folded table, walk | 2 | 5,003 | 62 | 6,714,627 | 0.91 (0.88) | 10.1 | 58.8 | 51.5× | 72.8× | 2.39× | ✓ | engineering |
| | balanced base, three summands | 3 | 20,501 | 251 | 477 | 0.74 (0.74) | 12.6 | 37.2 | 64.5× | 91.2× | 3.02× | ✓ | engineering |
| | **balanced base, two summands** | 2 | 20,501 | 251 | 1,454,383 | 1.00 (1.00) | **5.92** | 10.1 | **30.3×** | 42.8× | 1.42× | ✓ | engineering |

### 11.6 Reading the rungs: what moved, and what the round actually bought

**The base at the family optimum is the round's one clear gain, and only
where the rule of thumb was wrong.**  On the prime ladder the `⌈bits/3⌉`
rule gave 512 signed points where §11.2 puts the optimum at 222, and
moving there took the row from `24.1` to `14.3` — `6.14×` rho down to
`3.63×`, and within `1.28×` of the family's own floor.  On the binary
ladder the rule was already almost right (`512` against the law's `512`),
so the balanced row runs the bracket below it and gains almost nothing:
`51.3 → 50.7`.  The rule was not wrong on Koblitz either, where Round 2
had already sized that base.

**The cheap restart buys about a tenth, and only where restarts are
frequent.**  From the matched eight-repeat headline, geometric means with
95% intervals:

| instance | row | speedup | 95% interval | moved |
|:--|:--|--:|:--|:--|
| prime 24-bit | `mitm_m2_negfold_walk` | 1.098 | 1.010 – 1.193 | yes |
| binary `n = 27` | `mitm_m2_negfold_walk` | 1.210 | 0.944 – 1.550 | no |
| Koblitz `n = 37` | `mitm_m2_…_walk` | 1.041 | 0.653 – 1.660 | no |
| Koblitz `n = 41` | `mitm_m2_…_walk` | 0.959 | 0.710 – 1.296 | no |

One row of eight clears one.  The ladder, at three repeats over more
sizes, shows where the lever does pay: the binary `n = 21` two-summand
walk at `3.33×` (`1.47–7.51`) and `n = 24` at `1.77×` (`1.22–2.56`),
which are the rungs that restart most.

**And it costs something fixed, which the small rungs show.**  The pooled
offsets are drawn with the jumps at setup: sixteen extra `[c]G + [d]Q`,
thirty-two scalar multiplications, once per run whether or not the walk
ever restarts.  On rungs that barely restart that is a straight loss —
`K_1 / GF(2^11)` at `0.597`, `K_0 / GF(2^13)` at `0.650`, `bench-10bit`'s
three-summand walk at `0.656`.  Every one of those is below `2^20` and so
outside §1.6's window, but it is a real trade and not noise, and the fix
is obvious: draw the pool on the *first* restart rather than at setup, so
a run that never restarts never pays for it.  That is the next round's
first line, not this one's, because the code and the measurements have to
move together.

**The best row of each regime, against Round 2's:**

| regime | Round 2 | Round 3 | vs rho | vs family |
|:--|--:|--:|--:|--:|
| prime, `2^23.4` | 24.2 | **14.3** | 3.63× | 1.28× |
| binary, `2^24.4` | 71.2 | **50.7** | 21.9× | 1.22× |
| Koblitz, `2^39.0` | 5.12 | **5.92** | 30.3× | 1.42× |

The Koblitz row reads worse and did not move: the headline puts it at
`0.927` with an interval of `0.849` to `1.012`, inside its own spread.

### 11.7 A drift in the unit, found while comparing

111 of the 166 cross-round rows have **identical native counts** — same
trials, same relations, same group operations — because nothing this
round touched them.  Their cost ratio should be exactly one.  It is not:
it runs from `0.912` to `1.076`.

The cause is the unit itself.  Every non-addition is converted at a
factor measured on the host *at the start of each run*, and `ns_per_add`
came out `213 ns` for the Round-2 ladder and `146 ns` for this one.  A
row whose cost is mostly square roots or lookups is therefore repriced by
up to eight per cent between two runs that did identical work.

So the unit carries about a `±8%` run-to-run drift, and any cross-run
ratio inside that band is the unit moving rather than the method.  The
comparison file now records `native_counts_identical` per row, prints a
calibration-free ratio of group additions beside the converted one, and
states the drift's width; §11.6 claims a speedup only for rows whose
counts differ.  This is an **accounting** finding by the §3 test: no
algorithm changed, and it corrects how earlier rounds' cross-run
comparisons should be read.

### 11.8 The exponents, refitted

| regime | variant | α (r) | R² | sizes | α of the rung it was built on |
|:--|:--|--:|--:|--:|--:|
| prime | rho reference | 0.274 | 0.903 | 8 | — |
| prime | mitm_m2_negfold_walk | 0.487 | 0.958 | 8 | 0.681 |
| prime | mitm_m3_negfold_walk | 0.502 | 0.975 | 8 | 0.577 |
| prime | **mitm_m2_negfold_walk_balanced** | **0.405** | 0.961 | 8 | 0.487 |
| char2 | rho reference | 0.302 | 0.973 | 5 | — |
| char2 | mitm_m2_negfold_walk | 0.545 | 0.975 | 5 | 0.788 |
| char2 | mitm_m2_negfold_walk_balanced | 0.567 | 0.966 | 5 | 0.545 |
| koblitz | rho reference | 0.254 | 0.917 | 11 | — |
| koblitz | mitm_m2_…_frobfold_walk | 0.440 | 0.829 | 11 | 0.544 |
| koblitz | mitm_m2_…_frobfold_walk_balanced | 0.417 | 0.833 | 4 | 0.440 |

No fit is below the reference's, which is what §1.6 asks.  And the
balanced prime row's `0.405` must not be extrapolated, for the reason
§11.3 gives: the law predicts `2/3` and the row is measured inside its
convergence toward the law, not at its asymptote.

### 11.9 Against the targets of §1.6, and the `AGENTS.md` §8 gate

**Neither condition of §1.6 is met, and the round does not claim them.**
The closest row to the reference is the prime balanced walk at `3.63×`
rho on `r = 2^23.4`; no row is below `1×` at `r ≥ 2^20`.  No fitted total
exponent is below the reference's `0.25–0.30`.  No `yield/ceiling`
reaches `1.5` against the exact ceiling on a base outside a proper
subgroup.  Every logarithm on every row of all four runs was recovered
and verified as `[d]G = Q`, every rho run recovered and verified its own,
and `repeated_column_rows` and `pinned_by_repeated_row` are zero
everywhere.

**The §8 comparison** is `ic-boundary-round3-comparison-2026-09-21.json`:
166 cross-round rows paired with the same row of the previous round's
run, 116 within-run rows pairing each new row with the rung it was built
on, and 38 holdout rows on a fresh seed.  Within-run speedups run `0.092`
to `53.1` with a median of `1.327`; 22 are below one and are reported as
such, being the rungs where a lever costs more than it saves.

**The reference and the candidate were run**, as §8 asks, and for this
round that meant two binaries rather than two rows: the same command at
eight repeats on a binary built from a git worktree at the Round-2 commit
`acd192dc` and on this round's, on the same instances with the same seed,
run concurrently.  A round whose lever makes an existing row cheaper
cannot be measured any other way.

The frozen WDSat regression suite does not apply and was not run, for the
reason §10.7 gives: it measures one SAT-solver stage on sixty fixed
Weil-descended inputs and is not an adapter for another pipeline.  A
restart's step function, a base size and a derived model are not a solver
stage.  No runtime claim is made, so no paired wall-clock interval is
owed.

### 11.10 What does not count, and what this is not

- **Not a crossover.**  `3.63×` rho on the prime ladder is the closest
  any row has come, and §11.3 says the family cannot reach `1×` at any
  size.
- **Not an advance.**  Every rung is engineering by the §3 test: the base
  size is a parameter of the same counting bound, the restart is a
  cheaper way to draw a step function, and `yield/ceiling` against the
  exact ceiling does not move.
- **Not an exponent result**, and §11.3 explains why the balanced row's
  `0.405` is the opposite of one.
- **The permutation is not a measured improvement** (§11.1), and the
  pooled offsets are a measured *loss* on rungs that rarely restart
  (§11.6).  Both are on the table.
- **The `±8%` drift in the unit** (§11.7) is a correction to how every
  cross-run ratio in this note, including Round 2's, should be read.
- **Not a claim about a deployed curve**, and not a wall-clock benchmark.
  The extrapolations of §11.3 are extrapolations and are marked as such.

### 11.11 Reproducing

```bash
cargo build --release --bin ic
./target/release/ic boundary --repeats 3 --out round3.json
./target/release/ic boundary --seed 1213743172 --repeats 2 --prime-bits 22,24 \
    --char2-degrees 24,27 --koblitz-degrees 37,39,41 --s4-max-degree 24 --out holdout3.json
# The matched headline: the same command on two binaries.
git worktree add /tmp/r2 acd192dc && cargo build --release --bin ic --manifest-path /tmp/r2/Cargo.toml
for bin in ./target/release/ic /tmp/r2/target/release/ic; do
  $bin boundary --repeats 8 --prime-bits 24 --char2-degrees 27 --koblitz-degrees 37,41 \
      --s4-max-degree 20 --no-fold-max-degree 20 --out headline-$(basename $(dirname $bin)).json
done
python3 docs/ic/tools/boundary_round_compare.py round3.json \
    --across docs/ic/runs/ic-boundary-ledger-round2-2026-09-21.json \
    --baseline docs/ic/runs/ic-boundary-ledger-2026-09-21.json \
    --holdout holdout3.json --out comparison3.json
python3 docs/ic/tools/boundary_repeat_spread.py round3.json --min-repeats=3
python3 docs/ic/tools/boundary_round_note_tables.py round3.json --ladder
python3 docs/ic/tools/boundary_scoreboard_rows.py round3.json
cargo test --release --lib cryptanalysis::ic_
```

## 12. Round 4: the unit, pinned

**Frozen run:** `docs/ic/runs/ic-boundary-ledger-round4-2026-09-21.json`
— the Round-3 ladder rerun with nothing changed but the conversion.

§11.7 found the defect and this round fixes it.  It is **accounting** by
the `AGENTS.md` §3 test: the numbers move, the algorithm does not, and
no gain is claimed from it.  It is worth a round of its own because it
sets how much any later round can claim.

### 12.1 What was wrong with the unit

`S` is group-addition equivalents per `√r`, so every non-addition count
is converted through `ns_per_<unit> / ns_per_add`.  Through Round 3 both
factors were measured on the host at the start of each run.  The ratios
between them then drifted, for the *same* instance and the *same* unit,
across three ladders on one machine:

| | median | 90th | max |
|:--|--:|--:|--:|
| spread of a pinned ratio across three runs | 1.08 | 1.23 | 3.70 |

The worst are `ns_per_sqrt` on `bench-12bit` at `3.70` and
`ns_per_frobenius` on `K_0 / GF(2^41)` at `2.50`.  The consequence was
§11.7's: 111 of 166 cross-round rows had identical trials, relations and
group operations and still moved, by `0.912` to `1.076`.

`AGENTS.md` §6 keeps operation counts as the metric "because they
survive hardware".  A conversion re-measured per run does not survive
it, and an eight per cent floor is not a detail when the levers of §11.6
move rows by `1.10×` to `1.89×`.

### 12.2 What the unit is now

The ratios live in `docs/ic/calibration.json`, compiled into the binary
and applied per instance after the host measurement, so that
`Calibration::gae` — which divides by `ns_per_add` — yields exactly the
pinned ratio whatever the host's speed.  `boundary_pin_calibration.py`
regenerates the table as the **median** over a set of frozen runs, which
is what keeps one unlucky measurement from becoming the unit; it
reproduces the committed table exactly from the three ladders it was
built from.

Three properties this design keeps, and each matters:

- **The measurement survives.**  Every factor is still taken and still
  reported, as `calibration_measured`.  It is the wall-clock
  practicality note of §6, and it is how a host that no longer resembles
  the reference one becomes visible.  It prices nothing.
- **The fallback is loud.**  An instance the table does not carry — a
  new size, a freshly generated curve — keeps the host's factors, and
  `calibration_pinned` names every unit that fell back.  A row priced
  the old way says so.  Two consecutive quick runs agree to twelve
  digits on every pinned row, while the quick configuration's `n = 12`
  binary curve, which the table does not carry, still moves `1.09×`.
  That is the fallback working rather than failing quietly.
- **It prices and gates nothing else.**  Every consumer of a
  `Calibration` hands it to `price_phase`.  No native count and no
  branch depends on it, which is what makes the round accounting rather
  than a change to the method.

### 12.3 Re-pinning is itself a repricing

Regenerating the table reprices every row at once.  So the table records
the host it came from and the runs it was taken over, and a later round
that re-pins has to say so in the note and carry the delta the way this
one does — otherwise the drift returns as a step instead of as noise.
The reference host is the one named in the file.

### 12.4 The correction, measured

The ladder was rerun with nothing changed but the conversion, and
`boundary_repricing_check.py` compared every field of `group_ops` and
every key of `native`, in every phase, row by row and repeat by repeat.

| | |
|:--|--:|
| rows compared | 537 |
| identical in every native counter | **537** |
| rows whose counters moved | **0** |
| `S` after over before, median | 1.0000 |
| range | 0.9048 to 1.0766 |
| rows repriced by more than 2% | 21 of 537 |

So the round is accounting, as claimed: nothing computed anything
different, and the correction to individual rows reaches about `±9.5%`.
The checker's sensitivity is not assumed — run against two runs that
genuinely differ, the rotation ladder against the shuffle one, it finds
94 of 537 rows with moved counters and names the fields.

**Where the correction lands, and why.**  A row is repriced in
proportion to how much of its cost was *not* plain group additions,
because additions are the unit's numeraire and need no conversion.  The
share of each variant's cost that goes through a conversion, against the
largest repricing that variant saw:

| variant | converted share | worst ratio |
|:--|--:|--:|
| `semaev_s4_pairs_and_solve_m3` | 0.999 | 0.925 |
| `…_s4_pairs_and_solve_m3_signed_orbit_columns` | 0.996 | 0.964 |
| `semaev_s3_roots_m2` | 0.897 | 0.905 |
| `mitm_m2_negfold` | 0.591 | 0.995 |
| `mitm_m2_negfold_walk` | 0.247 | 0.996 |
| `mitm_m2_negfold_walk_balanced` | 0.196 | 0.997 |
| `mitm_m2_…_frobfold_walk_balanced` | 0.082 | 0.996 |
| `mitm_m3_signed_orbit_columns_negfold` | 0.011 | 1.0002 |
| `mitm_m3_signed_orbit_columns` | 0.006 | 1.0001 |

The two ends of that table are the whole story.  A row that is 99.9%
converted work moves by nine per cent; a row that is 0.6% converted work
moves by one part in ten thousand.

**And that is why §11's conclusions do not depend on the drift.**  The
best row of each regime is a walk over a pair table — additions, almost
all the way down — so the headline figures are unmoved:

| regime | best row | `S` before | `S` after | vs rho |
|:--|:--|--:|--:|--:|
| prime, `2^23.4` | balanced `m = 2` walk | 14.281 | 14.280 | 3.63× |
| binary, `2^24.4` | balanced `m = 2` walk | 50.656 | 50.652 | 21.9× |
| Koblitz, `2^39.0` | balanced `m = 2` folded walk | 5.924 | 5.926 | 30.3× |

The rows the drift *did* move are the ones this ledger already classes
baseline or relabelling — Semaev `S₃` roots and `S₄` pairs-and-solve —
and no conclusion ever rested on them.  So the defect was real, it was
worth removing, and it was not threatening the answers: it put a floor
under what a *future* round could resolve, which is the reason to fix it
rather than a correction to what earlier rounds said.

## 13. Round 5: the restart pool, drawn on first use

§11.6 ended by naming this round's first line, and this is it:

> The pooled offsets are drawn with the jumps at setup: sixteen extra
> `[c]G + [d]Q`, thirty-two scalar multiplications, once per run whether
> or not the walk ever restarts. … the fix is obvious: draw the pool on
> the *first* restart rather than at setup, so a run that never restarts
> never pays for it.

### 13.1 How much is actually lying there, counted before touching the code

The Round-4 ladder has 204 walk rows.  Their `walk_restarts` counter —
which counts the opening segment as the first — distributes like this:

| `walk_restarts` | rows | offsets the row used |
|--:|--:|--:|
| 1 | 110 | 0 |
| 2 | 10 | 1 |
| 3 | 7 | 2 |
| 4–15 | 29 | 3–14 |
| ≥ 16 | 48 | 16 |

So **110 of 204 rows drew sixteen offsets and took none**, and another 46
took fewer than they paid for.  Only the 48 rows that restart at least
sixteen times used the pool they bought.

What that waste is worth is arithmetic, not a guess, and it is worth
doing before the measurement so the measurement has something to
contradict.  One offset is two scalar multiplications and one addition;
a scalar multiplication by a uniform `k < r` is `⌈log₂r⌉ − 1` doublings
and about `log₂r/2 − 1` additions, so sixteen offsets cost about

```
ΔS_max  ≈  (48·log₂r − 48) / √r          group-addition equivalents
```

in the unit, and that is the *whole* saving on a row that takes no
offset.  It falls fast:

| `r` | `ΔS_max` |
|:--|--:|
| `2^10` | 10.4 |
| `2^15` | 3.4 |
| `2^20` | 0.89 |
| `2^27` | 0.11 |
| `2^39` | 0.0025 |

That table is the round's honest headline before a single run: this
removes a real cost, it removes it almost entirely from rungs below
`2^15`, and **it does essentially nothing inside §1.6's window**.  (That
last clause was measured and is corrected in §13.5: it is right about
`ΔS`, which is a function of `r` alone, and wrong about the *ratio*,
which also depends on the row's `S` — one row inside the window moves by
`7.2%`.)  The
three losses §11.6 recorded — `K_1 / GF(2^11)` at `0.597`,
`K_0 / GF(2^13)` at `0.650`, `bench-10bit`'s three-summand walk at
`0.656` — are all below `2^20`, and so is all of the gain.

There is one consequence that is not cosmetic.  §11.3 reads the ladder's
`S / S_family` column as the family law converging to its own optimum,
and the small rungs are where that column is measured; a fixed setup
cost that nobody uses inflates exactly those rows.  So the reason to
land this is the same as Round 4's reason for pinning the unit: it
raises the resolution of a future measurement rather than improving an
answer.

**The lever is also nearly exhausted by construction**, and that is
worth saying now rather than discovering it later.  Drawing offsets one
at a time on demand — which is what this implements — recovers the whole
of the 110 rows' waste and most of the 46's.  Nothing else remains in
the restart: what is left is sixteen jumps at setup, which every segment
uses, and one addition per restart, which is already the floor for
moving off a path.

### 13.2 The comparison is exact, by construction

A saving this small is below the run-to-run spread of §11.4, so a
before-and-after of two ladders would not resolve it.  Two changes make
it resolvable without any statistics at all:

- **The offsets come off their own random stream.**  `pool_rng` is
  seeded from the run's seed and is touched by nothing else, so moving
  the draw earlier or later cannot shift the walk's randomness.
- **A restart takes offset `k mod 16`**, not a random one, so selecting
  an offset consumes no randomness either.

With both, the eager and the lazy arm walk **bit-identical
trajectories**: same steps, same restarts, same trials, same relations,
same rows, same logarithm.  They differ in exactly one quantity — what
the offsets cost — and the run reports that quantity directly, as
`walk_pool_ops`, so the difference between the arms is a subtraction
rather than an inference.  Both arms come from one binary on one host;
`--eager-restart-pool` selects the baseline.

**The eager arm is Round 4's pool *policy*, not Round 4's run**, and the
distinction matters.  Round 4 chose an offset with a draw from the main
stream, which this round replaces with the fixed cycle; that alone
shifts every walk's randomness, so the Round-4 file is not a valid
baseline for a saving this small and is not used as one.  It stays
frozen for what it was.  Where a cross-round figure is wanted — §11.6's
three losses were measured against *Round 2* — it is quoted as a
cross-round figure, with the §12 repricing caveat, and never as the
paired result.

The test `drawing_the_restart_pool_lazily_changes_the_cost_and_nothing_else`
asserts the whole of that on two prime rungs and eight seeds.

### 13.3 The falsification target, declared before the ladders were read

The round succeeds only if **all** of these hold:

1. **The pairing holds.**  Every row of the two ladder arms is identical
   in every native counter except `walk_jumps`, `walk_pool_offsets` and
   `walk_pool_ops`; zero rows have a moved trajectory.
2. **The saving accounts for itself exactly.**  On every walk row,
   `gae_eager − gae_lazy` equals `walk_pool_ops_eager − walk_pool_ops_lazy`
   to the operation — no row saves more or less than the offsets it
   stopped drawing.
3. **The loss is gone and no new one appears.**  `S` falls on each of
   the three rows §11.6 named, and rises on none anywhere.
4. **Correctness is preserved.**  `all_verified` on both arms, every
   recovered logarithm equal to the planted one, and
   `repeated_column_rows` and `pinned_by_repeated_row` zero everywhere,
   on the holdout seed as well.
5. **It is engineering and is reported as engineering.**  Trials and
   yield against the exact ceiling must be *identical*, since the change
   touches no decision the search makes.  The fitted exponents may move,
   because the saving is larger at small `r` than at large: a steepening
   of the walk rows' fit is expected and must be reported as an artifact
   of removing a size-independent cost, never as a change in the
   method's growth.

**Inadmissible**: changing the factor base or the decomposition size;
changing the operation accounting or the pinned conversion; skipping
verification; choosing favourable seeds; quoting the relation phase
instead of the whole-pipeline `S`; or reporting a small-rung saving as
if it moved a conclusion inside §1.6's window.

**Abandon if** the pairing fails — if any trajectory counter moves, the
change did more than move a cost and this design is wrong.

### 13.4 The pairing, checked

Six runs, all from one binary (`blake3 ae6dd73f4404`) on one host: the
ladder, a holdout seed and the eight-repeat headline, each in both arms.
`boundary_pool_pairing_check.py` compares every native counter and every
group-operation count of every row, per repeat.

| | ladder | holdout | headline |
|:--|--:|--:|--:|
| rows compared | 537 | 106 | 240 |
| identical in every counter but the pool's | **537** | **106** | **240** |
| rows whose trajectory moved | **0** | **0** | **0** |
| walk rows | 204 | 48 | 112 |
| …that drew fewer offsets | 155 | 31 | 69 |
| rows whose saving ≠ the offsets they dropped | **0** | **0** | **0** |
| rows whose `S` rose | **0** | **0** | **0** |
| `all_verified` | yes | yes | yes |

So the design holds exactly: the two arms did the same work, row for
row and repeat for repeat, and the only thing that changed is what the
offsets cost.  The saving is a subtraction.

### 13.5 Where the saving lands, and a correction to §13.1

The closed form of §13.1 is a good predictor.  Per instance, the largest
saving any walk row realised against `(48·log₂r − 48)/√r`:

| instance | `log₂r` | predicted `ΔS` | measured `ΔS` | measured / predicted |
|:--|--:|--:|--:|--:|
| `K_0 / GF(2^15)` | 9.55 | 14.980 | 13.854 | 0.925 |
| `bench-10bit` | 9.68 | 14.531 | 13.490 | 0.928 |
| `K_1 / GF(2^11)` | 9.95 | 13.651 | 12.442 | 0.911 |
| `bench-14bit` | 13.99 | 4.893 | 4.717 | 0.964 |
| `bench-20bit` | 20.00 | 0.891 | 0.850 | 0.954 |
| `K_1 / GF(2^23)` | 22.00 | 0.492 | 0.478 | 0.972 |
| `generated-24bit` | 23.38 | 0.325 | 0.318 | 0.980 |
| `K_0 / GF(2^37)` | 27.78 | 0.085 | 0.082 | 0.969 |
| `K_0 / GF(2^41)` | 39.00 | 0.002 | 0.002 | 0.990 |

The formula is high by 1 to 9 per cent across the whole ladder, as it
should be: it prices sixteen offsets and the best row on each instance
still took a few.

**§13.1 overstated one thing and it is corrected here.**  It said the
saving "does essentially nothing inside §1.6's window".  In absolute
terms that is right — `ΔS ≤ 0.85` at `2^20` and `0.002` at `2^39` — but
`ΔS` is a function of `r` alone while the *ratio* depends on the row's
own `S`, and the Koblitz rows are the cheapest on the board because of
the `2n` fold.  Splitting the walk rows at the window's edge:

| | walk rows | best ratio | largest `ΔS` |
|:--|--:|--:|--:|
| below `2^20` | 120 | 0.5825 | 14.633 |
| at or above `2^20` | 84 | **0.9283** | 0.850 |

That best is `K_1 / GF(2^23)`'s two-summand folded walk at `2^22`: a
**7.2 per cent** cut, inside the window, on a row whose `S` is about
`6.7`.  Small, but not nothing, and the prediction should have said so.
The prime and binary rows at the same sizes move by 1 per cent or less,
because their `S` is five to eight times larger for the same `ΔS`.

### 13.6 The loss §11.6 recorded is closed

§11.6 measured the eager pool as a straight loss on the rungs that
barely restart, against Round 2.  The eager arm reproduces those losses
and the lazy arm removes them:

| row | Round 2 `S` | eager `S` | lazy `S` | lazy vs Round 2 | §11.6 quoted |
|:--|--:|--:|--:|--:|--:|
| `K_1 / GF(2^11)` `m = 2` folded walk | 19.22 | 31.66 | 19.21 | **1.001×** | 0.597 |
| `K_0 / GF(2^13)` `m = 2` folded walk | 18.54 | 28.54 | 18.53 | **1.000×** | 0.650 |
| `bench-10bit` `m = 3` walk | 26.31 | 39.80 | 26.31 | **1.000×** | 0.656 |

The eager arm's own ratios to Round 2 are `0.607`, `0.649` and `0.661`
against §11.6's `0.597`, `0.650` and `0.656` — an independent
reproduction of that measurement two rounds later, on a different random
stream.  The lazy arm lands on Round 2's figure to three or four
significant figures, which is what it should do: on a run that never
takes an offset, drawing the pool lazily means drawing nothing, and the
row pays exactly the sixteen jumps Round 2 paid.

These are cross-round comparisons and carry §12's repricing caveat.  It
is negligible here: all three are addition-dominated rows, where §12
measured the correction at under one part in a thousand.

### 13.7 The re-plumbing is a re-randomisation, and it is bigger than the lever

The eager arm is not the Round-4 file, and the gap between them is the
most useful number this round produced.  Against Round 4, the eager arm
(same policy, offsets moved to their own stream and selected by a cycle)
gives:

| | |
|:--|--:|
| rows compared | 537 |
| identical in every counter | 333 |
| rows whose counters moved | 204 |
| distinct variants among them | 68, **every one a walk row** |
| `S` after over before, median | 1.0000 |
| range | **0.5394 to 1.8308** |
| moved by more than 2% | 131 of 537 |

The structure is exactly right: every non-walk row is bit-identical to
Round 4, so nothing outside the walk changed; every walk row moved,
because its randomness did.  The median is `1.0000` and the geometric
mean over variant means is `0.9961`, so the re-seed is unbiased, as an
algorithm-preserving change must be.

But the *spread* is `0.54` to `1.83` at three repeats, while the lever
this round set out to measure is at most `7%` inside the window.  **A
before-and-after of two ladders could not have seen it.**  That is not a
retrospective justification for the paired design; §13.2 required it in
advance, and this is the measurement that says by how much.

It is also a warning about reading any single walk row across rounds.
The best row of each regime, Round 4 against this round:

| regime | best row | Round 4 | eager | lazy | vs rho | lazy / eager |
|:--|:--|--:|--:|--:|--:|--:|
| prime, `2^23.4` | balanced `m = 2` walk | 14.280 | 13.966 | **13.966** | 3.55× | 1.00000 |
| binary, `2^24.4` | balanced `m = 2` walk | 50.652 | 36.929 | **36.929** | 15.96× | 1.00000 |
| Koblitz, `2^39.0` | balanced `m = 2` folded walk | 5.926 | 5.117 | **5.115** | 26.13× | 0.99952 |

**None of those moves is this round's lever**, and the last column says
so: the lever is exactly `1.000` on all three, because a row that
restarts sixteen times or more draws the whole pool either way.  The
binary row reading `36.9` instead of `50.7` is the same algorithm on a
different random stream, and it sits inside the `0.54–1.83` spread above.
Classed by §3 that move is **accounting**: numbers changed, the algorithm
did not.  It is carried onto the page as the current measurement, with
Round 4's figure as its before mark and this paragraph as its
explanation, because hiding it would be worse — but it is not progress
and this ledger does not count it as any.

### 13.8 Against §13.3's target, and the class

| condition, as declared | result |
|:--|:--|
| 1. every row identical but for the three pool counters, zero trajectories moved | **met** — 883 rows over three run pairs, zero |
| 2. the saving equals the offsets dropped, to the operation | **met** — zero mismatches |
| 3. `S` falls on §11.6's three rows, rises nowhere | **met** — `0.607`, `0.649`, `0.661`; zero rows rose |
| 4. correctness preserved, holdout included | **met** — `all_verified` on all six runs, `repeated_column_rows` and `pinned_by_repeated_row` zero throughout |
| 5. trials and yield identical; any exponent move reported as an artifact | **met, with the move reported below** |

On the fifth: the trials and the yield against the exact ceiling are
inside "identical in every counter", so they did not move at all.  The
fitted exponents did, in the direction and for the reason §13.3 named —
a size-independent cost removed from every row makes the small rungs
cheaper by more than the large ones, so the slope rises:

| regime | row | eager `α` | lazy `α` |
|:--|:--|--:|--:|
| prime | `mitm_m2_negfold_walk` | 0.489 | 0.537 |
| prime | `mitm_m2_negfold_walk_balanced` | 0.408 | 0.447 |
| prime | `mitm_m3_negfold_walk` | 0.507 | 0.554 |
| binary | `mitm_m2_negfold_walk` | 0.533 | 0.542 |
| Koblitz | `mitm_m2_…_frobfold_walk` | 0.433 | 0.459 |
| Koblitz | `mitm_m2_…_frobfold_walk_balanced` | 0.411 | 0.431 |

Every walk row steepens, by `0.01` to `0.05`, at equal or better `R²`.
**This is an artifact of removing a fixed cost, not a change in the
method's growth**, and it must not be read as one: the underlying
algorithm is bit-identical between the two columns.  If anything the
lazy column is the more honest fit, because the eager one was measuring
a constant that no longer exists.

**The class, by §3's test: engineering.**  `S` fell, on 155 of 204 walk
rows; the trials, the yield against the ceiling and every other count
did not move at all; the ratio to the counting boundary is unchanged in
structure.  The re-plumbing that made the measurement possible is
**accounting** and is labelled separately, per §13.7.  Neither is an
advance and neither is reported as one.

### 13.9 What does not count, and what this is not

- The saving is confined to small `r`.  Nothing here moves the family
  law of §11.2, the crossover extrapolations of §11.3, or the verdict
  that no variant of this family crosses rho at any size.
- The binary regime's headline reading `36.9` where Round 4 read `50.7`
  is re-randomisation (§13.7), not a gain, and is classed accounting.
- The exponent steepening of §13.8 is an artifact of removing a
  size-independent constant.  It is not evidence about growth.
- The cross-round figures of §13.6 are cross-round and carry §12's
  caveat; only the eager-against-lazy columns are paired.
- `ΔS_max` in §13.1 is a derivation, checked against measurement in
  §13.5; the rows of §13.5's right-hand column are the measurement.

**One gap this round surfaced and did not close.**  The holdout's prime
22-bit curve is generated from the holdout seed, so it is a different
curve from the ladder's and is not in `docs/ic/calibration.json`; its
seven non-addition units price at factors measured on the host, and 32
of the holdout's rows therefore differ between two runs that did
identical work.  The pairing check reports those separately rather than
counting them, so no figure here rests on them.  The fix is to price an
unpinned instance at its regime's median ratio instead of the host's
measurement, which makes the unit host-independent everywhere rather
than only on the default ladder.  That is a repricing with its own
classification and its own before-and-after, so it is the next round's
line, not this one's.

### 13.10 Reproducing

```
# One binary, one host, two configurations.
cargo build --release --bin ic
ic boundary --out round5-eager.json --eager-restart-pool
ic boundary --out round5.json

# The holdout seed and the eight-repeat headline, both arms.
ic boundary --seed 1213743172 --repeats 2 --prime-bits 22,24 \
  --char2-degrees 24,27 --koblitz-degrees 37,39,41 [--eager-restart-pool] --out ...

ic boundary --repeats 8 --prime-bits 24 --char2-degrees 27 \
  --koblitz-degrees 37,41 [--eager-restart-pool] --out ...

# The pairing, per run pair.
python3 docs/ic/tools/boundary_pool_pairing_check.py round5-eager.json round5.json

# The re-plumbing against Round 4, which is not a pairing.
python3 docs/ic/tools/boundary_repricing_check.py \
  docs/ic/runs/ic-boundary-ledger-round4-2026-09-21.json round5-eager.json
```

The six runs record `git_commit` as `c730cccb`, the commit that was
checked out when they started: provenance is captured at start-up (§11.9)
and the code they ran was committed afterwards.  `binary_blake3`
(`ae6dd73f4404…`) is the identifier that matters and is the same on all
six.

## 14. The decomposition systems, in Petit–Quisquater's shape

§5 priced the decomposition oracles per target and §11.3 left the
algebraic ones as the question this ladder cannot reach: the family
shape law bounds the *tabled* family, and an algebraic oracle builds no
table, so its asymptotics sit outside everything §11 established.  This
section starts measuring them, in the shape Petit and Quisquater's
Table 2 uses (*On Polynomial Systems Arising from a Weil Descent*,
ASIACRYPT 2012, p. 461).

Their table reports, per `(curve family, n, n', m)` cell, the average
maximal degree a Gröbner basis reached, the average time and the peak
memory — and its point is not the timings.  It is that in every cell
the degree reached came out **below** the bound derived for a generic
system, because Semaev's polynomials are sparse and the bound is not.

### 14.1 This is a stage diagnostic, and says so first

Everything in this section prices **one decomposition oracle call on one
target**.  By §2 that is never a speed, and no row here may be quoted as
one.  The whole-pipeline unit `S`, the floor and the rho reference stay
where they are, in §3 and §10–§13.  What this section can do is give the
algebraic oracle a boundary of its own to be measured against, in the
way rho's `S ≈ 1.3` is the boundary for the whole method.

### 14.2 What is measured, and what is not reproduced

The descent is the plain one: `x_i = Σ_k v_{i,k} e_k` over a subspace
`V ⊂ F_{2^n}` of dimension `n'`, substituted into `S_{m+1}(x_1, …, x_m,
x_R)` and split into `n` boolean equations in `m·n'` unknowns.  `n' =
⌈n/m⌉` makes that square.  The engine is the repository's boolean-ring
Buchberger over `F_2[v]/(v²−v)`.

**Their numbers are not reproduced and are not claimed to be.**  Table 2
solves a *symmetrised* system in `mt + 1 = m² + 1` variables — five at
`m = 2`, ten at `m = 3` — from the block structure of their Section 4.
At `(n, n', m) = (11, 6, 2)` that is a five-variable system where this
one has twelve.  The degrees do not compare row by row.  What carries
over is the shape of the measurement: a derived degree bound in one
column, the degree reached beside it, and the ratio.

### 14.3 Two boundaries, both derived before measuring

**The degree bound** is the semi-regular degree.  For equations of
degrees `d_1, …, d_k` in `v` boolean variables, the Hilbert series of a
semi-regular quotient is `(1 + t)^v / Π_i (1 + t^{d_i})`, and the degree
of regularity is the index of its first non-positive coefficient
(Bardet–Faugère–Salvy).  It is what a system with no exploitable
structure reaches, so a measured degree below it is structure the solver
found.  The series is short enough to check by hand and a unit test does:
`(1+t)^4 / (1+t²)² = 1 + 4t + 4t² − 4t³ + …`, first non-positive at
three.

**The reference** is what §1 asks for — the cost of the best algorithm
that already solves the same problem, in the same unit.  At these sizes
that is exhaustive search over the subspace: evaluate every equation at
every point, `2^{m·n'} · Σ_i |terms_i|` monomial tests.  Any algebraic
oracle worth the name has to get below that line before its degree
behaviour matters at all.

### 14.4 The table

Eight targets per cell, seeded; `K` is the Koblitz curve
`y² + xy = x³ + a x² + 1` and `R` a random binary curve of the same
degree.

| E | n | n' | m | vars | eqs | D_av | D_pair | D_sr | D_av/D_sr | ops | enumerate | ops/enum | ms | KiB | no decomp |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| K | 7 | 4 | 2 | 8 | 7 | 3.0 | 3.0 | 3 | 1.00 | 1.159e4 | 1.830e4 | **0.6×** | 0.4 | 3 | 3/8 |
| K | 9 | 5 | 2 | 10 | 9 | 3.0 | 3.0 | 4 | 0.75 | 8.198e4 | 1.627e5 | **0.5×** | 3.6 | 9 | 2/8 |
| K | 11 | 6 | 2 | 12 | 11 | 3.0 | 3.0 | 4 | 0.75 | 4.995e5 | 9.564e5 | **0.5×** | 36.8 | 24 | 3/8 |
| K | 13 | 7 | 2 | 14 | 13 | 3.2 | 3.2 | 4 | 0.81 | 2.778e6 | 6.642e6 | **0.4×** | 303.0 | 65 | 4/8 |
| R | 7 | 4 | 2 | 8 | 7 | 3.0 | 3.0 | 3 | 1.00 | 1.098e4 | 1.907e4 | **0.6×** | 0.4 | 3 | 3/8 |
| R | 9 | 5 | 2 | 10 | 9 | 3.0 | 3.0 | 4 | 0.75 | 8.547e4 | 1.647e5 | **0.5×** | 3.9 | 8 | 3/8 |
| R | 11 | 6 | 2 | 12 | 11 | 3.0 | 3.0 | 4 | 0.75 | 4.387e5 | 9.728e5 | **0.5×** | 31.4 | 24 | 2/8 |
| R | 13 | 7 | 2 | 14 | 13 | 3.2 | 3.2 | 4 | 0.81 | 2.683e6 | 6.789e6 | **0.4×** | 298.0 | 65 | 3/8 |

Reading it:

- **The solving degree is flat in `n`.**  `D_av = 3.0` from `n = 7` to
  `n = 11`, `3.2` at `n = 13`, against a bound that rises from 3 to 4.
  Six of the eight cells sit below the bound: the phenomenon their table
  reports, on a different system and a different engine.  It is also,
  coincidentally, their own `D_av = 3.0` for the `m = 2` rows.
- **The Koblitz curve and the random one behave identically**, on every
  column.  Whatever the extra automorphism buys elsewhere in this
  ledger, it does not change the shape of this system.
- **The count is below enumeration's count and the margin grows**:
  `0.6×`, `0.5×`, `0.5×`, `0.4×`.  **Corrected in §15.4:** the two
  counts are in different units — a Buchberger monomial operation
  measures about `50 ns` on the calibration host and an enumeration
  test `0.07 ns` — and once each is converted at its own measured
  factor, as §2 requires, the engine costs about **`75×`** enumeration
  at `n = 13` rather than `0.4×` of it.  The degree reading stands; the
  cost reading was a ratio of unlike units and is withdrawn.

**`D_pair` is in the table to keep an earlier mistake visible**, and it
has a second story now.  The first version of this measurement reported
the highest degree of a pair *processed*, which is a property of
Buchberger's selection strategy and not of the system, and is not what
an F4 run reports.  Under the coprime criterion alone it read `5.0`
where the solving degree read `3.0`, and gave a ratio of `1.25` against
the bound — the opposite conclusion, from the wrong statistic.  Under
the chain criterion the two columns **coincide** at every cell: the
pairs the old engine was pushing to degree five were exactly the ones
that contribute nothing, so pruning them removes the gap the wrong
statistic was measuring.

**Three summands, where the verdict flips.**  The same measurement at
`m = 3`, four targets a cell under a 120-second per-target budget:

| E | n | n' | m | vars | eqs | D_av | D_pair | D_sr | D_av/D_sr | ops | enumerate | ops/enum | ms | KiB | no decomp |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| K | 7 | 3 | 3 | 9 | 7 | 7.0 | 7.0 | 7 | 1.00 | 3.633e7 | 4.291e5 | 84.7x | 19585.5 | 441 | 3/4 |
| K | 9 | 3 | 3 | 9 | 9 | 7.0 | 7.0 | 7 | 1.00 | 3.606e7 | 5.268e5 | 68.4x | 19481.0 | 430 | 2/4 |
| K | 11 | 4 | 3 | 12 | 11 | 9.0 | 9.0 | 8 | 1.12 | 4.313e9 | 2.281e7 | 189.1x | 120049.4 | 24567 | 3/4 |
| ^ | | | | | | — | — | — | — | — | — | — | — | — | 4 of 4 runs hit the budget: every figure on this row is a lower bound |
| R | 7 | 3 | 3 | 9 | 7 | 7.0 | 7.0 | 7 | 1.00 | 3.531e7 | 4.572e5 | 77.2x | 19451.3 | 441 | 3/4 |
| R | 9 | 3 | 3 | 9 | 9 | 7.0 | 7.0 | 7 | 1.00 | 3.876e7 | 5.728e5 | 67.7x | 18729.6 | 435 | 3/4 |
| R | 11 | 4 | 3 | 12 | 11 | 9.0 | 9.0 | 8 | 1.12 | 4.305e9 | 2.341e7 | 183.9x | 120031.9 | 24528 | 1/4 |
| ^ | | | | | | — | — | — | — | — | — | — | — | — | 4 of 4 runs hit the budget: every figure on this row is a lower bound |

At nine variables the solver finishes and costs **68× to 85× the
enumeration it is competing with** — the opposite of the two-summand
rows, which come in at `0.4×`.  At twelve variables it does not finish:
every run hit the budget, so those rows are lower bounds and are marked
as such.  Whatever the three-summand descent buys in relations per
target, this engine does not get it back on the decomposition, and the
`m = 2` result does not carry over.

### 14.5 The engine was the measurement, twice

Neither of the two corrections above was a tuning choice; both were
defects found by taking the measurement.

**The pruning that was documented but not implemented.**  The module
header of `pq_groebner_f2` claimed Buchberger with Gebauer–Möller
pruning; the code had only the coprime criterion.  Adding the chain
criterion — if another basis element's leading monomial divides the
pair's lcm and both of its pairs have left the queue, that S-polynomial
cannot contribute — gives, paired on the same seed and the same systems:

| cell | ops before | ops after | ratio | ops/enum before → after |
|:--|--:|--:|--:|:--|
| `K n = 7` | 1.193e5 | 1.159e4 | 0.097× | 6.5× → 0.6× |
| `K n = 9` | 1.584e6 | 8.198e4 | 0.052× | 9.7× → 0.5× |
| `K n = 11` | 2.654e7 | 4.995e5 | 0.019× | 27.8× → 0.5× |
| `K n = 13` | 2.884e8 | 2.778e6 | 0.010× | 43.4× → 0.4× |

**and it reverses the verdict.**  Under the coprime criterion alone the
Gröbner basis measured 6.5× to 43× *worse* than enumeration, widening
with `n`; under both it measures 0.4× to 0.6×, narrowing.  The first
reading was a statement about the engine and was withdrawn, not filed.
(Both ratios are in raw counts of unlike units; the tenfold-to-
hundredfold cut in the engine's own count is real, and §15.4 gives the
converted ratio to enumeration, which is `75×` the other way.)
The coprime-only ladder stays frozen beside the fixed one as the before
mark, because deleting it would hide the size of the correction.

That the pruning is exact is checked **on the systems being measured**
rather than inferred from the sixty-seven tests over `pq_*`,
`koblitz_groebner` and `polynomial_reuse` that also still pass: for six
targets the variety of the pruned basis is compared against every point
of the subspace evaluated in the original equations, and must agree,
with a guard that not every target was inconsistent.

**The cell that never returned.**  At `m = 3` and twelve variables the
engine ran for four hours and forty-nine minutes without finishing, and
was killed.  A Gröbner basis has no natural stopping point, so the
engine now takes a wall-clock budget: a run that exceeds it stops,
`timed_out` marks it, and the row says every figure on it is a lower
bound.  A table with an honest "did not finish" row is worth more than
one with a silently missing cell.

### 14.6 What does not count

- No row here is a speed, a crossover, or evidence about any deployed
  curve; the largest field is `GF(2^15)`.
- `ms` and `KiB` are practicality notes.  The metric is the monomial
  operation count, as §6 requires — and, as §15.4 found, a count in one
  engine's unit is not comparable with a count in another's until each
  is converted at a measured factor.
- The `ops/enum` column is a ratio to *this repository's* enumeration on
  *these* systems, **in raw counts of unlike units**; §15.4 gives the
  converted ratio, which is the other way round by two orders of
  magnitude.  It is not a statement about F4, about Magma, or about the
  symmetrised systems Petit–Quisquater solve.
- A row whose `timed_out` is non-zero is a statement about the engine's
  budget and not about the system's difficulty.
- The degrees are not comparable to Petit–Quisquater's row by row, for
  the reason §14.2 gives.

### 14.7 Reproducing

```
ic descent --cells 7:4:2,9:5:2,11:6:2,13:7:2 --targets 8 \
  --out docs/ic/runs/ic-descent-degrees-2026-09-22.json

# The before mark, on the coprime criterion alone, is frozen at
# docs/ic/runs/ic-descent-degrees-coprime-only-2026-09-22.json
```

## 15. The solver inside `S`: the framework's first frozen sweeps

§14 priced the algebraic oracle's systems on their own and said, first,
that no row of it was a speed.  The benchmarking framework
(`docs/ic/FRAMEWORK.md`) now lets a *whole run* choose its
polynomial-system engine, with the engine's work priced into the
decomposition phase and so into `S`.  This section freezes the first
two sweeps it ran, and one of them corrects §14.

### 15.1 The boundaries, stated first

- **The reference** is the pair-table oracle on the *same base*, same
  instance, same walk, same seed: `mitm[m=2]` on
  `binary-subspace[dimension=6]`.  It is the best thing already solving
  the same decomposition problem in the same unit, and the algebraic
  rows have to reach it before they are worth anything.
- **The method boundary** is rho, which sits at `S ≈ 1.3` on every rung
  of this ledger (§3).  The sweep did not run rho on this instance, so
  the `vs rho` column of the frozen table is empty and says so; no
  ratio to rho is quoted below.
- **The degree boundary** is §14.3's semi-regular degree, `D_sr = 4` for
  thirteen quadratics in twelve unknowns.

Everything here is a toy: `E(F_{2^13})`, `r = 4091`, `log₂ r = 12`.  It
is a measurement of a *capability* — the solver reaching `S` — and of
an accounting error, not of an attack.

### 15.2 One base, four oracles, whole runs

Frozen at `docs/ic/runs/ic-bench-solver-engines-2026-09-22.json`; two
planted targets per configuration, the same two for every row.
`random-binary-n13-b1503`, `binary-subspace[dimension=6]`: 60 signed
points, 30 columns, walk targets, `incremental-gauss`.  `S` and the
ratio to the reference are per target; the solver columns are totals
over the run.

| oracle | engine | trials | rows | hit | S | S / reference | solver calls | solver ops (unit) | ns / op | solver GAE | D_solve | D_sr | priced by | correct |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| pair table | — | 126 / 112 | 26 / 20 | 0.206 / 0.179 | 40.5 / 39.7 | 1 | — | — | — | — | — | — | — | yes / yes |
| descent-algebraic | exhaustive | 126 / 112 | 26 / 20 | 0.206 / 0.179 | 1,049 / 947 | 25.9× / 23.9× | 124 / 111 | 1.61e8 / 1.43e8 (monomial tests) | 0.074 | 6.62e4 / 5.98e4 | — | 4 | measured | yes / yes |
| descent-algebraic | buchberger-f2 | 126 / 112 | 26 / 20 | 0.206 / 0.179 | 77,765 / 70,202 | 1,920× / 1,768× | 124 / 111 | 1.70e7 / 1.59e7 (monomial operations) | 52.7 / 51.0 | 4.97e6 / 4.49e6 | 2.90 / 2.95 | 4 | measured | yes / yes |
| descent-algebraic | sat-cdcl | 135 / 112 | 28 / 20 | 0.207 / 0.179 | 188,708 / 164,302 | 4,659× / 4,139× | 133 / 111 | 1.25e5 / 1.10e5 (conflicts) | 17,446 / 17,142 | 1.21e7 / 1.05e7 | — | 4 | measured | yes / yes |

The phases behind the reference row, per target: base `3.9`, pair table
`1,830`, relations `737 / 687`, matrix `5.0 / 3.4`, verification `15`,
total `2,591 / 2,539` GAE.  The algebraic rows build no table (their
`setup` is `0`) and spend it all in the relation phase, which is the
solver: `4.97e6` of the Buchberger row's `4.97e6`.

Reading it:

- **The solver reaches `S`.**  That is the capability, and it is what
  §14.1 could not do.  Every algebraic row is a complete logarithm with
  the engine's work inside its `S`, and every row is verified.
- **On the same base the pair table wins by three to four orders of
  magnitude**, and it is not close: `1,900×` for the Gröbner engine,
  `4,000×` for CDCL, `24×` even for exhaustive search over the
  subspace.  A table of `C(60, 2)` additions is `1,830` GAE once; a
  system per target is `40,000` GAE per target for Buchberger.
- **Exhaustive search is the best algebraic engine here by `75×`**, in
  the unit — the point §15.4 makes.
- **The degree reading of §14 stands and does not help.**  `D_solve =
  2.9` against `D_sr = 4`, a ratio of `0.73`, the same "below the
  bound" §14.4 reports, on the engine that is `1,900×` the reference.
  The degree column and the `S` column disagree about which engine is
  good; the `S` column is the one that is a speed.
- **The oracles agree target by target.**  The Buchberger and
  exhaustive rows walked exactly the pair table's targets and found
  exactly its relations (`126 / 26`, `112 / 20`); the test suite checks
  the same thing on 150 targets directly.  CDCL's first repeat walked
  `135` targets to `28` rows: a target with two decompositions lifted
  to a different pair, and the matrix pinned nine rows later.  Same
  answer.
- **One system per run lifted to nothing** (`unliftable_systems = 1` on
  every algebraic row, `lift_failures` `1`–`3`).  A summation
  polynomial vanishes over the algebraic closure, so a solution may
  name a twist abscissa the base never holds.  It is counted, never
  folded into "did not decompose", and it is why the algebraic hit
  rates equal the pair table's rather than exceed them.

### 15.3 The relation matrix, on its own

Frozen at `docs/ic/runs/ic-bench-relation-matrix-2026-09-22.json`:
`bench-20bit`, `r = 1,046,999`, `mitm[negation_folded=1]`, walk
targets, two planted targets, both eliminations on the same relations,
pinned calibration.

| base | matrix | rows | rank | row_ops | matrix GAE | S | correct |
|:--|:--|--:|--:|--:|--:|--:|:--|
| 64 abscissae | incremental-gauss | 38 / 35 | 38 / 35 | 469 / 382 | 8.11 / 6.60 | 9.125 / 12.040 | yes / yes |
| 64 abscissae | structured-gauss | 38 / 35 | 38 / 35 | 433 / 358 | 7.49 / 6.19 | 9.124 / 12.039 | yes / yes |
| 256 abscissae | incremental-gauss | 151 / 110 | 151 / 110 | 2,105 / 1,040 | 36.4 / 18.0 | 67.086 / 66.880 | yes / yes |
| 256 abscissae | structured-gauss | 151 / 110 | 151 / 110 | 1,869 / 936 | 32.3 / 16.2 | 67.082 / 66.878 | yes / yes |

Same rows, same rank, same logarithm; pivoting on the lightest column
saves `6`–`11 %` of the `row_ops` on two-summand rows, and the matrix
is `0.05`–`0.09 %` of `S` at these sizes, so `S` moves in the fourth
digit.  **Engineering**, bounded, and not a finding — the matrix was
`0.002 %` of the cost at §11's sizes and is a lever only where the
relation count is.  (The 256-abscissa base is past the family optimum
on purpose, to give the matrix rows to work on; its pair table is
`65,792` of its `68,644` GAE.)

### 15.4 The accounting error, and what it corrects

The first freeze of the solver sweep priced the Buchberger engine's
"monomial operations" at the pinned word-XOR ratio, on the reasoning
that a monomial is a 64-bit mask.  Its calibration record says what
that would have done: a word XOR measures `0.44 ns` on the host, a
Buchberger monomial operation `51`–`53 ns` (`wall / ops`), an
enumeration test `0.074 ns`.  Pricing the middle one at the first one's
ratio would have made the Gröbner row `120×` cheaper than it is — §6's
"changing the unit", in one line of code, in the flattering direction.

The rule is now the narrow one: only a unit the calibration table
carries a ratio for is priced by count (`word XORs`, the matrix-F4 unit
of §5), every other engine unit is priced at the host's measured wall
time over its addition time and the row says `measured`, and the
report records `ns_per_op` so the count can be re-priced.  It is
host-dependent, as §12 says every measured factor is, and it is the
honest number.

It also corrects §14.4.  That table's `ops/enum` column divided a count
of Buchberger monomial operations by a count of enumeration tests and
read `0.4`–`0.6×`; the two are unlike units by a factor of `700`, and
converted at their measured factors the engine costs `4.97e6` GAE where
enumeration costs `6.62e4` — **`75×` enumeration, not `0.4×` of it**.
The degree statistic, which was the point of §14, is unaffected.  The
cost reading is withdrawn there and here.

By §3 this is **accounting**: numbers changed, no algorithm did, and no
gain is claimed.  The reading it replaces — a count below enumeration's
count while the converted cost was two orders of magnitude above it —
is the **relabelling** pattern, work moved into the price of an
"operation" where the count did not look.

### 15.5 Classified

| row | class | why |
|:--|:--|:--|
| descent-algebraic, any engine, at `n = 13` | measurement | new rows; `24×`–`4,700×` the pair table on the same base, no boundary crossed |
| the pair table charged to the oracle's `setup`, not the base | accounting | the same base read `1,830` beside one oracle and `0` beside another |
| solver units priced by measured factor, `ns_per_op` recorded | accounting | the correction of §15.4 |
| §14.4's `ops/enum` | relabelling → withdrawn | a ratio of unlike units |
| structured-gauss | engineering | `6`–`11 %` of `row_ops`, `10⁻⁴` of `S` |

### 15.6 What does not count

- Nothing here is a speed against rho: rho was not run on this instance
  and the table's `vs rho` column is empty.
- The algebraic rows' prices are `measured`, so their `S` is this
  host's and a ratio between two `measured` rows from different hosts
  means nothing (§12); the ratio to the pair-table reference is between
  a counted row and a measured one and is quoted as such.
- `xl-f2` is absent because it does not fit: one pass at degree
  `n_vars`, about `150 s` a call on this shape against Buchberger's
  `7 ms`, and a modelled rather than counted op figure.  It now
  declines above ten unknowns, and the sweep file says why.
- The descent these rows were built with was a truth table, capped at
  `n' = 8` for two summands; every algebraic row above sits at
  `n' = 6`.  §16 replaces the construction with a symbolic one and
  measures past the cap.  Nothing here extrapolates.

### 15.7 Reproducing

```
ic bench --sweep docs/ic/sweeps/solver-engines.json \
  --out docs/ic/runs/ic-bench-solver-engines-2026-09-22.json
ic bench --sweep docs/ic/sweeps/relation-matrix.json \
  --out docs/ic/runs/ic-bench-relation-matrix-2026-09-22.json
```

The frozen reports carry the calibration and which of its units were
pinned, so every GAE in them re-derives from its counts.

## 16. The descent, symbolic: past the truth-table cap

§14 and §15 built their boolean systems from a truth table — every one
of the `2^{m·n'}` points of `V^m` evaluated in the summation
polynomial and Möbius-transformed into algebraic normal forms — which
stopped the algebraic rows at `n' = 8` for two summands and `n' = 5`
for three.  This section replaces the construction with a symbolic
one, checks that nothing changed, and measures where the table could
not go.

### 16.1 The construction, and that it is the same one

The symbolic descent (`pq_descent_symbolic`) expands `S_{m+1}` in the
ring `F_{2^n}[v] / (v² − v)` term by term, with each abscissa the linear
form `x_i = Σ_k v_{i,k} e_k` over the subspace basis.  A polynomial is a
map from a boolean monomial (a 64-bit mask over the `v`'s) to its field
coefficient, and two facts of the ring do all the work: squaring is
linear, `(Σ c_M M)² = Σ c_M² M`, because `M² = M` and every cross term
carries a factor of two; and a product of monomials is their union.
The `S_3` and `S_4` formulas are mirrored operation for operation from
the evaluating code.  `S_3` comes out quadratic in `2n'` variables with
at most `n'² + 2n' + 1` field monomials whatever `n'`; `S_4`, a product
of two cubics in its resultant form, of degree at most six in `3n'`.
The only cap left is the mask: `m·n' ≤ 64`, so `n' ≤ 32` at two
summands and `n' ≤ 21` at three.

Three checks, all in the test suite:

- **Monomial for monomial against the truth table.**  The algebraic
  normal form of a boolean function is unique, so the two
  constructions of the same function must coincide exactly.  They do,
  for `S_3` and `S_4`, on random curves at `n = 7, 9, 11, 13` and every
  `n'` the table can reach.
- **Against the formulas past the table's reach.**  At `n' = 12`
  (24 variables) and at three summands with `n' = 7` (21 variables),
  every equation evaluated at random points of the subspace gives the
  bit of `S_{m+1}` it stands for, and the word formulas agree with the
  `F2mElement` ones.
- **The frozen cells reproduce.**  Re-running §14's `11:6:2` cell for
  both families through the symbolic path gives the frozen operation
  counts to the last digit (`499,453.75` and `438,743.75` monomial
  operations, `D_av = 3.0`, the same enumeration reference), since the
  target draw is unchanged.  A construction that reproduces the frozen
  table exactly has changed nothing but its reach.

### 16.2 The boundaries

As in §15.1: the reference is the pair-table oracle on the *same base*,
same instance, same walk, same two targets; rho was not run on this
instance and no ratio to it is quoted; the degree boundary is the
semi-regular degree of the system's shape — `D_sr = 5` for seventeen
quadratics in eighteen unknowns.  Everything here is a toy,
`E(F_{2^17})`, `r = 65,309`, `log₂ r = 16.0`, and a measurement of what
the engines do when the descent no longer stops them.

### 16.3 One base at `n = 17`, `n' = 9`: eighteen unknowns

Frozen at `docs/ic/runs/ic-bench-solver-engines-n17-2026-09-22.json`;
`random-binary-n17-b1ad26`, `binary-subspace[dimension=9]`: 476 signed
points, 238 columns, walk targets, `incremental-gauss`, two planted
targets, every row a verified logarithm.  `S` and the ratio are per
target; the solver columns are totals over the run.  The 120-second
per-call budget was never reached.

| oracle | engine | trials | rows | hit | S | S / reference | solver calls | solver ops (unit) | ns / op | solver GAE | D_solve | D_sr | wall | correct |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| pair table | — | 122 / 239 | 69 / 134 | 0.566 / 0.561 | 448.1 / 448.9 | 1 | — | — | — | — | — | — | 0.0 s | yes / yes |
| descent-algebraic | exhaustive | 122 / 239 | 69 / 134 | 0.566 / 0.561 | 19,793 / 36,432 | 44.2× / 81.2× | 122 / 239 | 2.83e10 / 5.32e10 (monomial tests) | 0.043 | 5.06e6 / 9.31e6 | — | 5 | 1.2 / 2.3 s | yes / yes |
| descent-algebraic | sat-cdcl | 183 / 222 | 107 / 124 | 0.585 / 0.559 | 4,017,177 / 4,670,150 | 8,966× / 10,404× | 183 / 222 | 6.42e6 / 7.62e6 (conflicts) | 38,853 / 38,023 | 1.03e9 / 1.19e9 | — | 5 | 249 / 290 s | yes / yes |
| descent-algebraic | buchberger-f2 | 122 / 239 | 69 / 134 | 0.566 / 0.561 | 45,377,987 / 87,187,764 | 101,276× / 194,235× | 122 / 239 | 9.96e9 / 1.92e10 (monomial operations) | 282.8 / 281.7 | 1.16e10 / 2.23e10 | 3.99 / 4.00 | 5 | 2,817 / 5,412 s | yes / yes |

The phases behind the reference row, per target: base `30`, pair table
`113,526` (`C(476, 2)` additions — this base is far past the family
optimum for `2^16`, on purpose, so that the descent has eighteen
unknowns to work on), relations `920 / 1,110`, matrix `8 / 25`,
verification `20 / 22`.  The algebraic rows build no table and spend
everything in the solver.

### 16.4 Reading it, against §15

| | `n = 13`, `n' = 6`, 12 unknowns (§15) | `n = 17`, `n' = 9`, 18 unknowns | growth |
|:--|--:|--:|--:|
| pair table, `S` | 40.5 / 39.7 | 448.1 / 448.9 | 11× |
| exhaustive, `S / reference` | 25.9× / 23.9× | 44.2× / 81.2× | |
| sat-cdcl, `S / reference` | 4,659× / 4,139× | 8,966× / 10,404× | |
| buchberger-f2, `S / reference` | 1,920× / 1,768× | 101,276× / 194,235× | |
| buchberger-f2 over exhaustive, in GAE | 75× | 2,293× / 2,394× | 31× |
| buchberger-f2, ns per monomial operation | 51–53 | 282 | 5.4× |
| `D_solve / D_sr` | 2.90 / 4 = 0.73 | 3.99 / 5 = 0.80 | |

- **The degree reading holds and keeps not helping.**  `D_solve = 4.0`
  against `D_sr = 5`: the system is still solved below the semi-regular
  bound, the phenomenon §14 set out to see, and the engine that sees it
  is now `10⁵` times the pair table on the same base.
- **The Gröbner engine's cost grew `2,300×` for six more unknowns**
  while the exhaustive engine's grew `77×` (`5.06e6` from `6.62e4` GAE,
  at almost the same call count), and the price of its own unit grew
  with it: a monomial operation that measured `52 ns` at twelve
  unknowns measures `282 ns` at eighteen, because the polynomials it
  operates on are five times longer.  "Monomial operations" is not a
  stable unit across sizes, which is one more reason the count is
  priced at a measured factor and the factor is recorded.
- **CDCL scales best of the three algebraic engines** (`9,000×` the
  reference from `4,000×`, against exhaustive's `2×` and Buchberger's
  `50×` worsening), and is still four orders of magnitude off the pair
  table.  Its per-conflict cost doubled (`17 µs` to `38 µs`).
- **The reference is not flattered.**  The pair table here is 99 % of
  the reference row's cost and the base is `10×` past the family
  optimum; a base at the optimum would make the reference *cheaper*
  and every ratio above larger.
- **No system lifted to nothing on this instance** save one in the
  first CDCL run; the hit rate is `0.56` on every row, and the CDCL rows
  walked more targets only because a target with several
  decompositions lifted to a different pair.  Same answers.

By §3 every row here is a **measurement**: new rows at a size the
table could not reach, no boundary crossed, no gain claimed.  What
they establish is negative and useful: with the descent no longer the
cap, the shipped engines are, and by a margin that widens with `n'` —
the direction §14.5 warned the reader to expect from "the engine was
the measurement".

### 16.5 The degree ladder past the cap

§14.4 stopped at `n = 13`, `n' = 7`, fourteen unknowns, because the
table did.  The same measurement — eight targets a cell, the Koblitz
curve `K` and a random curve `R` at each degree, `n' = ⌈n/2⌉` so the
system is square, a 120-second budget per target — continues here from
sixteen to twenty-two unknowns.  Frozen at
`docs/ic/runs/ic-descent-degrees-symbolic-2026-09-22.json`; the `K`
family has no instance at `n = 21` in the roster, so that cell is
recorded as skipped rather than filled from elsewhere.

| E | n | n' | m | vars | eqs | D_av | D_pair | D_sr | D_av/D_sr | ops | enumerate | ops/enum (raw counts) | ms | KiB | no decomp | budget hit |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| K | 15 | 8 | 2 | 16 | 15 | 3.8 | 3.8 | 4–5 | 0.94 | 1.976e7 | 3.703e7 | 0.5× | 3,514 | 395 | 2/8 | 0/8 |
| K | 17 | 9 | 2 | 18 | 17 | 4.0 | 4.0 | 5 | 0.80 | 9.345e7 | 2.103e8 | 0.4× | 22,293 | 706 | 3/8 | 0/8 |
| K | 19 | 10 | 2 | 20 | 19 | ≥ 4.0 | ≥ 4.0 | 5 | ≥ 0.80 | ≥ 3.387e8 | 1.226e9 | — | ≥ 114,059 | 1,696 | 6/8 | **6/8** |
| R | 15 | 8 | 2 | 16 | 15 | 3.6 | 3.6 | 4–5 | 0.91 | 1.194e7 | 3.749e7 | 0.3× | 1,705 | 171 | 1/8 | 0/8 |
| R | 17 | 9 | 2 | 18 | 17 | 4.0 | 4.0 | 5 | 0.80 | 9.763e7 | 2.124e8 | 0.5× | 22,704 | 706 | 5/8 | 0/8 |
| R | 19 | 10 | 2 | 20 | 19 | ≥ 4.0 | ≥ 4.0 | 5 | ≥ 0.80 | ≥ 3.640e8 | 1.233e9 | — | ≥ 114,096 | 1,657 | 3/8 | **7/8** |
| R | 21 | 11 | 2 | 22 | 21 | ≥ 4.0 | ≥ 4.0 | 5 | ≥ 0.80 | ≥ 1.053e9 | 6.671e9 | — | ≥ 120,003 | 3,620 | 4/8 | **8/8** |

A row with a budget hit carries lower bounds in every measured column
— a run that was stopped had not finished raising its degree or
spending its operations — and its `ops/enum` is left blank, since a
lower bound over a fixed reference says nothing.  `D_sr` reads `4–5`
where the bound differs between targets of one cell: a target whose
descended system carries a linear coordinate equation has a lower
bound than one whose equations are all quadratic.

Reading it:

- **The solving degree is `4` from eighteen unknowns on, against a
  bound of `5`**: every one of the thirty-two finished runs at
  `n' ≥ 9` reached exactly degree four, on both families, and so did
  every run the budget stopped.  §14's phenomenon — the descended
  system solved below the semi-regular degree — holds at every cell
  the symbolic descent reached, at `0.80` of the bound; at sixteen
  unknowns, where the bound is `4` or `5` by target, it holds at
  `0.91`–`0.94`.
- **`K` and `R` still behave identically**, on the degree and on the
  cost, as they did in §14.4.
- **The engine is the wall.**  Per target, the finished runs cost
  `1.7`–`3.5 s` at sixteen unknowns and `22 s` at eighteen (`8`–`51 s`
  by target); at twenty, thirteen of sixteen runs hit the 120-second
  budget, and at twenty-two all eight did.  That is a factor of six to
  thirteen per two unknowns, the growth §16.4 saw from the inside of a
  whole run, and it is the reason the algebraic rows stop where they
  do: not the descent any more, and not the degree, which is flat.
- **The `ops/enum` column is in raw counts of unlike units** and is
  kept only to line up with §14.4; §15.4 is the correction, and the
  converted ratio (about `2,300×` at eighteen unknowns, §16.4) is the
  one that means something.

By §3 these rows are **measurements**: the degree ladder extended by
three rungs, no boundary crossed, the timed-out cells marked as the
lower bounds they are.

### 16.6 What does not count

- Nothing here is a speed against rho; the `vs rho` column is empty.
- Every algebraic price is `measured` and host-dependent (§12); the
  ratio to the pair-table reference is between a counted row and a
  measured one, and is quoted as such.  The frozen file carries its
  calibration.
- The pair-table reference is on a base `10×` past the family optimum
  and is priced accordingly; it is the reference for *these* algebraic
  rows on *this* base, not the best index-calculus row at `2^16`.
- The symbolic descent lifts the descent's cap and nothing else: the
  engines' caps (26 unknowns for the enumerating extractions, 10 for
  XL, the budget for Buchberger) are where the rows now stop, and §8 of
  the framework manual says so.

### 16.7 Reproducing

```
ic bench --sweep docs/ic/sweeps/solver-engines-n17.json \
  --out docs/ic/runs/ic-bench-solver-engines-n17-2026-09-22.json

ic descent --cells 15:8:2,17:9:2,19:10:2,21:11:2 --targets 8 \
  --families K,R --budget-seconds 120 \
  --out docs/ic/runs/ic-descent-degrees-symbolic-2026-09-22.json

# The reproduce check: the frozen 11:6:2 cell through the symbolic path.
ic descent --cells 11:6:2 --targets 8 --families K,R
```

## 17. F4 and F5 behind the plug point

§16 ended on "the engine is the wall": the only Gröbner engine behind
the `SystemSolver` plug point was a Buchberger that takes pairs one at a
time, and it cost six to thirteen times more per two extra unknowns
until it hit the budget at twenty.  This round puts behind that plug
point the engines the thread did not have:

| engine | what it is | module |
|:--|:--|:--|
| `f4-f2` | Faugère's F4 over `F_2[v]/(v² − v)`: normal strategy, Gebauer–Möller criteria, the field products `v·g` as pairs, symbolic preprocessing, bit-packed elimination; a full reduced basis, solutions read off its linear elements | `pq_f4_f2` (new) |
| `matrix-f4` | the Koblitz oracle's hybrid: Macaulay matrices to a fixed degree, propagation, splitting | `koblitz_groebner` |
| `matrix-f5` | the same, with the rows the Boolean F5 criterion predicts to reduce to zero left out | `koblitz_groebner` |
| `inherited-f4` | the same, children specialising their parent's reduced basis | `koblitz_groebner` |
| `crossbred-f2` | Joux–Vitse: a Macaulay left kernel at degree `D`, then `2^k` bit-sliced linear solves | `crossbred` |
| `fes-f2` | fast exhaustive search, libfes-lite's Gray code: two word XORs per point | `mq_fes` |

and asks whether any of them moves an algebraic row.  §17.1 and §17.2
were written and committed before the suite ran; §17.3–§17.5 are checks
declared after reading parts of it, each committed before it ran; §17.6
onward reads the results.

### 17.1 The boundaries, stated before measuring

**At the stage** — one decomposition oracle call on one target, a stage
diagnostic and never a speed:

- **The reference is exhaustive search, not Buchberger.**  `fes-f2`
  wherever the system is quadratic, which is every two-summand descent
  (`S₃` descends to quadratics, §16.1), and the `exhaustive` evaluator
  where it is not (three summands, degree up to six).  It is the best
  algorithm that already decides these systems at these sizes, and a
  Gröbner engine that does not beat it on a cell has not earned its
  place in that cell.
- **The baseline is `buchberger-f2`**, the engine every frozen algebraic
  row so far used (§14–§16).  Beating it is engineering; beating the
  reference is the question.
- **The degree floor is `D_sr`**, the semi-regular degree of the
  system's shape (§14.3), derived per target: `5` from eighteen
  unknowns on for two summands (§16.5), `7` and `8` at nine and twelve
  unknowns for three.  The column it bounds is `D_learn`, the highest
  step degree at which an engine produced a new basis element.

**For the whole method** — `S`, every phase priced:

- **The method reference is counted Pollard rho on the same subgroup.**
  `ic bench` now runs it (`--rho-runs 16`: sixteen counted runs over the
  planted targets, with its spread), so the `vs rho` column §15 and §16
  left empty is filled from this round on.  At these toy orders rho's
  `S` is setup-dominated — a development run put it at `16.2` (range
  `10.1`–`33.1`) at `r = 4,091` — and the column is the ratio to the rho
  measured on the same instance, never to the asymptotic `1.3`.
- **The oracle reference is the pair table on the same base** (§15.1),
  and the **floor** is the generic one, `√(π/2A)` in `S` (§3).

**A known caveat on the baseline's degree, stated here because it
changes how §14–§16's degree column reads.**  The boolean Buchberger
processes S-pairs only, never the field products `v_i·g` for `v_i` in a
leading monomial.  Its output is therefore not guaranteed to be a
Gröbner basis of the boolean ideal: on the toy system
`{v₀v₁ + v₀ + v₁}` it stops with three standard monomials for one
solution, where `v₀·g = v₀` closes it.  A development check certified
its outputs as boolean bases on all thirty-two descent targets probed
(`n = 7`–`13`: every S-pair and every field product reduces to zero,
and the standard monomials number the solutions); the whole-pipeline
rows of §15 and §16.3 were verified end to end, logarithm and all; and
the suite below compares its verdict with the reference's on every
target, which is the check for the cells in between.  Its *degree* is
another matter: an engine that must
rediscover what a field product gives for free can need a higher
degree to do it, so a `D_solve` from the pair-only engine is an upper
bound on what the ideal needs, and the F4 engine, which processes field
pairs at sugar degree `deg LM + 1`, is the one whose `D_learn` is
compared with `D_sr` from here on.

### 17.2 The suite, and the target, declared before the suite ran

**The protocol.**  The WDSat regression suite of `AGENTS.md` §8 is a
SAT-counter protocol and says of itself that it is not an adapter for
other engines, so this round freezes an equivalent matched suite under
the parent accounting contract, at
`research/ic_framework_engines_20260922/`:

- **Stage cells**, `ic descent --solver …`, the descent of §14/§16 with
  the same target draw: two summands at `n:n' = 7:4, 9:5, 11:6, 13:7,
  15:8, 17:9, 19:10, 21:11` (eight to twenty-two unknowns, square
  systems) and three summands at `7:3, 9:3, 11:4`; families `K` and `R`
  (`K` has no instance at `n = 21`); eight targets a cell; seed
  `0x0DE5CE47`; three repetitions, interleaved per target with the
  engine order rotated each repetition; a 120-second budget per call.
  Engines: `buchberger-f2` (baseline), `f4-f2`, `matrix-f4`,
  `matrix-f5`, `inherited-f4`, `crossbred-f2` at its module defaults
  (`D = 3`, `k = 8`), `fes-f2`, `exhaustive`, and `sat-cdcl` on the
  separate first-solution leaderboard.
- **Past Buchberger's reach**, `23:12` and `25:13` (twenty-four and
  twenty-six unknowns), the same protocol without `buchberger-f2`, whose
  every call at twenty-two unknowns already hit the budget (§16.5), and
  without `sat-cdcl`, which is on a separate leaderboard and was
  `5,000×` the reference at twenty-two in development.
- **A holdout**: seed `0x0F4F5EED` — new targets, and for `R` a new
  curve — at `15:8, 17:9, 21:11, 23:12` and `9:3:3`, same engines.
- **The whole method**, `ic bench`: the §15 base at `n = 13` and a
  square base at `n = 15` (`dimension = 8`), each on three curves
  (seeds `20260922`, `20260923`, `20260924`), and the §16.3 base at
  `n = 17` on its curve; two planted targets each; the pair table, then
  the descent oracle once per engine including `buchberger-f2`; rho on
  every instance.
- **Checks**: every system fingerprinted with blake3 and every cell's
  reference answers digested, so a later run proves it saw the same
  inputs; every engine's answers compared with the reference's on every
  target and repetition; a call that exhausts its budget is recorded,
  never repeated on that target and never read as "no solution";
  every whole-pipeline logarithm verified against the planted one.
- **Not in it**: `xl-f2`, for §15.6's reasons (one pass at degree
  `n_vars`, no budget hook, a modelled count; it declines above ten
  unknowns).

What was known when this was written, disclosed so that the target can
be judged against it: a one-target development probe per cell, run to
size the budget and to check that no engine hangs, put `f4-f2` at
`0.76 s` a call at twenty-two unknowns (Buchberger: over the 120-second
budget on every §16.5 target), and `crossbred-f2` at `4.2×` and `1.8×`
`fes-f2`'s wall time at twenty and twenty-two unknowns on one target
each.  None of the probe's numbers is used below.

**Success, all three:**

1. **Reach.**  An F4-family engine (`f4-f2`, `matrix-f4`, `matrix-f5`,
   `inherited-f4`) decides every two-summand target at sixteen to
   twenty-two unknowns, main and holdout, within the budget, with zero
   disagreements against the reference.
2. **Engineering gate.**  The contract's promotion gate, on the whole
   method: at `n = 13` and `n = 15` on all three curves, the best
   F4-family row costs at least `20 %` less than the `buchberger-f2`
   row per verified relation, all phases priced, with a `95 %` paired
   bootstrap interval that excludes no improvement.
3. **The question.**  Some algebraic engine's per-call wall time falls
   below the reference's on a two-summand cell of the stage suite: the
   median over the cell's targets of `engine / fes-f2` (each target's
   median over repetitions) below one, with its `95 %` bootstrap
   interval below one too.  A ratio that falls from sixteen to
   twenty-six unknowns without crossing is recorded, and the crossover
   it implies is extrapolated and marked as extrapolation; it does not
   meet this target.

**Abandon** the algebraic engines for two-summand descents at these
sizes if no algebraic engine's ratio to `fes-f2` at twenty-six unknowns
is below its ratio at sixteen: the engines are then not gaining on
enumeration over the range this repository can reach, and more
engineering on them is engineering on the wrong side of the boundary.

**Classes, fixed in advance.**  A lower wall time than the baseline's
is **engineering**, however large.  A `D_learn` below the baseline's
`D_solve` is **accounting** on the degree column, not an advance: the
ideal did not change, only whether the engine closed it (§17.1).  An
**advance** at the stage is target 3's first clause; for the method it
would be a whole-pipeline `S` below the same instance's pair-table
row, which nothing in development came near.

**Inadmissible**: changing a cell, a seed, a target count, the budget or
an engine's parameters after reading the suite's numbers; dropping a
budget-exhausted or declined call from a row; pricing an
elimination-only count at the word-XOR ratio (every candidate engine's
count is qualified and priced by measured wall time, §15.4); ranking
the first-solution engine against the complete enumerators; and
tuning `crossbred-f2`'s `D` and `k` per cell in this round — a tuned
crossbred is a separate round with its own declaration.

### 17.3 A confirmation, declared after the stage parts were read

The frozen stage parts met target 3 on exactly one cell:
`crossbred-f2` below `fes-f2` at twenty-six unknowns (`25:13:2`, the
`R` curve — `K` has no instance at `n = 25`), a median ratio of `0.77`
over eight targets with a `95 %` interval of `[0.72, 0.79]`; at
twenty-four unknowns, over thirty-two targets on both seeds and both
families, the same ratio is `1.04 [1.03, 1.07]`.  One curve and one
seed at the one size that crosses is thin, so before anything is read
into it this confirmation was written down, committed, and then run
(`confirm_20260923.json`, results in `results/confirm_v1/`):

- the holdout seed `0x0F4F5EED` at `25:13:2` — a new `R` curve and new
  targets;
- one size further, `27:14:2` (twenty-eight unknowns), on both seeds,
  and `29:15:2` (thirty) on the main seed;
- the `B` part's engines, and the suite's protocol unchanged.

**It confirms** if `crossbred-f2 / fes-f2` is below one with its
interval below one on the holdout at twenty-six unknowns and on each
seed at twenty-eight.  A miss on any of the three and the crossing is
reported as seed-dependent.  Thirty unknowns is recorded, not graded.
This was declared after the numbers it checks were read, and says so:
it can weaken the finding or support it, not create it.

### 17.4 Is the reference strong enough?  Declared before the vector search ran

The crossing of §17.3 is against this repository's fast exhaustive
search: a faithful port of libfes-lite's scalar kernel — every equation
one bit of a 64-bit word, sixteen Gray-code steps unrolled per focus
advance, two word XORs a point — measured at `0.58 ns` a point on this
host, about 1.2 cycles at its 2.1 GHz base clock.  That is close to what
a scalar walk can do, and it is not the best exhaustive search there
is.  libfes-lite's vector kernels fix a few variables differently in
each lane of a SIMD register and walk every lane's sub-cube with one
Gray code: the quadratic part of the free variables is the same in
every lane, so one broadcast second derivative and one vector XOR
advance them all.  With at most 32 equations a 32-bit lane holds a
system, so AVX-512 walks sixteen sub-cubes a step.  §1 asks for the best
algorithm that already solves the problem, and a crossing against a
reference that a known technique makes several times faster may be a
crossing of the implementation rather than of the algorithm.  So the
reference is strengthened before the crossing is read as anything:

- **`fes-f2-wide`** (new): libfes-lite's unrolled kernel with the
  first-derivative table held as vectors of 32-bit lanes — sixteen on
  AVX-512, eight on AVX2 — and the second derivatives broadcast;
  `2^{n−4}` vector steps for `2^n` points.  It is checked against the
  scalar search on random systems with planted roots from the smallest
  size the unrolled chunk takes to 32 equations, and it joins the
  registry's agreement tests with every other engine.
- **The check** (`strength_20260923.json`, results in
  `results/strength_v1/`): `crossbred-f2`, `fes-f2`, `fes-f2-wide` and
  `inherited-f4` on `23:12:2`, `25:13:2`, `27:14:2` and `29:15:2`
  (twenty-four to thirty unknowns), both seeds, both families where an
  instance exists, the suite's protocol otherwise.  `fes-f2-wide` is the
  reference wherever it applies.

**The crossing stands** if `crossbred-f2 / fes-f2-wide` is below one
with its interval below one at some size on both seeds.  **If it does
not**, §17.3's crossing is recorded as a crossing of the scalar
reference only, and the size at which the strengthened reference would
be crossed is extrapolated from the measured ratios and marked as
extrapolation.  Either way the vector search replaces the scalar one as
the stage reference from here on, and the frozen suite's rows against
`fes-f2` stay as they are, labelled with the reference they were read
against.

### 17.5 The hybrids past thirty unknowns, declared after §17.4 was read

§17.4's check settled the question it was declared for (§17.9 reports
it) and showed something it was not looking for: against the vector
reference, `inherited-f4` falls from `14.8×` at twenty-eight unknowns
to `4.55×` at thirty, on both seeds and both families, because its
cost per call barely grows there (`0.67 s` to `0.77 s`) while the
reference's quadruples.  A ratio that falls three-fold in one step is a
trend or an accident of two sizes, and the only way to tell is to
measure further out.  So, declared after reading §17.4 and before
running:

- `31:16:2` and `33:17:2` (thirty-two and thirty-four unknowns), both
  seeds, both families where an instance exists
  (`hybrids_20260923.json`, results in `results/hybrids_v1/`);
- engines `inherited-f4`, `matrix-f4`, `matrix-f5`, `crossbred-f2` (to
  record where its defaults stand), `fes-f2` (which stops at
  thirty-two) and `fes-f2-wide`, whose lanes reach thirty-six and which
  is the reference; the suite's protocol otherwise.  At thirty-four
  unknowns no scalar search checks the reference, so there the
  hybrids' agreement with it, by an unrelated algorithm, is the check.

**A stage crossing against the strengthened reference** is some hybrid
below `fes-f2-wide` with its interval below one at thirty-two or
thirty-four unknowns on both seeds.  If none is, the trajectory is
recorded, and any crossing it implies is extrapolated from the six
sizes 24–34 and marked as extrapolation.  Either way it is a stage
diagnostic: no whole-method `S` exists at these sizes.

**Amendment, before the thirty-four-unknown cell was rerun.**  The
declaration above says the vector search's lanes reach thirty-six
unknowns; they do, but a lane holds 32 equations and `33:17:2` has 33,
so on the first run (`results/hybrids_v2/`) the search declined that
cell and it had no exhaustive reference at all — the comparison fell
back to `inherited-f4` checking itself.  The search now walks the
first 32 equations and filters its candidates with the rest, as libfes
does for more equations than a word (about `2^{n−32}` spurious
candidates a call), tested against the scalar search on systems of 33
to 55 equations.  The cell is rerun on both seeds, unchanged otherwise
(`hybrids_34_20260923.json`, `results/hybrids_v3/`); the first run's
rows are kept and read only at thirty-two unknowns.

### 17.6 Reach: target 1, met

**The answer first.**  F4 and F5 move the algebraic rows, and move them
a long way: the F4-family engines decide the two-summand targets the
pair-only Buchberger could not reach past eighteen unknowns — `f4-f2`
through twenty-eight, the matrix hybrids through thirty-four — and take
the Gröbner row's whole-method `S` down by `12×` to `1,535×` — engineering,
with no boundary crossed, since the best algebraic row is still at
least `24×` the pair table on its base and `333×` rho.  The one stage
crossing the suite found, crossbred below fast exhaustive search at
twenty-six and twenty-eight unknowns, was against a scalar reference and
did not survive its vector form.  What survives is a lead: the inherited
hybrid gaining on the vector reference at half a ratio per two unknowns,
`1.20×` at thirty-four, with a crossing extrapolated near thirty-five.

Frozen at `research/ic_framework_engines_20260922/results/baseline_v1/`
(parts `A1`–`A3`, `B`, `H1`, `H2`) and `results/confirm_v1/`; every
system is fingerprinted in the reports and `manifest.json`.  Two-summand
cells, both families and both seeds pooled where they exist: targets
decided within the 120-second budget in all three repetitions, and the
median over targets of each target's median wall time.

| unknowns | `D_sr` | `buchberger-f2` | `f4-f2` | `matrix-f4` | `matrix-f5` | `inherited-f4` | `crossbred-f2` | `sat-cdcl` (first) | `fes-f2` | `exhaustive` |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 12 | 4 | 16/16 · 36 ms | 16/16 · 2.1 ms | 16/16 · 1.6 ms | 16/16 · 1.6 ms | 16/16 · 1.4 ms | 16/16 · 0.98 ms | 16/16 · 3.9 ms | 16/16 · 0.049 ms | 16/16 · 0.097 ms |
| 14 | 4 | 16/16 · 372 ms | 16/16 · 4.5 ms | 16/16 · 5.6 ms | 16/16 · 5.8 ms | 16/16 · 3.0 ms | 16/16 · 1.5 ms | 16/16 · 68 ms | 16/16 · 0.085 ms | 16/16 · 0.49 ms |
| 16 | 4–5 | 32/32 · 2.1 s | 32/32 · 9.9 ms | 32/32 · 11.7 ms | 32/32 · 11.9 ms | 32/32 · 5.4 ms | 32/32 · 2.1 ms | 32/32 · 160 ms | 32/32 · 0.15 ms | 32/32 · 2.2 ms |
| 18 | 5 | 32/32 · 21 s | 32/32 · 21 ms | 32/32 · 30 ms | 32/32 · 30 ms | 32/32 · 10 ms | 32/32 · 4.1 ms | 32/32 · 1.2 s | 32/32 · 0.28 ms | 32/32 · 9.2 ms |
| 20 | 5 | **2/16** · 80 s | 16/16 · 39 ms | 16/16 · 83 ms | 16/16 · 83 ms | 16/16 · 29 ms | 16/16 · 4.2 ms | 14/16 · 6.7 s | 16/16 · 0.78 ms | 16/16 · 45 ms |
| 22 | 5 | **0/16** | 16/16 · 694 ms | 16/16 · 210 ms | 16/16 · 207 ms | 16/16 · 38 ms | 16/16 · 5.0 ms | 6/16 · 6.6 s | 16/16 · 2.6 ms | 16/16 · 197 ms |
| 24 | 6 | — | 32/32 · 2.8 s | 32/32 · 130 ms | 32/32 · 137 ms | 32/32 · 77 ms | 32/32 · 10 ms | — | 32/32 · 9.7 ms | 32/32 · 947 ms |
| 26 | 6 | — | 16/16 · 7.7 s | 16/16 · 295 ms | 16/16 · 311 ms | 16/16 · 178 ms | 16/16 · 29 ms | — | 16/16 · 39 ms | 16/16 · 4.5 s |
| 28 | 6 | — | 16/16 · 17.5 s | 16/16 · 710 ms | 16/16 · 738 ms | 16/16 · 620 ms | 16/16 · 105 ms | — | 16/16 · 154 ms | declines |
| 30 | 6 | — | **5/16** · 30 s | 16/16 · 1.7 s | 16/16 · 1.8 s | 16/16 · 760 ms | **0/16** | — | 16/16 · 612 ms | declines |

Every decided answer on every target and repetition equals the
reference's; the suite's `disagreements` list is empty in every part.

- **Target 1 is met by all four F4-family engines**, on both seeds:
  every two-summand target from sixteen to twenty-two unknowns decided,
  inside the budget, in every repetition.  Where the pair-only
  Buchberger decided 2 of 16 targets at twenty unknowns and none at
  twenty-two, `f4-f2` takes 39 ms and 0.69 s, and the hybrids 29–210 ms.
- **The wall moved ten unknowns for F4 and sixteen for the hybrids.**
  `f4-f2` decides every target through twenty-eight unknowns (17.5 s a
  call) and stops at thirty, where eleven of sixteen systems need a
  degree-five matrix past its 1 GiB cap; the matrix hybrids decide
  everything through thirty-four and `inherited-f4` through thirty-two,
  leaving two of sixteen at thirty-four (§17.9).  At the sizes both
  decide, the paired per-call ratio of
  `f4-f2` to `buchberger-f2` is `0.04–0.05` at twelve unknowns, `0.004–
  0.006` at sixteen and `0.0007–0.0014` at eighteen: 20× to 1,400×,
  growing with the size, on every cell's interval.
- `sat-cdcl`, on its own first-solution leaderboard, runs out of its
  conflict budget on 2 of 16 targets at twenty unknowns and 10 of 16 at
  twenty-two.  `crossbred-f2` at its module defaults (`D = 3`, `k = 8`)
  decides everything through twenty-eight and nothing at thirty, where
  those parameters stop fitting the system: a budget verdict on every
  target, never a wrong one.

### 17.7 The degree column, re-read: accounting

The same runs, the engine's own degree statistics (mean over the decided
targets; `D_learn` the highest step degree that produced a new element,
`D_reach` the highest processed):

| unknowns | `D_sr` | `f4-f2` `D_learn` | `f4-f2` `D_reach` | `buchberger-f2` `D_solve` |
|--:|:--|--:|--:|--:|
| 12 | 4 | 3.00 | 3.69 | 3.00 |
| 14 | 4 | 3.00 | 3.56 | 3.25 |
| 16 | 4–5 | 3.00 | 3.81 | 3.84 |
| 18 | 5 | 3.00 | 3.56 | 3.97 |
| 20 | 5 | 3.00 | 3.44 | 4.00 (2 targets) |
| 22 | 5 | 3.50 | 4.25 | — |
| 24 | 6 | 3.66 | 4.66 | — |
| 26 | 6 | 3.62 | 4.62 | — |
| 28 | 6 | 3.94 | 4.94 | — |

§16.5 read the descended systems as solved at `0.80` of the semi-regular
bound from eighteen unknowns on, with the pair-only engine's `D_solve =
4`.  An engine that closes the boolean ideal — the field products as
pairs, §17.1 — learns everything it needs at degree **three** through
twenty unknowns: `0.60` of the bound.  The ideal did not change; the
engine that measured it did, and the extra degree was the pair-only
engine rediscovering what a field product gives for free.  By §3 that
is **accounting** on the degree column: §14–§16's degree readings stand
as statements about that engine, and the statement about the systems is
now `0.60`, with `D_learn` rising to `3.5`–`3.9` (`0.60`–`0.70` of a
bound that itself rises to six) from twenty-two unknowns on.  Petit and
Quisquater's phenomenon is stronger here than §16 said, and — as §14.5
warned — it is still a statement about degrees, not costs: `f4-f2` is
the engine that measures it and the one whose cost grows fastest past
twenty unknowns, because its processed degree reaches five.

### 17.8 The stage against the scalar reference: target 3, met, and confirmed

Wall time over `fes-f2`, per target (each target's median over
repetitions), median over the targets with a 95 % percentile bootstrap
interval, both families and both seeds pooled:

| engine | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| `buchberger-f2` | 14,834 | 71,649 | — | — | — | — | — | — |
| `f4-f2` | 70.2 [60, 74] | 71.1 [56, 81] | 47.7 [45, 78] | 260 [43, 282] | 288 [284, 293] | 180 / 206 | 114 [110, 121] | 49.2 (5 of 16 decided) |
| `matrix-f4` | 84.5 | 105 | 103 | 81.0 | 13.4 | 7.76 / 7.78 | 4.56 | 2.75 |
| `matrix-f5` | 82.3 | 105 | 105 | 79.8 | 14.2 | 8.01 / 8.09 | 4.81 | 2.87 |
| `inherited-f4` | 37.3 | 35.7 | 38.4 | 14.4 | 8.07 | 4.83 / 4.50 | 4.10 | 1.24 |
| `crossbred-f2` | 14.5 [14.2, 14.8] | 14.5 [14.1, 14.8] | 5.33 [5.0, 5.6] | 1.91 [1.8, 2.0] | 1.05 [1.03, 1.07] | **0.77** [0.72, 0.79] / **0.735** [0.69, 0.77] | **0.686** [0.66, 0.74] / **0.67** [0.65, 0.70] | fails |

(Sixteen to twenty-four unknowns pool the frozen suite's main and holdout seeds; `26` is the frozen `B` part / §17.3's holdout; `28` pools §17.3's two seeds, crossbred shown per seed; `30` is §17.3's main seed.)

- **Target 3 is met, against this reference.**  `crossbred-f2` falls
  below `fes-f2` at twenty-six unknowns (`0.77 [0.72, 0.79]`, the `R`
  curve; `K` has no instance at `n = 25`), and §17.3's confirmation,
  declared after that number was read, holds on all three of its cells:
  `0.735` on the holdout at twenty-six, `0.686` and `0.67` on the two
  seeds at twenty-eight, every interval below one.  At twenty-four
  unknowns it is `1.05 [1.03, 1.07]` over thirty-two targets.
- **The abandon condition is not triggered.**  `matrix-f4`, `matrix-f5`,
  `inherited-f4` and `crossbred-f2` all have a lower ratio at twenty-six
  unknowns than at sixteen; `f4-f2` does not (`70 → 180`), because its
  processed degree reaches five from twenty-two unknowns on.
- **The crossing is a window, not an asymptote.**  At thirty unknowns
  `crossbred-f2`'s default parameters stop fitting, and a larger `k` or
  `D` is a tuned crossbred, which §17.2 put out of this round.

### 17.9 Against the vector reference: the crossing was the implementation's, and the hybrids close in

Frozen at `results/strength_v1/` (§17.4), `results/hybrids_v2/` (§17.5,
read at thirty-two unknowns only) and `results/hybrids_v3/` (§17.5's
amended thirty-four-unknown cell).  `fes-f2-wide` answers exactly as
the scalar search on every target both run, and is `3.0`–`4.0×` faster:
`3.4`–`3.5 ms` against `10.5`–`10.7 ms` at twenty-four unknowns, `0.65`–
`0.68 s` against `2.4`–`2.6 s` at thirty-two, `2.6`–`2.7 s` at thirty-four
where the scalar walk stops.  Sixteen lanes buy four, not sixteen: the
walk is bound by its dependent table updates and the zero test, not by
the XORs.  Ratios to it, as in §17.8:

| engine | 24 | 26 | 28 | 30 | 32 | 34 |
|:--|--:|--:|--:|--:|--:|--:|
| `crossbred-f2` | 3.00 [2.81, 3.10] | 2.54 [2.46, 2.59] | 2.45 [2.36, 2.55] | defaults fail | defaults fail | defaults fail |
| `inherited-f4` | 22.6 [21.7, 23.2] | 15.4 [14.7, 15.8] | 14.8 [14.5, 15.8] | 4.55 [4.47, 4.62] | 2.09 [2.03, 2.11] | 1.20 [1.20, 1.23], 14 of 16 decided |
| `matrix-f4` / `matrix-f5` | | | | | 5.40 / 5.42 | 2.45 / 2.48 |
| `fes-f2`, the scalar walk | 3.06 | 3.60 | 3.86 | 3.95 | 3.6–3.9 | declines |

Every cell pools both seeds and both families where an instance exists
(`K` has none at `n = 25`, `27` or `33`); no engine disagreed with the
reference on any target.

- **§17.4's check: the crossing does not stand.**  `crossbred-f2` is
  `2.45`–`3.0×` the vector form of the reference wherever its default
  parameters apply, on both seeds.  As declared, §17.3's crossing is
  recorded as a crossing of the scalar implementation only.  There is
  nothing to extrapolate at these defaults: the ratio flattens (`3.00 →
  2.54 → 2.45`) and the parameters stop fitting at thirty unknowns.
- **§17.5's check: no hybrid crosses at thirty-two or thirty-four.**
  `inherited-f4` is `2.09×` at thirty-two and `1.20×` at thirty-four
  on the targets it decides; at thirty-four it leaves two of the
  holdout's eight undecided (its node budget), so it is not below the
  reference on that cell whatever the ratio on the rest.  The matrix
  hybrids decide everything and are `2.45`–`2.48×`.
- **Against the scalar walk it would have been a second crossing.**  At
  thirty-two unknowns `inherited-f4` takes `0.55×` the scalar search's
  time (every one of the thirty-two targets between `0.48` and `0.62`),
  a cleaner crossing than crossbred's; against the vector form it is
  `2.09×`.  Two engines crossed the scalar reference in this round and
  neither crosses the strengthened one, which is the reason §17.4
  strengthened it before reading either.
- **The trajectory, extrapolated and marked as such.**  `inherited-f4`'s
  ratio to the vector reference falls `×0.54` per two unknowns on a
  log-linear fit over the six sizes 24–34, which puts the crossing at
  `34.9` unknowns.  The fit is poor in shape (flat from twenty-six to
  twenty-eight, residual `+0.55` there, steep after), the engine starts
  leaving targets undecided at thirty-four, and thirty-six unknowns is
  where both the vector reference (thirty-two walked, four in the lanes)
  and a practical curve constructor (an exhaustive `2^{35}` point count)
  stop.  It is the most specific lead this round leaves — a hybrid F4
  gaining on the strongest exhaustive search at about half a ratio per
  two unknowns — and it is not a result.

### 17.10 Three summands

Reference `exhaustive` (the `S_4` descents are of degree six, so no
fast exhaustive search applies); nine unknowns at `7:3:3` and `9:3:3`,
twelve at `11:4:3`, both families, and the holdout at `9:3:3`:

| unknowns | `D_sr` | `exhaustive` | `matrix-f4` | `matrix-f5` | `inherited-f4` | `f4-f2` (`D_learn`) | `sat-cdcl` | `buchberger-f2` |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| 9 | 7 | 0.04 ms | 2.0–2.2 ms (44–62×) | 2.2–2.4 ms (47–65×) | 8.5–9.6 ms | 7–18 ms (7.0) | 10–23 ms | 20 s |
| 12 | 8 | 1.2–1.3 ms | 14.8–14.9 ms (12×) | 14.8–14.9 ms (13×) | 116–118 ms | 244–416 ms (7.0) | 0.6–0.75 s | over budget, 16/16 |

Nothing algebraic comes near enumeration at these sizes; the matrix
hybrid's ratio falls from `44`–`62×` to `12×` between nine and twelve
unknowns, two sizes and no fit.  `F4` learns at degree seven on both
shapes: at the semi-regular bound on nine unknowns, one below it on
twelve.  `crossbred-f2`'s defaults do not apply to degree-six systems.

### 17.11 The whole method: target 2, met; `S` against rho

Frozen at `results/baseline_v1/W13-*`, `W15-*` and `W17-*`: complete,
verified logarithms, every phase inside `S`, two planted targets per
curve, counted rho on every instance.  Medians over curves and targets.

| row | `n = 13`, 12 unknowns: `S` · / pair table · / rho | `n = 15`, 16 unknowns | `n = 17`, 18 unknowns |
|:--|--:|--:|--:|
| counted rho, `S` (per curve) | 16.2 · 18.7 · 14.6 | 13.1 · 9.2 · 8.9 | 5.2 |
| pair table (the oracle reference) | 51.1 · 1 · 3.5× | 267 · 1 · 29× | 448 · 1 · 86× |
| descent + `fes-f2` | 370 · 7.2× · 22× | 367 · 1.37× · 40× | 763 · 1.70× · 146× |
| descent + `exhaustive` | 827 · 16× · 50× | 6,512 · 24× · 710× | 28,199 · 63× · 5,413× |
| descent + `crossbred-f2` | 8,514 · 167× · 509× | 6,597 · 25× · 722× | 10,929 · 24× · 2,098× |
| descent + `inherited-f4` | 5,564 · 109× · 344× | 17,861 · 67× · 1,949× | 28,719 · 64× · 5,513× |
| descent + `f4-f2` | 5,211 · 102× · 333× | 27,143 · 102× · 2,961× | 42,647 · 95× · 8,186× |
| descent + `matrix-f4` / `matrix-f5` | 6,440 / 6,663 · 126–130× | 36,406 / 37,356 · 136–140× | 77,482 / 77,954 · 173–174× |
| descent + `sat-cdcl` (first solution) | 149,661 · 2,927× | 1.11×10⁶ · 4,157× | 4.45×10⁶ · 9,932× |
| descent + `buchberger-f2` (the baseline) | 63,622 · 1,244× · 3,832× | 7.74×10⁶ · 29,029× · 844,785× | 6.55×10⁷ · 145,974× · 1.26×10⁷× |

`n = 13` is §15's base (`dimension = 6`, thirteen equations in twelve
unknowns) on three curves, `n = 15` a square base (`dimension = 8`)
on three curves, `n = 17` §16.3's base and curve; two planted targets
each, every logarithm verified, 140 rows.  The fresh `buchberger-f2`
rows at `n = 17` reproduce §16.3's frozen ones to 3 % (`46.2` and
`84.7` million against `45.4` and `87.2`), on a host that was running
three other parts beside them.

**Target 2 — the contract's gate — is met** at `n = 13` and `n = 15`,
by every F4-family engine, on all three curves: total cost per verified
relation over `buchberger-f2`'s, median over the six (curve, target)
pairs with a 95 % bootstrap interval,

| | `f4-f2` | `inherited-f4` | `matrix-f4` | `matrix-f5` |
|:--|--:|--:|--:|--:|
| `n = 13` | 0.088 [0.074, 0.102] | 0.091 [0.082, 0.108] | 0.132 [0.107, 0.155] | 0.115 [0.106, 0.157] |
| `n = 15` | 0.0035 [0.0034, 0.0036] | 0.0022 [0.0021, 0.0024] | 0.0045 [0.0044, 0.0047] | 0.0046 [0.0045, 0.0050] |
| `n = 17` (one curve, two targets; not graded) | 0.00065 | 0.00043 | 0.0012 | 0.0012 |

against a gate of `0.8`.  The Gröbner row's `S` falls `12×` at
`n = 13`, `285×` at `n = 15` and `1,535×` at `n = 17` (medians),
the factor growing with the size as the per-call ratios of §17.6
do.

Reading it:

- **Nothing is below rho, and nothing moved toward it that a generic
  algorithm could not.**  The best row on every instance is still the
  pair table, `3.5×`, `29×` and `86×` rho on bases chosen past the
  family optimum so the descent has unknowns to work on (§16.3); every
  algebraic row is at least `333×` rho.  (Against the matched rho of
  §18: `21.0×`, `133×` and `257×`, and at least `2,024×`.)  F4 took the Gröbner row from
  `1,244×` the pair table to `102×` at `n = 13` and from `145,974×`
  to `95×` at `n = 17`, and left it two orders of magnitude above the
  exhaustive-search row on the same base.
- **The stage ranking does not carry over unchanged**, which is why §2
  makes only this table a speed.  `crossbred-f2` is the fastest
  algebraic engine per call on every square stage cell from twelve to
  twenty-eight unknowns (§17.6), and the slowest F4-family row at
  `n = 13` (`167×` the pair table against `f4-f2`'s `102×`).  The
  §15 base descends thirteen equations in twelve unknowns and four in
  five of its systems have no decomposition; `f4-f2` refutes an
  inconsistent system at about half what a consistent one costs it
  (`0.61 ms` a call in the run; `0.8`–`2.4 ms` against `2.1`–`5.1 ms`
  on the twelve- and fourteen-unknown stage cells), while crossbred's
  `2^k` linear solves and exhaustive search cost the same either
  way.  At `n = 15` and `n =
  17`, square bases, crossbred is again the best algebraic row
  (`25×` and `24×` the pair table).
- **The exhaustive-search row is the reference's, not a candidate's.**
  Descent plus `fes-f2` sits at `1.37×` the pair table at `n = 15`
  (below it on two of six targets) and `1.70×` at `n = 17`.  Its cost
  is `2^{2n'}` points per target, over a number of targets that grows
  with the base, against a table of about `2^{2n'}/2` additions built
  once; a word XOR costs `1/415`–`1/557` of an addition on this host
  (the reports' calibrations), which is why the two meet near `n' = 8`
  on these deliberately large bases.  It is the boundary every
  algebraic engine here still has to reach inside a whole run;
  crossbred-descent sits `14×`–`18×` above it at `n = 15` and `17`.

### 17.12 Classified

| change | class | why |
|:--|:--|:--|
| F4-family engines in place of `buchberger-f2` in the descent oracle | **engineering** | `S` falls `12×`, `285×` and `1,535×` at `n = 13`, `15`, `17`, the contract's gate met with room; no boundary crossed: every algebraic row stays at least `333×` rho and `24×` the pair table on its base |
| every two-summand target decided from sixteen to twenty-two unknowns; `f4-f2` to twenty-eight, the matrix hybrids to thirty-four | **engineering** | the same systems decided within the budget where the baseline stopped at eighteen |
| `D_learn = 3` where the pair-only engine read `4` | **accounting** | the ideal did not change; the engine that measured it did (§17.7) |
| `exhaustive` counts the tests it performs | **accounting** | a modelled count replaced by the work done; wall-priced rows unaffected |
| the `vs rho` column filled | measurement | the reference §15 and §16 left empty, sixteen counted runs per instance; the plain walk, not the negation-aware one (§17.13) |
| `crossbred-f2` below `fes-f2` at twenty-six and twenty-eight unknowns | **advance against the declared reference, withdrawn** | against the vector form of the same algorithm it is `2.45`–`3.0×` (§17.9) |
| `fes-f2-wide` | **accounting**: the boundary moved | the stage reference is `3`–`4×` stronger; every stage ratio is now read against it |
| the hybrids at thirty-two and thirty-four unknowns | measurement | `2.09×` and `1.20×` the vector reference; a crossing extrapolated near thirty-five, marked as such |

### 17.13 What does not count

- Nothing in §17.6–§17.10 is a speed.  They price one oracle call; the
  speed is §17.11, and there nothing is below rho.
- §17.8's crossing is against the scalar reference only (§17.9).
- Every price here is `measured`: wall time on one host (an Intel Xeon
  at 2.1 GHz with AVX-512, four cores), with parts running concurrently
  on separate cores.  The stage cells are paired within one process and
  interleaved per target; the whole-method rows are not interleaved, and
  their host was running three other parts beside them.  Ratios between
  rows from different hosts mean nothing (§12).
- The vector reference is `3`–`4×` the scalar one, not the sixteen its
  lanes suggest; a better vector kernel would raise every stage ratio
  read against it and move any crossing further out.
- `crossbred-f2` ran at its module defaults only; §17.2 put a tuned
  crossbred in a round of its own.
- `n = 17` is one curve and two targets; its gate ratios are reported,
  not graded.
- Engines whose solutions come back in a different order walk different
  target sequences (a target with two decompositions lifts to a
  different pair), so rows on one instance are not identical workloads;
  cost per verified relation is the paired statistic for that reason,
  and the trials column shows how far the walks diverged.
- The hybrids' `D_reach` is their fixed Macaulay degree (three, or six
  on three summands), not a solving degree; only `f4-f2` reports
  `D_learn`.
- The `K` family has no instance at `n = 21`, `25`, `27` or `33`;
  those cells are `R` only.  `xl-f2` is not in the suite (§15.6).
- The first run of the thirty-four-unknown cell had no reference and is
  not read (§17.5's amendment).
- **Every `buchberger-f2` row here is the pair-only engine's.**  While
  this suite ran, `main` closed the repository's Buchberger under the
  field equations (34154ed9, 794496b4), the defect §17.1 described.
  The frozen rows measure the engine their provenance names; a rerun of
  the baseline on current code measures a different engine, and would
  be a new row.
- **The `vs rho` column is the plain walk's.**  `rho_reference` runs
  Teske's walk with no automorphism (`A = 1`).  On these random binary
  curves the eligible automorphism is negation, and the contract's
  matched reference is the negation-aware walk, about `√2` cheaper; so
  every ratio to rho here is understated by up to that factor, and the
  verdict — nothing below rho — only strengthens.  On a Koblitz curve
  the gap would be about `√(2n)`, and none of this round's whole runs
  was on one.  *§18.6 re-prices this table against the matched walk,
  and the gap was not `√2` but `3.0`–`6.1×`, mostly the plain walk's
  set-up: the pair table reads `21.0×`, `133×` and `257×` rho.*

### 17.14 Reproducing

```
cargo build --release --bin ic
S=research/ic_framework_engines_20260922

python3 $S/run.py --run-id baseline_v1                                  # §17.2
python3 $S/compare.py --run baseline_v1
python3 $S/run.py --run-id confirm_v1  --suite confirm_20260923.json    # §17.3
python3 $S/run.py --run-id strength_v1 --suite strength_20260923.json   # §17.4
python3 $S/run.py --run-id hybrids_v2  --suite hybrids_20260923.json    # §17.5
python3 $S/run.py --run-id hybrids_v3  --suite hybrids_34_20260923.json # its amendment
# and compare.py --run <id> --suite <file> for each

# A new engine against the frozen baseline, paired on this host:
python3 $S/run.py --run-id candidate_x --add-engine x
python3 $S/compare.py --run candidate_x --manifest $S/manifest.json
```

Each part's provenance names the commit, the binary's hash, the host
and the suite file's hash it ran with; `manifest.json` holds every
system's fingerprint and every cell's verdict digest from the frozen
run.

## 18. The matched rho: the reference the contract asks for

Through §17, every `vs rho` ratio on the prime and random binary
curves, and every one `ic bench` printed, divided by `rho_reference`:
Teske's r-adding walk on points, with no automorphism (`A = 1`).  The
accounting contract asks for more (`comparison_contract.rho`): "the
eligible automorphism-aware rho reference on the same subgroup and
operation unit".  On a generic prime curve and on a binary curve not
defined over a proper subfield the eligible automorphism is negation
(`A = 2`); on a Koblitz curve it is negation and the Frobenius
(`A = 2n`), which the Koblitz regime has used since Round 1 and `ic
bench` never did.  §17.13 recorded the gap.  This section closes it.
§18.1 and §18.2 were written and committed before anything below ran.

**What was known when this was written.**  Reading the frozen Round-5
headline's rho records
(`docs/ic/runs/ic-boundary-ledger-round5-headline-2026-09-22.json`,
64 runs per instance) showed the frozen walk is not at its own floor:
its walk steps over `√(πr/2)` average `2.31` on the 24-bit prime
curve, `1.48` on the `n = 27` binary curve and `1.70`–`1.72` on the
Koblitz curves' plain rows (medians `1.40` on the prime curve).  Its
walks average `846` steps against the `512` its distinguished-point
rate predicts.  The explanation offered here is a prediction, not a
finding: with walks of about `θ = √r/8` steps, a walk closes a cycle
with no distinguished point with probability about `θ²/r ≈ 2 %`, and
each such walk pays the whole `20θ` cap.  The frozen walk also pays two
scalar multiplications per walk start and about `560` operations for
its table, which is most of `S` at `r ≈ 2^{12}` (`S = 16.2` there,
§17.11).  A unit test of the new walk on `Z/(10^9 + 7)` asserts that
the plain-to-negation walk ratio lies in `[1.2, 1.65]` and each walk
within `[0.8, 1.3]` of its own floor; it passes, and its values were
not printed.  One smoke test of `ic bench` touched an evaluation
instance before this was committed: eight runs on the `W13` curve of
seed `20260922` with the provisional rule below, negation walk
`S = 2.70` against the frozen walk's `18.4` on the same eight seeds;
and one on `K_1 / GF(2^17)`, which no evaluation uses, where the
negation walk (`2.15`) came in under the signed-Frobenius walk
(`3.51`).  Neither is used below.

### 18.1 The boundaries, stated before measuring

- **The floor** is the generic one, `√(πr/2A)` operations, with the
  `A` of the curve: `S_floor = 0.886` for `A = 2`, `√(π/4n)` for a
  Koblitz curve (§3).  The walk has its own version of it: walk
  operations over `√(πr/2A)`, reported per run as
  `steps_over_expected`.
- **The reference for the reference** is what is known about adding
  walks: a random-mapping walk needs `√(πr/2A)` steps; an adding walk
  with `J` jumps needs about `1/√(1 − 1/J)` times more (Teske;
  Bernstein–Lange); the negation map's look-ahead adds about `1/(2J)`
  additions a step.  So the matched walk should sit about `1.03`–`1.1`
  above its own floor at `J = 16`, more at small `J`.
- **Everything is inside `S`**: the table, the start stride, one
  addition per walk start, the walk with its rejected look-ahead
  additions and its cycle-escape doublings, and the verification
  `[d]G = Q` of every candidate.  Not charged, as in the frozen walk:
  the hash of each point, the comparisons with recent points, and,
  under negation, one field negation and one key comparison per
  canonicalisation, which the report counts
  (`canonicalisations_uncharged`).

### 18.2 The walk, the calibration and the targets, declared before any of them ran

**The walk** (`ic_boundary::rho_reference_walk`), the same code for
both automorphism counts so that a paired run isolates the `√2`:

- van Oorschot–Wiener with distinguished points every `2^{⌊bits(r)/4⌋}`
  steps, walks of about `r^{1/4}`;
- walk `k` starts at `[k]T` for one random `T`, one addition a start;
- `J` jumps `[a_j]G + [b_j]Q`, coefficients of `min(12, 3·log₂J)` bits,
  by joint double-and-add;
- under negation (`RhoWalk::negation`): the walk moves between
  canonical representatives of `{P, −P}` (the smaller key), with the
  Wiener–Zuccherato look-ahead against the fruitless 2-cycle (a jump
  whose sum flips sign and would take the same jump again is replaced
  by the next jump), every step compared with the last sixteen points,
  and a detected cycle left by doubling its smallest-key point — a
  choice every walk entering that cycle makes identically — with a
  walk that falls back into the cycle it just left abandoned;
- the tuned walk on points (`RhoWalk::plain`) is the same code without
  the canonicalisation and the look-ahead.

**The calibration.**  The jump count `J` is the walk's one free
parameter: it trades table set-up (about `1.75·w·J` operations) against
the adding walk's randomness and, under negation, the rate of
fruitless cycles.  It is fixed on curves no evaluation uses: the
roster prime curves at `12`, `14`, `16`, `18` and `20` bits, prime
curves generated at `22`, `24` and `26` bits, and random binary curves
at `n = 13, 15, …, 27`, both from seed `0xCA11B`, `128` runs each,
`J ∈ {4, 8, 16}` for both tuned walks.  The rule: at each size the `J`
whose negation walk has the lowest mean `S`, with anything within
`3 %` of the lowest going to the larger `J`; the size thresholds between
the chosen counts go into `rho_jumps_for` and are not changed after any
evaluation run.  The code as committed with this declaration carries a
provisional rule (`4` to 18 bits, `8` to 23, `16` above) that the
calibration either confirms or replaces.

**The evaluation.**  (a) A ladder on fresh curves, seed `0xE7A1`: prime
curves generated at `12, 14, …, 26` bits (`--generated-primes`, so none
is a roster or calibration curve) and random binary curves at `n = 13,
15, …, 27`, `128` runs each, the three walks paired on the same planted
logarithms and seeds.  (b) The re-pricing of every frozen report whose
`vs rho` the scoreboard or this note quotes on a prime or random binary
curve: the Round-5 headline (`prime` 24-bit, `char2` `n = 27`) and the
§17.11 whole-method runs (`W13-*`, `W15-*`, `W17-*` in
`research/ic_framework_engines_20260922/results/baseline_v1/`).  Each
instance is rebuilt from what its report recorded, the frozen walk is
re-run on the recorded seeds first, and the matched walk then runs on
the same seeds and targets.  The index-calculus rows are not re-run;
their counts did not change.

*Amendment, written while the calibration ran and before any of its
output was read.*  The scoreboard's boundary panel, its verdict
(`3.55×`, `16.0×`) and its exponent marks for rho cite the Round-5
ladder, `docs/ic/runs/ic-boundary-ledger-round5-2026-09-22.json`, not
the headline.  So (b) covers that ladder too — eight prime rungs from
10 to 24 bits and five binary rungs from `n = 15` to `27` — and its
holdout, `…-round5-holdout-2026-09-22.json`, and target 2 covers every
one of their recorded runs.  The ladder's prime rungs from 12 to 20
bits are the roster curves the calibration also runs, on other seeds
and targets.  The rule the calibration fixes is a function of the
subgroup's bit length alone; the overlap is disclosed here rather than
avoided.

**Targets:**

1. **Correct.**  Every run of every walk, calibration and evaluation,
   recovers the planted logarithm and verifies `[d]G = Q`; walks
   abandoned in a cycle stay under `1 %` of walks.
2. **Identity.**  The frozen walk re-run on the recorded seeds
   reproduces every recorded run exactly: steps, walks, distinguished
   points, additions, doublings, scalar multiplications, the logarithm
   and `S`, in the headline and in all seven `W` reports.
3. **At its own floor.**  At every evaluation size with `r ≥ 2^{20}`,
   the negation walk's mean walk operations over `√(πr/4)` lie in
   `[0.9, 1.2]`.
4. **The `√2`, paired.**  The tuned walk on points over the negation
   walk, in walk operations, lies in `[1.25, 1.55]` pooled over the
   evaluation sizes with `r ≥ 2^{20}`, with a `95 %` bootstrap interval
   inside that range.
5. **Flat.**  The negation walk's whole `S` is at most `1.5 × 0.886`
   at every evaluation size with `r ≥ 2^{20}`, and its total operations
   fit `r^α` with `α` within `0.5 ± 0.05` over the ladder's sizes from
   `2^{20}` up, at least four of them.
6. **Matched.**  Every instance priced has no eligible automorphism
   beyond negation: prime curves with `a, b ≠ 0` or outside the
   `j = 0`, `j = 1728` congruences; binary curves whose `b` lies in no
   proper subfield.

**Abandon** the new reference, and leave every ratio against the frozen
walk, if target 1 or 2 fails; if target 3 or 4 fails the walk is
mis-built and nothing is re-priced with it; target 5 failing at the
small sizes is reported and does not stop the re-pricing, since the
reference is then still the best counted walk the repository has.

**On Koblitz curves in `ic bench`** both eligible walks run on the
same seeds, the repository's signed-Frobenius walk (`A = 2n`, priced
as the Koblitz regime prices it) and the negation walk, and the one
with the lower mean `S` prices the column; both are in the report.
That is a choice of the cheaper of two references, which can only make
the reference stronger.

**Class, fixed in advance: accounting.**  The reference moves; no
index-calculus row's counts change.  Every `vs rho` ratio is expected
to rise, and a ratio that rises is not a regression of any method.

**Inadmissible:** re-tuning `J`, the distinguished-point rate, the
window or the look-ahead after an evaluation run; dropping or re-seeding
a failed run; leaving set-up, starts or verification out of `S`;
dividing an index-calculus row by a rho mean from other seeds, targets
or curves than the row's own; and reading the frozen walk's before mark
as a property of rho rather than of that walk.

### 18.3 The calibration, and the rule it fixed

Frozen at `research/ic_rho_reference_20260923/calibration/`: sixteen
curves, 128 runs each, every run of every walk verified, `J = 4, 8, 16`
for both tuned walks.  Mean `S` of the negation walk, and what the
declared rule chose:

| bits(r) | curve | log₂ r | frozen walk `S` | negation `S`, J = 4 | J = 8 | J = 16 | chosen | its walk / own floor | tuned plain at that J, walk / own floor | J = 4 walks capped |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 12 | char2 `random-binary-n13-ba3f` | 12.0 | 14.77 | 2.56 | 3.52 | 6.31 | **4** | 1.56 | 1.28 | 1 |
| 12 | prime `bench-12bit` | 11.9 | 15.35 | 2.59 | 3.62 | 6.71 | **4** | 1.54 | 1.37 | 0 |
| 13 | char2 `random-binary-n15-b2493` | 13.0 | 12.20 | 2.17 | 2.83 | 4.93 | **4** | 1.44 | 1.20 | 1 |
| 14 | prime `bench-14bit` | 14.0 | 8.90 | 1.99 | 2.32 | 3.83 | **4** | 1.46 | 1.29 | 4 |
| 16 | prime `bench-16bit` | 16.0 | 5.66 | 1.90 | 1.85 | 2.41 | **8** | 1.33 | 1.18 | 15 |
| 17 | char2 `random-binary-n17-bbe86` | 16.0 | 5.68 | 1.87 | 1.85 | 2.41 | **8** | 1.33 | 1.15 | 14 |
| 18 | prime `bench-18bit` | 18.0 | 4.34 | 1.68 | 1.41 | 1.71 | **8** | 1.17 | 1.12 | 13 |
| 19 | char2 `random-binary-n19-b5f57e` | 18.0 | 4.05 | 1.62 | 1.47 | 1.72 | **8** | 1.23 | 1.11 | 21 |
| 19 | char2 `random-binary-n21-b10cc7b` | 19.0 | 3.58 | 1.61 | 1.29 | 1.44 | **8** | 1.14 | 1.04 | 34 |
| 20 | prime `bench-20bit` | 20.0 | 3.05 | 2.09 | 1.22 | 1.37 | **8** | 1.16 | 1.10 | 100 |
| 22 | prime `generated-22bit-3914087` | 21.9 | 2.56 | 2.01 | 1.09 | 1.14 | **8** | 1.10 | 1.14 | 173 |
| 23 | char2 `random-binary-n23-b3a296b` | 22.0 | 2.38 | 1.92 | 1.14 | 1.13 | **16** | 1.04 | 1.06 | 177 |
| 24 | char2 `random-binary-n25-b1448ace` | 24.0 | 2.36 | 3.26 | 1.13 | 1.04 | **16** | 1.05 | 1.04 | 626 |
| 24 | prime `generated-24bit-11921101` | 23.5 | 2.28 | 3.02 | 1.10 | 1.05 | **16** | 1.05 | 1.03 | 471 |
| 25 | char2 `random-binary-n27-b3bfa9c2` | 24.4 | 1.99 | 2.96 | 1.15 | 1.06 | **16** | 1.09 | 1.13 | 623 |
| 26 | prime `generated-26bit-38313677` | 25.2 | 2.71 | 3.14 | 1.06 | 0.99 | **16** | 1.03 | 1.00 | 888 |

The rule, now in `rho_jumps_for`: **four jumps to 15 bits, eight from
16 to 22, sixteen from 23.**  It replaces the provisional one at 16–18
bits (four there became eight) and at 23 (eight became sixteen); the
three sizes where two counts came within `3 %` (16, 17 and 23 bits) went
to the larger, as declared; there is no curve at 15 bits, and the
eight-over-four ratio falls from `1.17` at 14 bits to `0.97` at 16, so
15 stays with four.

What the calibration shows beyond the rule, read before any evaluation
ran:

- **Four jumps break down as walks lengthen.**  The four-jump negation
  walk runs `3.0`–`3.5×` its own floor at 24–26 bits, with 470–890 of
  its walks running to the cap: fruitless cycles longer than the
  sixteen-point window go undetected there.  Where the rule uses it
  (at most 15 bits) it caps at most four walks in 128 runs and is
  `1.4`–`1.6×` its floor, the price of a table cheap enough to be worth
  it when `√r` is 64.
- **The tuned walks sit near their floors where the table is not the
  cost**: `1.03`–`1.17×` for the negation walk from 18 bits up, and
  `1.00`–`1.14×` for the tuned walk on points.  Their paired walk ratio
  is the evaluation's target 4, not read here.
- **The frozen walk's `S` is mostly not rho.**  At 12 bits it is
  `14.8`–`15.4` where the negation walk is `2.6`; at 24–26 bits `2.0`–`2.7`
  where the negation walk is `1.0`–`1.1`.

### 18.4 The targets, graded

Frozen at `research/ic_rho_reference_20260923/` (`evaluation/`,
`reprice/`); `analyse.py` grades them and writes `analysis.json`, which
every number below comes from.

| target | declared | measured | |
|:--|:--|:--|:--|
| 1. correct | every run verified; abandoned walks under `1 %` | 21,096 runs (calibration, evaluation, re-pricing), none unverified; 88 of 821,867 tuned walks abandoned in a cycle, `0.011 %` | **met** |
| 2. identity | the frozen walk reproduces every recorded run | all 616 recorded runs of the ten reports, exactly: steps, walks, distinguished points, additions, doublings, scalar multiplications, logarithm, `S` | **met** |
| 3. at its own floor | `[0.9, 1.2]` from `2^{20}` | `1.02`–`1.20` at the six sizes from `2^{20}` (the 22-bit prime curve at `1.197`, eight jumps) | **met**, one size at the edge |
| 4. the `√2`, paired | `[1.25, 1.55]`, interval inside | pooled `1.362`, `95 %` interval `[1.296, 1.430]`; per size `1.19` at eight jumps, `1.37`–`1.45` at sixteen | **met** |
| 5. flat | `S ≤ 1.33` from `2^{20}`; `α` within `0.5 ± 0.05` over at least four sizes from `2^{20}` | `S` `0.98`–`1.20`; `α = 0.452` pooled over the six sizes from `2^{20}`, the only fit with four or more (each regime has three there: prime `0.429`, binary `0.473`) | **met on the pooled fit, at its lower edge** |
| 6. matched | no automorphism beyond negation | none on the sixteen evaluation curves or the 26 re-priced instances: every prime curve has `a, b ≠ 0`, no binary `b` lies in a proper subfield | **met** |

Target 5 is the one to read with care.  The exponent is below one half
because set-up is still visible at these sizes: `S` falls from `1.20` at
`2^{21}` to `0.98` at `2^{25}`.  The table, the stride and the
distinguished-point tail are a few hundred operations against `√r` of a
few thousand.  Over every size of the ladder the fit is `0.394` (prime)
and `0.397` (binary).  So "flat" holds in the sense declared, `S` inside
`1.5×` the floor from `2^{20}`, and the pooled exponent clears its bound
by `0.002`.  Neither is an asymptotic statement; the asymptote the walk
is heading for is its floor times the adding-walk penalty, about
`0.886 × 1.03`–`1.07`.

### 18.5 The evaluation ladder

Seed `0xE7A1`, sixteen curves no other run used, 128 runs a walk, all
three walks on the same planted logarithms and seeds:

| regime | curve | log₂ r | J | frozen walk `S` | tuned, points `S` | **negation `S`** (median) | negation walk / own floor | points / negation, walk ops | frozen / negation, `S` |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| prime | `generated-12bit-3529` | 11.8 | 4 | 15.05 | 3.00 | **2.80** (2.57) | 1.727 | 1.104 | 5.37× |
| prime | `generated-14bit-14389` | 13.8 | 4 | 9.46 | 2.37 | **2.09** (1.95) | 1.540 | 1.160 | 4.52× |
| prime | `generated-16bit-53411` | 15.7 | 8 | 6.21 | 2.28 | **1.82** (1.75) | 1.235 | 1.388 | 3.41× |
| prime | `generated-18bit-212099` | 17.7 | 8 | 4.42 | 1.77 | **1.49** (1.42) | 1.216 | 1.246 | 2.97× |
| prime | `generated-20bit-563153` | 19.1 | 8 | 4.14 | 1.62 | **1.38** (1.29) | 1.272 | 1.205 | 3.00× |
| prime | `generated-22bit-2434423` | 21.2 | 8 | 3.01 | 1.41 | **1.20** (1.19) | 1.197 | 1.193 | 2.50× |
| prime | `generated-24bit-16483309` | 24.0 | 16 | 2.56 | 1.49 | **1.06** (1.04) | 1.075 | 1.445 | 2.42× |
| prime | `generated-26bit-42652151` | 25.3 | 16 | 2.45 | 1.34 | **0.98** (0.91) | 1.024 | 1.391 | 2.50× |
| char2 | `random-binary-n13-b168f` | 10.4 | 4 | 21.10 | 3.53 | **3.17** (3.11) | 1.408 | 1.219 | 6.66× |
| char2 | `random-binary-n15-b157a` | 14.0 | 4 | 9.06 | 2.40 | **2.08** (1.93) | 1.550 | 1.180 | 4.35× |
| char2 | `random-binary-n17-b10103` | 14.4 | 4 | 9.03 | 2.25 | **2.06** (1.97) | 1.590 | 1.096 | 4.38× |
| char2 | `random-binary-n19-b770f7` | 18.0 | 8 | 4.39 | 1.73 | **1.38** (1.32) | 1.141 | 1.322 | 3.18× |
| char2 | `random-binary-n21-b968e` | 18.0 | 8 | 4.28 | 1.78 | **1.47** (1.40) | 1.235 | 1.263 | 2.91× |
| char2 | `random-binary-n23-b480f9a` | 22.0 | 16 | 2.87 | 1.51 | **1.13** (1.11) | 1.047 | 1.395 | 2.53× |
| char2 | `random-binary-n25-b1a9b772` | 22.4 | 16 | 2.27 | 1.51 | **1.14** (1.10) | 1.082 | 1.374 | 1.99× |
| char2 | `random-binary-n27-b227199e` | 26.0 | 16 | 2.28 | 1.46 | **1.06** (1.06) | 1.127 | 1.392 | 2.15× |

Reading it:

- **The reference moved `2.0`–`6.7×`.**  Most of that is the frozen
  walk's tuning, not the automorphism.  Frozen over tuned-on-points is
  `1.5`–`6.0×`, and tuned-on-points over negation is `1.10`–`1.45×` in
  walk operations.  The `√2` is the second factor; the first is the
  frozen walk's per-walk scalar multiplications, its `560`-operation
  table and its trapped walks (§18 opening).
- **The `√2` is paid for by the jump count.**  At sixteen jumps the
  paired walk ratio is `1.37`–`1.45`, the `√2` less the look-ahead's
  `1/32`.  At eight it is `1.19`–`1.39`, and at four `1.10`–`1.22`,
  where the look-ahead's `1/8` and the fruitless cycles it cannot stop
  eat most of it.  The calibration chose four jumps at these sizes
  anyway, because a table costs more there than the `√2` saves.
- **The fitted exponents** of the whole cost over all eight sizes:
  frozen walk `0.311` (prime) and `0.287` (binary), tuned on points
  `0.417` and `0.420`, negation `0.394` and `0.397`, every `R²` above
  `0.96`.  The frozen walk's exponent was its set-up's; the tuned walks'
  are closer to one half and still under it, for the set-up reason
  §18.4 gives.

### 18.6 The frozen reports, re-priced

Every recorded run of the frozen walk replayed exactly (target 2), so
these ratios differ from the frozen ones only in the reference.  The
matched walk's mean `S`, on the same curves, seeds and targets:

| report | curve | log₂ r | frozen rho `S` | matched rho `S` | frozen / matched |
|:--|:--|--:|--:|--:|--:|
| Round-5 ladder | `bench-10bit` … `bench-20bit` | 9.7–20.0 | 22.66 · 18.58 · 8.91 · 6.36 · 4.20 · 2.68 | 3.60 · 2.85 · 2.03 · 1.74 · 1.45 · 1.29 | 2.1–6.5× |
| | `generated-22bit-3290411` · `-24bit-10935329` | 21.7 · 23.4 | 2.98 · 3.93 | 1.22 · 1.04 | 2.45× · 3.79× |
| | `random-binary-n15` · `n18` · `n21` · `n24` · `n27` | 14.0–24.4 | 8.36 · 7.56 · 3.18 · 2.52 · 2.31 | 2.01 · 2.01 · 1.36 · 1.18 · 0.97 | 2.1–4.2× |
| Round-5 headline | `generated-24bit-10935329` · `random-binary-n27-b845462` | 23.4 · 24.4 | 3.28 · 2.18 | 1.05 · 1.03 | 3.12× · 2.12× |
| Round-5 holdout | two prime, two binary | 21.0–26.0 | 2.96 · 2.66 · 2.43 · 1.99 | 1.12 · 0.88 · 0.96 · 0.90 | 2.2–3.0× |
| §17.11 `W13` | three curves | 11.0–12.0 | 16.2 · 18.7 · 14.6 | 2.65 · 3.10 · 2.44 | 6.0–6.1× |
| §17.11 `W15` | three curves | 12.4–14.0 | 13.1 · 9.2 · 8.9 | 2.21 · 2.06 · 1.76 | 4.5–5.9× |
| §17.11 `W17` | one curve | 16.0 | 5.21 | 1.74 | 2.99× |

**The boundary panel's rows** (the Round-5 ladder's largest prime and
binary rungs, means over three targets), `vs rho` frozen → matched:

| variant | 24-bit prime | `n = 27` binary |
|:--|--:|--:|
| Semaev `S₃` roots / `S₄` pairs-and-solve | 1,300× → 4,922× | 49,906× → 119,465× |
| direct subtraction, `m = 2` | 518× → 1,962× | — |
| meet in the middle, `m = 2` / `m = 3` | 68.3× → 258× · 14.5× → 54.7× | — / 89.1× → 213× |
| … negation-folded table | 63.2× → 239× (`m = 2`) · 9.41× → 35.6× (`m = 3`) | 711× → 1,702× (`m = 2`) · 82.7× → 198× (`m = 3`) |
| … and walk targets | 6.32× → 23.9× (`m = 2`) · 8.06× → 30.5× (`m = 3`) | 18.0× → 43.1× (`m = 2`) · 73.2× → 175× (`m = 3`) |
| … and a base at the family optimum | **3.55× → 13.5×** | **16.0× → 38.2×** |

The best row on every other rung, frozen → matched: prime `1.01×`,
`0.84×`, `1.47×`, `1.84×`, `2.00×`, `3.99×`, `4.28×` from 10 to 22 bits
become `6.3×`, `5.5×`, `6.4×`, `6.7×`, `5.8×`, `8.3×`, `10.5×`; binary
`3.15×`, `2.82×`, `7.66×`, `10.8×` at `n = 15`–`24` become `13.1×`,
`10.6×`, `18.0×`, `23.0×`.  **The frozen ladder had a row below its
rho**, the 12-bit prime rung's `mitm_m2_negfold_walk` at `0.84×`.  No
earlier round quoted it: the page showed the largest rungs, and its
legend warned that the frozen `vs rho` "flatters index calculus at
small sizes".  It was the plain walk's set-up, `S = 18.6` on a
subgroup of 3,889 elements, and not a crossing.  Against the matched
walk it is `5.5×`.  The headline's best rows go from `4.19×` and
`21.4×` to `13.1×` and `45.3×`; the holdout's from `4.80×`, `4.74×`,
`11.8×` and `16.8×` to `12.7×`, `14.4×`, `29.8×` and `37.1×`.

**§17.11 re-priced**, medians over curves and targets as there:

| row | `n = 13`: / rho, frozen → matched | `n = 15` | `n = 17` |
|:--|--:|--:|--:|
| pair table (the oracle reference) | 3.5× → **21.0×** | 29.1× → **133×** | 86.1× → **257×** |
| descent + `fes-f2` | 22.4× → 135× | 40.0× → 178× | 146× → 437× |
| descent + `exhaustive` | 50.1× → 302× | 710× → 3,163× | 5,413× → 16,163× |
| descent + `crossbred-f2` | 509× → 3,093× | 722× → 3,216× | 2,098× → 6,264× |
| descent + `inherited-f4` | 344× → 2,089× | 1,949× → 8,676× | 5,513× → 16,460× |
| descent + `f4-f2` | 333× → 2,024× | 2,961× → 13,185× | 8,186× → 24,443× |
| descent + `matrix-f4` / `matrix-f5` | 362× / 387× → 2,180× / 2,335× | 3,972× / 4,075× → 17,685× / 18,146× | 14,873× / 14,963× → 44,410× / 44,680× |
| descent + `sat-cdcl` (first solution) | 9,035× → 54,427× | 1.04×10⁵× → 5.64×10⁵× | 8.55×10⁵× → 2.55×10⁶× |
| descent + `buchberger-f2` | 3,832× → 23,091× | 8.45×10⁵× → 3.76×10⁶× | 1.26×10⁷× → 3.75×10⁷× |

**The rho exponent on the Round-5 ladder**, the mark the page's exponent
panel draws, moves from `0.274` to `0.372` (prime, eight rungs) and
from `0.302` to `0.392` (binary, five).  The reference still carries
set-up at the ladder's small rungs, so it is still under one half.
Every index-calculus exponent on that panel is fitted to the IC rows'
own totals and does not move.

### 18.7 Classified

| change | class | why |
|:--|:--|:--|
| the matched walk prices every prime and random binary `vs rho`, `ic bench`'s included | **accounting** | the reference moved and no index-calculus count changed; every ratio rises, `2.1`–`6.5×` on the re-priced reports; nothing changes class, and the verdict — nothing below rho — now holds on every rung, the frozen `0.84×` included |
| the tuned walk's set-up: stride starts, distinguished points every `r^{1/4}`, a table sized to `r` | engineering of the reference | frozen over tuned-on-points `1.5`–`6.0×` on the evaluation ladder; a better-built reference, not a property of rho |
| the negation map | the automorphism the floor always credited | `1.36×` pooled and paired (§18.4, target 4); `√2` less the look-ahead |
| `ic bench` on a Koblitz curve: the cheaper of the signed-Frobenius and negation walks | accounting | the negation walk is cheaper to `r ≈ 2^{18}`, the signed walk from `2^{20.5}` (§18.8), and the reference is whichever a generic attacker would run |

### 18.8 An open item: the Koblitz reference does not price its canonicalisation

The Koblitz diagnostic
(`research/ic_rho_reference_20260923/diagnostic/koblitz.json`, seven
curves, 64 runs each, all verified) was run to support `ic bench`'s
rule, and found two things about the signed-Frobenius walk that the
Koblitz regime has used as its reference since Round 1.

- **Set-up at small rungs.**  The walk starts 32 parallel walks and a
  16-jump table by scalar multiplication, about 97 of them, so below
  about `2^{20}` it is far weaker than the negation walk:

  | curve | log₂ r | signed-Frobenius `S` | negation `S` |
  |:--|--:|--:|--:|
  | `K_0 / GF(2^13)` | 11.0 | 13.5 | **2.78** |
  | `K_1 / GF(2^17)` | 16.0 | 3.56 | **1.80** |
  | `K_1 / GF(2^19)` | 18.0 | 2.10 | **1.38** |
  | `K_0 / GF(2^31)` | 20.5 | **1.10** | 1.19 |
  | `K_1 / GF(2^23)` | 22.0 | **0.88** | 1.17 |
  | `K_0 / GF(2^37)` | 27.8 | **0.42** | 0.93 |
  | `K_0 / GF(2^41)` | 39.0 | **0.16** | 1.03 |

  So the Koblitz ladder's `vs rho` below about `2^{20}` flatters index
  calculus by up to `4.8×`, as the prime and binary ladders' did.
- **Canonicalisation is free.**  Every step canonicalises by `n`
  Frobenius maps.  The counters read 13.0, 17.0, 19.1, 23.2, 31.2, 37.4
  and 41.1 maps a step, and the unit charges none of them.  The
  index-calculus rows on the same curves pay for theirs: the folded
  tables charge a canonicalisation per probe, `0.063`–`0.074` additions
  at the pinned ratio.  Priced as implemented, at the pinned Frobenius
  ratio, the maps would add `1.7 %` (`n = 13`) to `101 %` (`n = 41`) to
  the walk's cost.  Priced with the primitive the index-calculus rows
  are charged for, one canonicalisation a step, they would add about
  `7 %`.  A walk that needs no per-step canonicalisation would pay less
  than that: Bailey et al.'s `P + φ^j(P)`, with `j` read from an
  orbit invariant, commutes with the Frobenius and canonicalises only
  at distinguished points.  It is not built here.

So every Koblitz `vs rho` on the ledger and the page is off in both
directions: generous below `2^{20}`, harsh above.  At the boundary
panel's `n = 41` row the second effect is about `1.07×` on the affine
unit.  The collection thread's `1.17×` (`K_0 / GF(2^41)`, a
batched-addition unit, reference `S` = rho iterations over `√r`)
carries the same asymmetry.  A batched addition shares its inversion,
so a canonicalisation is a larger fraction of it than of an affine
addition, and that correction is larger than `7 %`, possibly enough to
take `1.17×` below one.  That is not measured here and is not a result.
Until a declared round prices the Koblitz reference, `1.17×` is an
upper bound, not a margin.  The instruction-count panel
(`ic-crossover-20260922`) is not affected: an instruction count
charges rho's canonicalisation like everything else.

*Settled in §19.*  Priced in the batched unit, a canonical step costs
`2.74` units, so the correction guessed at above is real and large.  But
a larger error ran the other way: `1.17×` compared 32 targets solved
together with rho solving one.  Against batch rho at `k = 32` it reads
`6.51×`, and `6.38×` with rho's step and the index-calculus build both
priced (§19.5).  Nothing went below one.

### 18.9 What does not count

- The matched walk's `S` at toy sizes is still partly set-up: `3.6` at
  `2^{9.7}`.  It is the cheapest counted walk the repository has, not
  the asymptotic `0.886`.
- Target 5 is met on the pooled fit only, by `0.002`; §18.4 says why.
- The calibration's roster curves are also rungs of the re-priced
  Round-5 ladder (§18.2's amendment).  The rule depends on bit length
  alone.
- The matched walk leaves the same native work uncharged as the frozen
  walk — the hash per point — plus, under negation, one field negation
  and one key comparison per canonicalisation and up to sixteen key
  comparisons a step against recent points.  All are counted, none is
  a group operation, and each is a small fraction of an addition.
- The re-pricing re-uses the frozen index-calculus rows.  Their counts
  are the frozen ones, which is why this is accounting.
- `ic bench` on a Koblitz curve takes the cheaper of two walks by mean
  `S`.  That is a choice between references, and it can only make the
  reference stronger.
- The Koblitz reference in `ic boundary` is unchanged, with both of
  §18.8's defects.

### 18.10 Reproducing

```
cargo build --release --bin ic
sh research/ic_rho_reference_20260923/run.sh      # calibration, evaluation, re-pricing
python3 research/ic_rho_reference_20260923/analyse.py
```

`ic rho --koblitz-degrees 13,17,19,23,31,37,41 --runs 64 --seed 59297`
is §18.8's diagnostic.  Every report names the binary's hash and the
commit it ran from; each re-pricing report names its source file's
blake3.

## 19. The Koblitz references: what the collection thread's `vs rho` compares

§18.8 left one item open: the Koblitz reference prices its per-step
canonicalisation at zero.  Reading how the collection thread builds its
`vs rho` turned up a larger mismatch as well, and this section tests
both.  The sources are `scripts/ic_e2e_benchmark.py::rho_s_of` and the
frozen `docs/ic/runs/koblitz-*` runs.  §19.1 was written and committed
before anything below it ran; the probe it discloses is the only
measurement made first.

**Three facts, read from the code and the frozen runs.**

1. **Thirty-two targets against one.**  Every collection-thread figure
   is `total operations / (32·√r)`: `1.17×`, `1.33×`, `8.65×`, `19.3×`
   and the `n = 61` panel.  The 32 targets share one factor-base
   selection, one table build and one collection.  The descent adds
   about two operations a target (`koblitz-select-packed-20260922.json`:
   54 over 32).  The reference beside each figure is single-target rho,
   run separately for each target and averaged.  The repository's own
   rule is `src/ecc_safety.rs::check_multi_target_margin` and §5 of
   `RESEARCH_ECC2K130_RR_SOLVER_PANEL.md` (Kuhn–Struik;
   Galbraith–Lin–Scott).  It says `k` logarithms in one group cost about
   `√(k·r)` in total, so per-target rho falls as `1/√k`, and "a batch
   win over `k` independent rho runs measures the baseline, not the
   algorithm."
2. **The reference's step is counted as its addition.**
   `rho_S = rho_group_additions / (targets·√r)`, while each step also
   canonicalises over the `2n` conjugates.
3. **The index-calculus build is counted as its additions.**  One
   addition per stored pair, while each stored pair is also
   canonicalised: its key is the normal-basis canonical form of the
   sum's abscissa.

**The probe**, run before this was written: `examples/koblitz_rho_price.rs`,
16 targets.  The thread's unit is one batched affine addition, `42.8 ns`
at `n = 41` and `54.1 ns` at `n = 53` on this host.  In that unit:

- The implemented walk's step costs `20.9` and `23.3` units.  Of that,
  `18.3` and `21.1` are its canonicalisation, a chain of about `1.5n`
  dependent squarings.
- A normal-basis rotation with `y` lifted by squarings costs `7.0` and
  `8.1`.
- Priced by time, the implemented walk is `S = 3.95` at `n = 41`,
  twenty-one times the `0.185` the thread counts for it.

That is a slow canonicalisation, not rho.  A walk that canonicalises by
table pays a fraction of an addition a step, and so does one that never
canonicalises at all (Bailey et al.'s `P ↦ P + φ^j(P)`, with `j` an
orbit invariant).  Pricing rho as implemented would manufacture a
crossing.

### 19.1 Declared before anything below ran

**Boundaries.**  For one target the floor is `√(π/4n)` in `S`.  For `k`
targets the generic bound falls as `1/√k` per target (Yun's `Ω(√(kr))`
in total).  The reference for a `k`-target figure is batch rho at the
same `k`, measured.

**Measurements.**

- **(a) Batch rho.**  A batch mode of the tuned walk (§18.2) on the
  signed-Frobenius classes (`A = 2n`), with counted operations per
  target and every logarithm verified:
  - jumps in `G` only;
  - one table of distinguished points shared by every target;
  - targets solved in sequence, so a later walk can finish on an
    earlier target's trail.

  It runs on the thread's two rung curves, `K_0 / GF(2^41)` and
  `K_0 / GF(2^53)`, at `k = 1, 4, 16, 32`, with 16 batches each on
  fresh targets.  On the `n = 61` panel's curve it runs at `k = 1` and
  `32`, with 8 batches.
- **(b) Step price,** in the thread's unit, measured in the same
  process.  Two well-built walks, each timed by its per-step primitives:
  - a canonical walk: one batched addition plus a table-driven
    canonicalisation (normal-basis coordinates, least rotation, and two
    Frobenius powers applied by table);
  - Bailey et al.'s walk: one batched addition plus the coordinates, a
    popcount and two table-applied Frobenius powers.
- **(c) Build price,** in the thread's unit, measured in the same
  process: the folded table's build time per stored pair, on the
  thread's `n = 41` base (`|F| = 15,744`).

**Targets:**

1. **Correct.**  Every target of every batch is recovered and verified.
2. **The batch law holds for this walk.**  Per-target cost at `k = 32`
   over `k = 1` lies in `[0.09, 0.35]`, a factor of two either side of
   `1/√32`.
3. **Re-read.**  Every 32-target Koblitz figure on the page is re-read
   against batch rho at `k = 32` on its own curve.  That is measured at
   `n = 41` and `53`, and at `n = 61` too if (a) completes there.
   Beside each re-read go (b) and (c) as measured corrections, and
   beside those the one-target (cold) figure.

**Class: accounting.**  No index-calculus count changes; the reference
is matched to the problem the figures solve.

**Inadmissible:**

- pricing the reference at its implemented canonicalisation, the
  probe's `20.9` units a step;
- using the `1/√k` formula where the batch measurement exists;
- comparing a 32-target figure with batch rho at a different `k`;
- dropping a failed batch;
- using the batch reference for a cold figure.

**Abandon** the re-pricing if target 1 fails.  If target 2 fails, the
batch mode is mis-built, and the page gets the formula's reading,
marked as a model, instead of a measurement.

### 19.2 What ran

The walk of §18 is now generic over the classes it moves between
(`RhoClasses`: points, `{P, −P}`, and on a Koblitz curve `{±φ^t(P)}`).
Its step, look-ahead, cycle escape and cap are one piece of code, shared
by the single-target walk and the batch.  The single-target walk replays
all 3,072 runs of eight instances of the frozen §18 evaluation ladder
exactly — steps, operations and every counter — so §18's numbers stand as
measured.

- **(a)** `rho_batch_with`, run by `ic rho --batch-koblitz`: the tuned
  walk on `SignedFrobeniusClasses`.  A point's representative is the
  conjugate whose abscissa has the least normal-basis rotation, carried
  there by table-applied Frobenius powers (`FrobeniusPowers`), then the
  sign with the smaller ordinate.  It is `[±λ^t]P`, which a test checks on
  every member of the class.  Sixteen jumps `[c_j]G` with full-size `c_j`,
  distinguished points every `2^{⌊bits(r)/4⌋}` steps, one table for the
  batch.
- **(b) and (c)**: `examples/koblitz_reference_prices.rs`, one thread,
  seven interleaved rounds, every quantity divided by the unit measured in
  its own round; medians with the range over rounds.

Everything is in `research/ic_rho_koblitz_20260923/`, run from commit
`3e8dd352` with a clean tree (`provenance.txt`) and seed `0xBA7C4`.  Three
things ran that §19.1 did not declare, and each is labelled in its file:

- the build price on the base of every other quoted figure, so that no
  base carries another's price.  One of these, `n = 41, |F| = 16,400`, was
  added after the first pass and before any analysis;
- a lone affine addition, to price the canonicalisation in the
  three-regime ledger's unit;
- `add_pairwise`, the batched addition a parallel walk actually runs.

### 19.3 (a) Batch rho: the batch law holds — targets 1 and 2 met

| curve | k | batches | S per target (95 % CI) | over k = 1 (95 % CI) | batch law | over floor | own / earlier trail | ok |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|
| K_0 / GF(2^41) | 1 | 16 | 0.1429 [0.1030, 0.1828] | 1.000 [1.000, 1.000] | 1.000 | 1.032 | 16 / 0 | 16/16 |
| K_0 / GF(2^41) | 4 | 16 | 0.0742 [0.0664, 0.0819] | 0.519 [0.409, 0.705] | 0.547 | 0.536 | 28 / 36 | 64/64 |
| K_0 / GF(2^41) | 16 | 16 | 0.0396 [0.0367, 0.0425] | 0.277 [0.220, 0.372] | 0.280 | 0.286 | 40 / 216 | 256/256 |
| K_0 / GF(2^41) | 32 | 16 | 0.0300 [0.0287, 0.0312] | 0.210 [0.168, 0.280] | 0.199 | 0.216 | 43 / 469 | 512/512 |
| K_0 / GF(2^53) | 1 | 16 | 0.1052 [0.0682, 0.1422] | 1.000 [1.000, 1.000] | 1.000 | 0.864 | 16 / 0 | 16/16 |
| K_0 / GF(2^53) | 4 | 16 | 0.0626 [0.0552, 0.0700] | 0.595 [0.435, 0.830] | 0.547 | 0.515 | 27 / 37 | 64/64 |
| K_0 / GF(2^53) | 16 | 16 | 0.0340 [0.0313, 0.0367] | 0.323 [0.239, 0.447] | 0.280 | 0.279 | 42 / 214 | 256/256 |
| K_0 / GF(2^53) | 32 | 16 | 0.0252 [0.0240, 0.0263] | 0.239 [0.177, 0.329] | 0.199 | 0.207 | 41 / 471 | 512/512 |
| K_0 / GF(2^61) | 1 | 8 | 0.1151 [0.0503, 0.1798] | 1.000 [1.000, 1.000] | 1.000 | 1.014 | 8 / 0 | 8/8 |
| K_0 / GF(2^61) | 32 | 8 | 0.0219 [0.0205, 0.0233] | 0.190 [0.131, 0.320] | 0.199 | 0.193 | 21 / 235 | 256/256 |

- **Target 1 met.**  All 1,960 targets were recovered and verified, and
  no batch was dropped.
- **Target 2 met on all three curves.**  Per-target cost at `k = 32` over
  `k = 1` is `0.210` `[0.168, 0.280]`, `0.239`
  `[0.177, 0.329]` and `0.190` `[0.131, 0.320]` (batches
  resampled), against Kuhn–Struik's `0.199` and the declared
  `[0.09, 0.35]`.  The sizes in between follow the law as well: `0.519`
  and `0.277` at `n = 41`, where the law gives `0.547` and `0.280`.
- **Later targets finish on earlier trails.**  At `k = 32`, 469 of 512
  targets at `n = 41` did.
- **Fruitless cycles are rare on these classes.**  A Frobenius shift
  `t ≠ 0` never undoes the jump just taken; only the `t = 0` negation
  does, and the look-ahead catches it.  Between 0 and 22 short cycles
  were detected per size, and 4 walks in the whole run reached the cap,
  all charged.
- **One target alone** costs `1.03`, `0.86` and `1.01` times its floor,
  with the wide intervals that sixteen (at `n = 61`, eight) single rho
  runs carry.  The thread's single-target references were the
  repository's implemented walk, counted the same way.  They sat at
  `1.20–1.23`, `1.78` and `1.56` times the same floors.

### 19.4 (b) and (c): the step and the stored pair, priced

A step, in units of one batched addition (`add_many` over 1,024 points):

| curve | unit (ns) | canonical step | Bailey step | add_pairwise | affine add | canonicalisation in affine units |
|:--|--:|--:|--:|--:|--:|--:|
| n = 41 | 43.8 | 2.74 [2.59, 2.89] | 1.39 [1.38, 1.54] | 1.06 | 14.8 | 0.117 |
| n = 53 | 50.6 | 2.83 [2.75, 2.87] | 1.38 [1.36, 1.42] | 0.98 | 18.7 | 0.098 |
| n = 61 | 54.0 | 2.92 [2.72, 3.04] | 1.41 [1.37, 1.44] | 1.04 | 21.0 | 0.091 |

A stored pair of the folded table, in the same unit:

| base | orbits | stored pairs | units per stored pair | range | declared |
|:--|--:|--:|--:|--:|:--|
| n = 41, F = 5,248 | 64 | 170,560 | 4.47 | [3.80, 4.96] | no |
| n = 41, F = 15,744 | 192 | 1,519,296 | 6.14 | [5.48, 6.76] | yes |
| n = 41, F = 16,400 | 200 | 1,648,200 | 5.69 | [4.69, 6.95] | no |
| n = 53, F = 15,264 | 144 | 1,106,640 | 5.87 | [5.23, 6.71] | no |
| n = 61, F = 6,832 | 56 | 194,712 | 4.96 | [3.80, 6.77] | no |
| n = 61, F = 9,760 | 80 | 395,280 | 5.18 | [4.70, 5.64] | no |
| n = 61, F = 12,688 | 104 | 666,120 | 5.56 | [5.34, 5.92] | no |
| n = 61, F = 18,544 | 152 | 1,418,616 | 5.86 | [5.12, 7.81] | no |

- **A canonical step costs `2.74–2.92` units, not one.**  The table-driven
  canonicalisation costs more than the addition it follows (`1.74–1.92`
  units), because the least rotation was, at `3e8dd352`, a serial scan
  over `n` rotations (see §19.8 for what main has changed since).
- **Bailey et al.'s step costs `1.38–1.41` units**, because that walk
  never canonicalises.  Its step *count* was not measured here, so any
  figure priced at its step is a model, and is marked as one.
- **The unit is fair to rho.**  The batched addition of 1,024 independent
  walks (`add_pairwise`) costs `0.98–1.06` units.
- **A stored pair costs `4.5–6.1` units, not one**, and
  `6.14` `[5.48, 6.76]` on the declared base.  The folded
  build makes two passes over every row, one to count and one to fill.
  Each pass computes a batched addition and a canonical key per pair,
  then scatters the pair into its bucket.  A one-pass build would roughly
  halve that; it was not measured.
- **In the three-regime ledger's unit**, a lone affine addition with its
  own inversion (`14.8–21.0` units here), the canonicalisation costs
  `0.09–0.12` of an addition.  §18.8 had estimated about 7%.  Those rows
  are single-target, so they are not re-read here (§19.7).

### 19.5 The re-read — target 3 met

Every 32-target Koblitz figure on the page, from its frozen file
(`analyse.py` lists the file and field of each).  The columns are:

- **quoted**: as the page had it, against single-target rho;
- **re-read**: against batch rho at `k = 32` on the same curve, both sides
  counted as the thread counts, one unit per rho step and per stored pair;
- **+ (b) canonical**: rho's step priced at the canonical walk it ran;
- **+ (b) Bailey, model**: priced at Bailey's step, holding the step count
  (a model, see §19.4);
- **+ (c)**: the index-calculus build priced per stored pair, where the
  table is folded and its base was measured;
- **+ (b) and (c)**: both, with the canonical step;
- **cold**: one target alone — the whole pipeline less the other 31
  targets' descents, against the `k = 1` measurement.

`≥` marks rows whose `S` leaves selection and the linear algebra unpriced.

| panel | row | n | quoted | re-read (k = 32) | + (b) canonical | + (b) Bailey, model | + (c) | + (b) and (c) | cold | cold + (b), (c) |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| select-packed-20260922 | aimed, BigUint base build | 41 | 1.33× | **7.41×** | 2.70× | 5.32× | 18.4× | 6.71× | 49.7× | 45.0× |
| select-packed-20260922 | aimed, packed base build | 41 | 1.17× | **6.51×** | 2.37× | 4.67× | 17.5× | 6.38× | 43.7× | 42.8× |
| collection-aim-20260922 | original rung, full scan | 41 | 8.50× | **47.3×** | 17.2× | 33.9× | 48.1× | 17.5× | 317× | 118× |
| collection-aim-20260922 | + collection window 164 | 41 | 3.48× | **19.4×** | 7.06× | 13.9× | 20.2× | 7.36× | 130× | 49.4× |
| collection-aim-20260922 | + wider base 16,400 | 41 | 2.78× | **15.4×** | 5.63× | 11.1× | 26.3× | 9.60× | 104× | 64.4× |
| collection-aim-20260922 | + lean relation target | 41 | 2.04× | **11.3×** | 4.13× | 8.13× | 22.2× | 8.10× | 76.0× | 54.3× |
| collection-aim-20260922 | + stop when determined | 41 | 1.75× | **9.73×** | 3.55× | 6.98× | 20.7× | 7.55× | 65.2× | 50.7× |
| collection-aim-20260922 | + aim at uncovered columns (relabelling, not kept) | 41 | 1.71× | **9.53×** | 3.47× | 6.84× | 20.5× | 7.48× | 63.9× | 50.2× |
| collection-aim-20260922 | + aim at least-mentioned | 41 | 1.33× | **7.42×** | 2.70× | 5.32× | 18.4× | 6.71× | 49.8× | 45.0× |
| collection-aim-20260922 | + aim at least-mentioned (holdout, seeds 900-931) | 41 | 1.33× | **7.42×** | 2.70× | 5.32× | 18.4× | 6.71× | 49.8× | 45.0× |
| phase-prices-20260921 | k0n41-subgroup-rung | 41 | 8.75× | **49.8×** | 18.2× | 35.7× | 50.6× | 18.5× | 334× | 124× |
| phase-prices-20260921 | k0n53-subgroup-rung | 53 | 19.4× | **167×** | 59.1× | 121× | 169× | 59.7× | 1,257× | 449× |
| collection-window-20260921 | baseline_full_scan | 53 | 19.1× | **165×** | 58.4× | 119× | 167× | 59.0× | 1,242× | 443× |
| collection-window-20260921 | window_1908 | 53 | 7.67× | **66.2×** | 23.4× | 47.9× | 67.7× | 23.9× | 485× | 175× |
| collection-window-20260921 | window_954 | 53 | 7.43× | **64.2×** | 22.7× | 46.4× | 65.7× | 23.2× | 469× | 170× |
| collection-window-20260921 | window_477 | 53 | 7.12× | **61.4×** | 21.7× | 44.4× | 62.9× | 22.2× | 448× | 162× |
| collection-window-20260921 | window_477_holdout | 53 | 6.99× | **60.4×** | 21.4× | 43.7× | 61.8× | 21.9× | 448× | 162× |
| collection-window-20260921 | baseline_full_scan | 41 | 8.50× | **47.3×** | 17.2× | 33.9× | 48.1× | 17.5× | 317× | 118× |
| collection-window-20260921 | window_164 | 41 | 3.48× | **19.4×** | 7.06× | 13.9× | 20.2× | 7.36× | 130× | 49.4× |
| probe-volume-20260921 | k0n41-subgroup-rung compact | 41 | ≥ 11.9× | **≥ 67.5×** | ≥ 24.6× | ≥ 48.4× | — | — | ≥ 453× | — |
| probe-volume-20260921 | k0n41-subgroup-rung folded | 41 | ≥ 8.65× | **≥ 49.2×** | ≥ 17.9× | ≥ 35.3× | ≥ 50.1× | ≥ 18.3× | ≥ 330× | ≥ 122× |
| probe-volume-20260921 | k0n53-subgroup-rung compact | 53 | ≥ 23.7× | **≥ 205×** | ≥ 72.5× | ≥ 148× | — | — | ≥ 1,548× | — |
| probe-volume-20260921 | k0n53-subgroup-rung folded | 53 | ≥ 19.3× | **≥ 167×** | ≥ 59.0× | ≥ 121× | ≥ 168× | ≥ 59.5× | ≥ 1,254× | ≥ 447× |
| tier-crossover-20260921 | 6832 points full | 61 | ≥ 12.1× | **≥ 98.0×** | ≥ 33.6× | ≥ 69.4× | — | — | ≥ 505× | — |
| tier-crossover-20260921 | 6832 points compact | 61 | ≥ 12.8× | **≥ 104×** | ≥ 35.5× | ≥ 73.4× | — | — | ≥ 529× | — |
| tier-crossover-20260921 | 6832 points folded | 61 | ≥ 15.2× | **≥ 124×** | ≥ 42.3× | ≥ 87.4× | ≥ 124× | ≥ 42.4× | ≥ 601× | ≥ 206× |
| tier-crossover-20260921 | 9760 points full | 61 | ≥ 12.5× | **≥ 101×** | ≥ 34.6× | ≥ 71.5× | — | — | ≥ 552× | — |
| tier-crossover-20260921 | 9760 points compact | 61 | ≥ 12.2× | **≥ 98.4×** | ≥ 33.7× | ≥ 69.6× | — | — | ≥ 535× | — |
| tier-crossover-20260921 | 9760 points folded | 61 | ≥ 13.9× | **≥ 113×** | ≥ 38.7× | ≥ 79.9× | ≥ 113× | ≥ 38.8× | ≥ 600× | ≥ 206× |
| tier-crossover-20260921 | 12688 points full | 61 | ≥ 13.0× | **≥ 105×** | ≥ 36.1× | ≥ 74.5× | — | — | ≥ 605× | — |
| tier-crossover-20260921 | 12688 points compact | 61 | ≥ 12.2× | **≥ 99.3×** | ≥ 34.0× | ≥ 70.2× | — | — | ≥ 573× | — |
| tier-crossover-20260921 | 12688 points folded | 61 | ≥ 13.3× | **≥ 108×** | ≥ 37.1× | ≥ 76.5× | ≥ 108× | ≥ 37.2× | ≥ 608× | ≥ 209× |
| tier-crossover-20260921 | 18544 points full | 61 | ≥ 15.4× | **≥ 125×** | ≥ 42.8× | ≥ 88.4× | — | — | ≥ 745× | — |
| tier-crossover-20260921 | 18544 points compact | 61 | ≥ 14.4× | **≥ 117×** | ≥ 40.0× | ≥ 82.5× | — | — | ≥ 695× | — |
| tier-crossover-20260921 | 18544 points folded | 61 | ≥ 13.6× | **≥ 110×** | ≥ 37.7× | ≥ 77.8× | ≥ 111× | ≥ 37.9× | ≥ 651× | ≥ 225× |

**The headline.**

| | `vs rho` |
|:--|--:|
| quoted, against single-target rho | 1.17× |
| re-read, against batch rho at `k = 32` | **6.51×** |
| rho's step priced at the canonical walk | 2.37× |
| rho's step priced at Bailey's step (a model) | 4.67× |
| the build priced per stored pair | 17.5× |
| both, canonical step | **6.38×** |
| both, Bailey's step (a model) | 12.6× |
| cold, one target | 43.7× |
| cold, with both | 42.8× |

With both corrections, the build becomes 75% of that pipeline rather
than 33%, and collection falls from 54% to 20%.

**Where the factor comes from.**  The thread's reference over batch rho at
`k = 32` is `5.56×` at `n = 41`, `8.63×` at `n = 53` and `8.10×` at
`n = 61`.  The batch accounts for `4.6–5.2×` of it: the single-target floor
over the `k = 32` measurement.  The rest comes from the thread's
single-target walk sitting `1.20–1.78×` above that floor.

**The two corrections pull in opposite directions, and on the headline row
they nearly cancel.**  Pricing rho's step lowers every ratio, by
`2.74–2.92×` at the canonical walk and `1.38–1.41×` at Bailey's step.
Pricing the build raises each row by as much as its build weighs.  The
n = 41 rows on the 15,744- and 16,400-point bases rise `1.7–2.7×`, the
5,248-point rows `2–4%`, and the n = 53 and n = 61 rows `0.1–2.4%`.  The
lowest reading anywhere in the table is `2.37×`, and it is the lopsided
one: the headline row with rho's step priced and the build still at one
unit a pair.  Priced on both sides, the closest figure is **`6.38×`**.

**Cold, one target costs from `43.7×` rho** (the headline row; `42.8×`
with both corrections) **to `≥ 1,548×`**.  The `k = 1` measurement it
divides by carries a 95% interval of `±28%` at `n = 41`, `±35%` at
`n = 53` and `±56%` at `n = 61`.

### 19.6 Classification

**Accounting.**  No index-calculus count changed.  The reference was
matched to the problem the figures solve, and the thread's own two
uncharged costs were priced.  The thread's engineering steps stand as
measured, because they are ratios of index-calculus totals on the same
instances and do not depend on the reference: `1.139×` for the packed
selection, `1.311×` for the least-mention aim, `7.26×` from the rung it
started from.  What falls is the distance to rho, which those steps were
read against.  The closest figure moves from `1.17×` to `6.5×` batch rho
as counted.  It moves to `6.4×` with every rho step and stored pair priced
by time, or to `12.6×` if Bailey's step is taken with the canonical walk's
count.  On the curves, bases and target counts these figures cover,
nothing measured end to end on an elliptic curve is below rho.

### 19.7 What stays open

- **Bailey's walk, counted.**  Its step costs half a canonical step.
  Until its step count is measured, the `4.67×` and `12.6×` readings are
  models.
- **The one-target references are thin.**  Sixteen (eight) runs put
  `±28–56%` on every cold figure.
- **The thread's single-target walk is `1.56–1.78×` its floor at
  `n = 53` and `61`.**  Any single-target figure priced against it there
  is flattered by as much.
- **The three-regime ledger's Koblitz rows** are single-target and in the
  affine unit.  Priced at the table-driven canonicalisation, their
  reference's step rises by `9–12%` (§19.4), so those rows would read
  `8–11%` lower.  They are not re-read here.
- **The build's two passes** are the thread's largest phase once priced,
  and the obvious next lever for it.

### 19.8 After the merge with main (2026-09-25)

Both prices in §19.4 were measured at `3e8dd352`.  Since then main has
made the Koblitz primitives faster (#689 and the AVX-512 kernels that
followed it), with bit-identical outputs:

- the least rotation searches the longest zero runs, about `log₂ n`
  steps instead of `n`;
- the folded build keys its pairs in bulk;
- the scan forms its rests without their ordinates.

Every frozen count therefore stands.  On the merged code, the walk
replays the frozen runs it was checked against exactly: the eight
instances of §18's ladder and the first batches of §19's `n = 41` run.
But on current main a canonical rho step, a stored pair and the
index-calculus scan all cost less than §19.4 measured.  The re-read is
internally consistent, because the frozen figures, their build price and
the reference's step all date from that commit.  It is not a price of
current main.  Re-pricing there means re-running both sides, which is
not done here.

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
