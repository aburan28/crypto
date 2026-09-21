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
