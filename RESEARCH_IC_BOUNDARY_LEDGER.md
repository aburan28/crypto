# The index-calculus boundary ledger: three regimes, one unit, every phase priced

**Harness:** `src/cryptanalysis/ic_boundary.rs` (pipelines, rho reference, ledger),
`src/cryptanalysis/ic_oracle_pricing.rs` (decomposition oracles per target),
`src/cryptanalysis/ic_corpus.rs` (Trimoska-style instance corpus)
**Command:** `cargo build --release --bin ic && ./target/release/ic boundary --oracles --out <file>`
**Frozen run:** `docs/ic/runs/ic-boundary-ledger-2026-09-21.json` (every number below is
read from it; the Markdown tables are the report's own `markdown` fields)
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

## 2. One table, one unit

<!-- FILL:MAIN_TABLE -->

## 3. Reading the table

<!-- FILL:READING -->

## 4. Every phase priced: the exponents

<!-- FILL:EXPONENTS -->

## 5. The decomposition oracles, per target

<!-- FILL:ORACLES -->

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
| F4 over `F_p` | `f4_fp::f4`, `f4_fp::solve` | degree-bounded F4 with dense modular elimination | not on a pipeline here; used by the `E(F_{p³})` Gaudry work (`RESEARCH_RESIDUAL_WALKS.md`) |
| XL, crossbred | `pq_xl`, `crossbred` | linearisation at a degree bound; guess-and-linearise | not priced here |
| CDCL SAT with native parity rows | `sat::Solver`, `semaev_sat` | Gauss–Jordan on XOR rows at the propagation fixpoint, Macaulay rows as implied clauses | oracle cells `cdcl_sat_native_xor` |
| meet in the middle (pair table) | `ic_boundary::PairTable` (this note), `koblitz_index_calculus::PairSumTable` (the production one, with compact and folded tiers) | all pair sums keyed by abscissa; one probe (`m = 2`) or `|F|` subtractions and probes (`m = 3`) | every regime |
| enumeration | `ic_oracle_pricing::enumerate_counted`, `koblitz_index_calculus::enumerate_decompose` | `|F|^{m−1}` group additions | oracle cells `enumerate` |
| dense Gauss–Jordan over `Z/rZ` | `ic_boundary::IncrementalGauss` | incremental reduced echelon form, multiply-subtracts counted | every regime (`linear_algebra` phase) |
| relation filtering + block Wiedemann | `koblitz_sparse_la` | singleton and clique removal, structured elimination, Krylov sequence with a matrix Berlekamp–Massey | not needed at these matrix sizes; measured at `n^{0.68}` on `E(F_{p³})` in `RESEARCH_RESIDUAL_WALKS.md` |
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
`n19l6-19-U` is satisfiable despite its name (`RESEARCH_TRIMOSKA_BENCHMARKS.md`),
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
