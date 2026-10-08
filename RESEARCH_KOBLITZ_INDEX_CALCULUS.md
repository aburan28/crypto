# Index calculus on Koblitz curves with Frobenius-invariant factor bases

**Modules:** `src/cryptanalysis/koblitz_index_calculus.rs`,
`src/cryptanalysis/koblitz_groebner.rs`, `src/cryptanalysis/koblitz_bench.rs`
**Scaling target:** `RESEARCH_KOBLITZ_SCALING_TARGET.md` — the scoped,
measurable follow-on thread (primary metric, correctness gate, five
pre-registered hypotheses)
**Demos:**   `cargo run --release --example koblitz_index_calculus_demo`,
`cargo run --release --example koblitz_semaev_groebner_demo`
**Paper:**  S. D. Galbraith, R. Granger, S.-P. Merz, C. Petit,
*On index calculus algorithms for subfield curves*, SAC 2020
(ePrint 2020/1315).

## The setting

A Koblitz curve `K_a : y² + xy = x³ + a x² + 1`, `a ∈ {0, 1}`, is
defined over `F_2` but used over `F_{2^n}`.  Because the coefficients
lie in the base field, the `2`-power Frobenius `π(x, y) = (x², y²)` is
an endomorphism of `E(F_{2^n})`; it satisfies `π² − tπ + 2 = 0` with
`t = (−1)^{1−a}`, and on the large prime-order subgroup `⟨G⟩` it acts
as multiplication by a scalar `λ` (a root of `λ² − tλ + 2 mod r`).

Pollard rho already exploits this: walking on `π`-orbits gives the
Gallant–Lambert–Vanstone `√(2n)` speed-up.  GGMP ask what the *same*
structure does for index calculus, and answer: `≈ n` in relation
collection and `≈ n²` in the linear algebra.

## Where the speed-up lives

Take the factor base `F` **Frobenius invariant**, `π(F) = F`.  Then `F`
splits into `π`-orbits of size dividing `n`, and for `P = π^k(Q)`

```
log_G P = λ^k · log_G Q   (mod r).
```

So the `|F|` unknowns collapse to `≈ |F| / n` orbit unknowns:

| Stage               | Classical | Frobenius-invariant | Saving |
|---------------------|-----------|---------------------|--------|
| relations needed    | `|F|`     | `|F| / n`           | `n`    |
| matrix (both dims)  | `|F|`     | `|F| / n`           | `n²`   |
| decomposition search| `m!` copies of each relation | ordered tuples only | `m!` |

The `m!` saving is the paper's "distribute the summands over
`F_1 = F, F_2 = π(F), …, F_m = π^{m−1}(F)`, then rewrite the relation so
that all `P'_i ∈ F_1`" step.  Because `π(F) = F`, the shifted bases
coincide; what survives is that each summand carries a known
`λ^{k_i}` weight into the column of its orbit representative.

## Building an invariant factor base (linearised polynomials)

The construction implemented here is the paper's first one:

1. `x^n − 1 = (x − 1) f_1 ⋯ f_s` over `F_2`, each `f_j` irreducible of
   degree `ℓ = ord_n(2)`.
2. `F_j(X) = Σ_k f_{j,k} X^{2^k}` — the linearised polynomial of `f_j`.
   Under the isomorphism `F_2[x] ≅ {linearised polys, composition}`,
   `x^n − 1 ↦ X^{2^n} − X`, so all `2^ℓ` roots of `F_j` lie in
   `F_{2^n}`; and since `f_{j,k} ∈ F_2`, the root set is closed under
   squaring.
3. `F := { P ∈ E(F_{2^n}) : F_j(x(P)) = 0 }`, of size `≈ 2^ℓ`.

The root space is computed as the kernel of the `F_2`-linear map
`v ↦ F_j(v)` written in the polynomial basis — self-checking, since the
kernel dimension must come out equal to `deg f_j`.

## Measured on the toys

`cargo run --release --example koblitz_index_calculus_demo`:

```
 n  a    ℓ   #E(F_2^n)          r    h     |F|   orbits
 7  1    3         142         71    2      15        3
 9  0    6         508        127    4      55        7
 9  1    6         518         37   14      73        9
11  1   10        1982        991    2     991       91
13  0   12        8012       2003    4    4005      309
```

For `K_0 / F_2^9`: `x^9 − 1` has the single non-trivial factor
`x^6 + x^3 + 1`, giving `F(X) = X + X^8 + X^64`, a 64-element invariant
subspace, `|F| = 55` points in `7` orbits, and `λ = 22` on the
order-127 subgroup.  Every random `(a, b)` decomposed on the first try;
`d` is recovered from an `8 × 8` system instead of a `56 × 56` one.

Speed-ups, measured as `|F| / #orbits`:

```
 n  a     |F|   orbits   relations   lin. alg.       rho
 9  0      55        7       7.86×      61.73×     4.24×
11  1     991       91      10.89×     118.59×     4.69×
13  0    4005      309      12.96×     167.99×     5.10×
```

which is the paper's `≈ n` / `≈ n²` against rho's `√(2n)` — and the
reason the paper's conclusion still reads *"security of Koblitz curves
remains strong"*: index calculus gains **more** from the subfield
structure than rho does, while remaining the worse attack outright at
the sizes anyone deploys.

## What is real here and what is not

Real:

- the Semaev/Weil-restriction decomposition oracle solved by matrix-F4,
  cross-checked target by target against exhaustive search;
- the invariant factor base from linearised polynomials, verified
  `π(F) = F` point by point;
- the orbit decomposition and the `λ^k` relation rewriting, verified
  against brute-forced discrete logs;
- the cofactor-corrected linear system `h·a + h·b·d ≡ Σ λ^{k_i} x_{o_i}`
  and its solution mod `r`, with the recovered `d` re-checked as
  `[d]G = Q` before it is returned.

Not real:

- **anything at deployed sizes.**  `n ≤ 41`: the factor base is
  materialised and `#E` is factored by trial division.
- **the Gröbner engine's constants** — see below.

## The decomposition oracle: Semaev, solved three ways

The question `is R = P_1 + … + P_m with all P_i ∈ F?` is answered three
ways, and the tests assert all three always agree:

- `DecompositionStrategy::Enumerate` — ordered-tuple search over the
  materialised factor base, `|F|^{m−1}` group operations per target.
- `DecompositionStrategy::Groebner` (default) — the real algorithm's
  oracle, in `koblitz_groebner.rs`.
- `DecompositionStrategy::Sat` — the same system as CNF, solved by the
  repo's CDCL solver (see "A third oracle" below).

The algebraic path, for `m = 2`:

1. `S₃(x₁, x₂, x(R)) = (x₁+x₂)²x(R)² + x₁x₂x(R) + (x₁x₂)² + b = 0`
   is the vanishing condition for `±P₁ ± P₂ ± R = O`.
2. Write `x_i = Σ_t u_{i,t} b_t` over an `F_2`-basis of the invariant
   subspace `V`.  Restricting the unknowns to `V` is what makes the
   system solvable — and it keeps the degree low: **quadratic** for
   `m = 2`, because
   squaring is `F_2`-linear in characteristic 2 (and in the Boolean
   quotient `p² = p`, so a symbolic square is a permutation of
   coordinates), while `x₁x₂` is bilinear.
3. Weil-restrict: one `F_{2^n}`-equation becomes `n` Boolean equations
   in `m·ℓ` unknowns.  For `m ≥ 3`, `S₃` is *chained* rather than
   resolved (`S_k` has degree `2^{k−2}` per variable) — `m` summands
   become `m−1` links over `m−2` intermediate field unknowns, and since
   those intermediates are unknown too the term `x₁x₂·e` makes the
   chained system **cubic** rather than quadratic.  That is why `m ≥ 3`
   costs so much more in every engine.
4. Solve by matrix-F4 with splitting: reduce the Macaulay matrix at
   degree `2, 3, …`; a reduction yielding the constant `1` refutes the
   branch; rows collapsed to `v_i (+1)` are propagated; where the
   algebra stalls, split on the lowest free variable.  Roots are lifted
   back to points — `S₃` fixes the summands only up to sign — and the
   group identity re-checked, so no spurious relation can escape.

For `K_0 / F_2^9`: 12 unknowns, 9 equations, degree 2.

```
  matrix-F4 at degree 2:    8 reduced rows in  33.02µs
  matrix-F4 at degree 3:  104 reduced rows in 311.54µs
  matrix-F4 at degree 4:  584 reduced rows in   4.03ms
```

### Two engines

`SolverEngine::Buchberger` runs the repo's textbook Gröbner engine
(`pq_groebner_f2`) at every node; `SolverEngine::MatrixF4` (default) runs
the Macaulay linear algebra this module adds.  Same verdicts, measured
on the same end-to-end DLP:

| engine | time to solve `Q = [53]G` on `K_0 / F_2^9` |
|--------|--------------------------------------------|
| matrix-F4 | 21.5 ms |
| Buchberger | 85.5 s |

Buchberger chokes on the dense 12-variable system (~8 s for a single
basis); F4's one linear-algebra pass over all products at bounded degree
is the fix, and it is also what the real algorithm uses.

### A third oracle: SAT

`DecompositionStrategy::Sat` hands the *same* system to the repo's CDCL
solver (`cryptanalysis::sat`) through
`semaev_sat::encode_boolean_system`: one Tseitin auxiliary per monomial
of degree ≥ 2, one parity constraint per equation, and a blocking clause
per rejected model so roots can be enumerated.  A final UNSAT is a
refutation in exactly the sense the F4 constant-`1` row is.

This is the Soos–Nohl–Castelluccia comparison the `semaev_sat` module was
written for, now runnable on identical inputs: clause learning and degree
growth blow up on different systems.  Measured over 24 targets:

| instance | search | SAT | verdict |
|---|---|---|---|
| `K_0/F_2^9`, m=2 | 0.58 ms | 79 ms | 24/24 agree |
| `K_1/F_2^9`, m=2 | 0.46 ms | 30 ms | 24/24 agree |
| `K_0/F_2^7`, m=2 | 0.10 ms | 4.5 ms | 24/24 refuted (UNSAT) |
| `K_0/F_2^9`, m=3 | 0.58 ms | 3.44 s | 24/24 agree |

Zero spurious models in every run — each model is re-checked against the
original equations before it is used, and the lifted points are
re-checked in the group.

The `m = 3` row is the interesting one: the chained system is cubic, and
SAT pays for that much as F4 does.

Fixing the encoder was part of this: the chain-auxiliary count for wide
parity constraints was `⌈(w+1)/2⌉ − 2`, which under-allocates from
`w = 6` on and panicked with "ran out of XOR-chain auxiliaries".  The
fold consumes three literals per auxiliary, so the count is `⌊w/3⌋`;
`xor_parity_encoding_is_exact_at_every_width` pins it by counting models
at every width up to 8.  That bug was in the merged `semaev_sat`, not in
new code.

### Where the algebra wins, and where it does not

It refutes.  On `K_0 / F_2^7` the invariant subspace has 8 elements but
only one is the abscissa of a curve point, so **no** target decomposes.
All 28 targets are closed by an F4 infeasibility certificate, with zero
branches searched — an answer exhaustive search cannot give in kind, only
by exhausting the space.

It does not, yet, win on wall-clock:

```
 [k]G        search     matrix-F4     verdict
    3       19.10µs        1.86ms   both: yes
   17       42.12µs        1.53ms   both: yes
   53       17.64µs        2.08ms   both: yes
  101       62.04µs        2.21ms   both: yes
```

At `|F| = 55` the search is `55` point additions and the Macaulay matrix
already has 2325 columns.  The crossover is past what this module can
materialise a factor base for (`|F| = 2^ℓ` points in memory, `ℓ = ord_n(2)`),
so the honest statement is: the algebraic oracle is implemented, correct,
and cross-checked, and it is the one that scales — cost governed by the
degree of regularity rather than by `|F|` — but at `n ≤ 24` the brute
oracle is faster and remains available.

## Which invariant factor bases exist at all

A Frobenius-stable `F_2`-subspace of `F_{2^n}` is an
`F_2[x]/(x^n − 1)`-submodule — a binary cyclic code of length `n` — so
it is a divisor of `x^n − 1`, and the irreducible factors correspond one
for one with the 2-cyclotomic cosets mod `n`.  That is a **complete
classification**: the achievable dimensions are exactly the subset sums
of the coset sizes (`available_subspace_dimensions`).

```
   n  cosets  coset sizes                 available dimensions
   7       3  [1, 3, 3]                   0,1,3,4,6,7
   9       3  [1, 2, 6]                   0,1,2,3,6,7,8,9
  15       5  [1, 2, 4, 4, 4]             0,1,2,3,4,5,6,7,8, …
  31       7  [1, 5, 5, 5, 5, 5, 5]       0,1,5,6,10,11,15,16,20, …
 127      19  [1, 7 × 18]                 0,1,7,8,14,15,21,22,28, …
 131       2  [1, 130]                    0,1,130,131
 163       2  [1, 162]                    0,1,162,163
```

The last two rows are the point: where `2` is primitive mod `n` there
are only two cosets, so the construction is **empty exactly at the sizes
that matter**.  That was asserted in this repo before; it is now
computed and tested.

### Divisor bases: sizing the base to the instance

`build_frobenius_factor_base` uses one irreducible factor, so its
dimension is stuck at `ord_n(2)`.
`build_frobenius_factor_base_from_divisor` takes any *product* of
factors, so `|F| ≈ 2^dim` is tunable across the whole classification.

This matters more than it sounds.  The summand count a decomposition
needs is `m ≈ n/dim`, and `m ≥ 3` is what forces the chained system and
its `(m − 2)·n` extra unknowns.  A big enough invariant subspace buys
`m = 2` — **no chaining, a quadratic system, `≈ 2·dim ≈ n` unknowns**.
Measured on toy curves:

| curve | dim | \|F\| | orbits | m | vars | search | F4 | SAT |
|:------|----:|------:|-------:|--:|-----:|-------:|---:|----:|
| K_1/2^7 | 3 | 15 | 3 | 2 | 6 | 11 µs | 198 µs | 144 µs |
| K_1/2^7 | 6 | 71 | 11 | 2 | 12 | 4.7 µs | 1.4 ms | 448 µs |
| K_0/2^9 | 8 | 253 | 29 | 2 | 16 | 8.1 µs | 3.9 ms | 5.3 ms |
| K_1/2^15 | 10 | 1057 | 73 | 2 | 20 | 64 µs | 17 ms | 32 ms |

### The cofactor class decides which `m` can work

A separate constraint, and not a size one.  `x = 0` lies in every
invariant subspace, so the 2-torsion point is always in `F`; more
generally the base can sit entirely off `⟨G⟩`.  A target `R ∈ ⟨G⟩` has
`[r]R = O`, so a decomposition needs the summands' `h`-torsion classes to
cancel.  On `K_1 / F_2^7` (cofactor 2) **no** factor-base point is in
`⟨G⟩`, so odd `m` decomposes *nothing* — at any `|F|`:

```
dim 3, |F| = 15   m=2: 11/11   m=3: 0/11   m=4: 11/11
dim 6, |F| = 71   m=2: 11/11   m=3: 0/11   m=4: 11/11
```

`FrobeniusFactorBase::admissible_summand_counts` computes this up front
from one scalar multiplication per point, instead of searching for
decompositions that cannot exist.  It was found by a test written to
assert the opposite.

## Where this goes next

`RESEARCH_KOBLITZ_SCALING_TARGET.md` turns the open end of this work
into a measurable thread.  The short version: every oracle here works at
`m = 2` and `m = 3`, a useful attack needs `m ≈ n/ℓ`, and what stops us
is not the algebra but the chaining representation —

```
    unknowns(n, ℓ, m) = m·ℓ + (m − 2)·n
```

which at `n = 63, ℓ = 6` wants 633 unknowns against a 64-variable
budget.  The first fall degree, measured over 16 target draws per
instance across the whole ladder to `n = 63`, never exceeds 3 — so the
systems themselves stay benign; it is the variable count that bites.

## End-to-end optimisation — 2026-09-10

**Modules:** `koblitz_factor_base_search.rs` (search),
`koblitz_relation_solver.rs` (incremental linear algebra), additions to
`koblitz_index_calculus.rs` (pair table, orbit restriction, symmetry
breaking) and `sat.rs` (`Solver::add_vars`).
**Tool:** `ic search`, `ic run --factor-base`, `ic run --solver pair-table`.
**Docs:** `docs/ic/README.md`.

Four changes to the pipeline as `ic` runs it, each cross-checked against
the reference oracle it replaces or the dense solver it supersedes.

### 1. Choosing the factor base by what the pipeline pays for

The run cost is `trials × (cost per trial) + linear algebra`, and the
expected trial count is `(U + 1 + extra) / p_m` where `U` is the number of
relation columns and `p_m` the fraction of subgroup points that decompose
into `m` base points.  Neither `|F|` nor the dimension is the objective;
`U / p_m` is, and it was never measured by the selection tools — `ic
compare` timed the legacy single-factor family and the census examples
counted coverage per *orbit* without a column model.

`koblitz_factor_base_search` measures `p_m` **exactly**: a pair-sum table
of `F` (every `P_i + P_j`, sorted by a packed point key) enumerates every
witness `R = P_{i_1} + … + P_{i_m}` of every target, over the whole
subgroup when `r − 1 ≤ 4096` and over a seeded sample otherwise.  Cofactor
classes need no separate treatment — a sum lands in `⟨G⟩` or it does not.
Candidates are the complete divisor lattice of `x^n − 1` in a dimension
window, Frobenius unions of random seeds, and the 2-torsion saturation of
either; each is then **pruned** greedily, dropping a signed orbit while
that lowers `(U + 1 + extra)/p_m`.  The witness list makes each pruning
step an exact recount, not a re-search: removing orbit `o` loses exactly
the targets all of whose witnesses use `o`.

The result is a serialisable recipe (`FactorBaseSpec`) that rebuilds the
identical base, so `ic run --factor-base` replays it and `ic search`
validates the best few by real child runs on fresh known-answer fixtures
before selecting the fastest one that verified every holdout.

What it finds.  On `K_1 / F_2^15` (r = 211, h = 154) the legacy base —
factor index 0, 31 points — yields **no relation at all** in 20 000
trials: `ic run --degree 15 --curve-a 1` was simply unsolvable.  The
search scores 86 candidates on all 210 targets in 8 s and validates a
30-point, one-column pruned union that covers every target; the
end-to-end run then takes 21 ms.  The prior review's numbers reproduce
exactly (`the_prior_review_numbers_reproduce`: divisor `[0, 2]` 61 points
150/210 at `m = 3`, union seeded by 3468, 4413 210/210).

Pruning is where the algorithmic content is.  A signed orbit is a column
the linear algebra has to determine and a set of summands that make
targets reachable; the trade is explicit in `(U + 1 + extra)/p_m`, and at
`n = 15` the greedy pass takes a 301-point, five-column saturated union
to 30 points and one column at unchanged coverage — the entire subgroup
is a sum of two points of *one* signed orbit of the projected subgroup.
That is the extreme case; typical steps drop the singleton 2-torsion
orbit (`x = 0`, which every invariant subspace contains and which never
helps) and the low-participation orbits of a saturation.

### 2. A meet-in-the-middle oracle

`DecompositionStrategy::PairTable` answers "is `R` a sum of `m` base
points" from the same table: one lookup at `m = 2`, `|F|` at `m = 3`,
`|F|²` at `m = 4`, against `|F|^{m−1}` group operations for enumeration.
It is exact and complete like enumeration and re-adds every answer in
the group before it becomes a relation; it costs `|F|(|F|+1)/2` entries
of 16 bytes once per run (the packed key separates `P` from `−P` by
comparing `y` with `x + y`, exact to `n ≤ 62`).  At `m = 3` this is what
makes the raised degree cap usable: `n = 31` needs `m = 3` on a
dimension-11 base, and enumeration would pay `|F|² ≈ 4 · 10^6` group
operations per target.

### 3. Incremental modular linear algebra

The driver kept every relation and, in early-stop mode, re-ran a dense
`BigUint` Gaussian elimination after each batch — `O(U³)` big-integer
work per check that could not tell "not yet determined" from
"determined" and had to verify a partial solution in the group.
`IncrementalRelationSolver` keeps the matrix in reduced row echelon form
over `Z/rZ` with `u64` arithmetic and updates it per row in `O(U²)`; a
pivot in the `d` column *is* the answer, so the run stops the moment the
scalar is pinned.  Dependent relations are counted rather than padded,
an inconsistent relation fails the run closed, and `verification_failures`
would flag a pinned-but-wrong scalar, which only a wrong relation could
produce (it has not happened).  Cross-checked on random planted systems
against `gaussian_eliminate_mod_n` at five moduli up to 32 bits.

Two accounting changes ride along in `ic run`'s default mode: columns
are signed Frobenius orbits merged by cofactor projection (both public
identities were already implemented but disabled in `ic`), and relation
batches are decomposed in parallel for every oracle, not only SAT.
`--control` restores the previous accounting for matched comparisons.

### 4. What did not help: symmetry breaking in SAT

The SAT oracle refutes each root `m!` times because the summands are
interchangeable.  `SatDecompositionOptions::symmetry_breaking` installs
exact lexicographic `code(x_1) ≤ … ≤ code(x_m)` constraints
(`Solver::add_vars` grows the encoding after the fact); the constraint
admits exactly the sorted assignments (`lex_leader_constraint_admits_exactly_the_sorted_pairs`
counts them) and changes no verdict.  It also does not pay:

| instance | targets | conflicts off | conflicts on |
|:---------|--------:|--------------:|-------------:|
| `K_0/2^9`, m = 3 | 8 | 5 512 | 6 412 |
| `K_1/2^9`, m = 2 | 16 | 428 | 823 |
| `K_0/2^7`, m = 2 | 16 | 70 | 79 |
| `K_1/2^15`, m = 3 | 8 | 0 | 0 |

The degree-2 Macaulay rows and the trace row already leave these
systems nearly propagation-closed, and the auxiliaries add decisions
without removing search.  Off by default, kept as a control
(`koblitz_sat_symmetry_ablation`).

### Measured, before and after

`ic run` on the same fixtures, wall time of the whole process, this
branch against `main` built from the same tree; 4 cores.  Columns are
relation unknowns; the new default merges signed orbits by cofactor
projection.

| curve | solver | columns main → branch | relations main → branch | wall main | wall branch |
|:------|:-------|:---------------------:|:-----------------------:|----------:|------------:|
| `K_0/2^9` | enumerate | 7 → 3 | 11 → 4 | 0.010 s | 0.009 s |
| `K_0/2^9` | groebner | 7 → 3 | 11 → 4 | 0.026 s | 0.012 s |
| `K_0/2^9` | sat | 7 → 3 | 11 → 4 | 0.014 s | 0.009 s |
| `K_1/2^11` | enumerate | 91 → 45 | 95 → 16 | 0.035 s | 0.028 s |
| `K_1/2^11` | groebner | 91 → 45 | 95 → 16 | 0.84 s | 0.071 s |
| `K_1/2^11` | sat | 91 → 45 | 95 → 32 | 0.216 s | 0.063 s |
| `K_0/2^13` | enumerate | 309 → 77 | 313 → 12 | 0.140 s | 0.119 s |
| `K_0/2^13` | groebner | 309 → 77 | 313 → 12 | 6.78 s | 0.223 s |
| `K_0/2^13` | sat | 309 → 77 | 313 → 20 | 2.17 s | 0.249 s |
| `K_1/2^17` | enumerate | 13 → 6 | 17 → 7 | 0.063 s | 0.032 s |
| `K_1/2^17` | groebner | 13 → 6 | 17 → 7 | 1.33 s | 0.222 s |
| `K_1/2^17` | sat | 13 → 6 | 17 → 7 | 0.555 s | 0.095 s |
| `K_1/2^23` | enumerate | 91 → 45 | 95 → 23 | 4.31 s | 0.495 s |
| `K_1/2^23` | groebner | 91 → 45 | 95 → 27 | 206 s | 20.6 s |
| `K_1/2^23` | sat | 91 → 45 | 95 → 23 | > 900 s (killed) | 137 s |
| `K_1/2^23` | pair-table (new) | — | 23 | — | 6.2 s |

Three things the table says.  The relation count is what moved: with
signed, projected columns and the exact stopping rule the run needs
16 relations at `n = 11` where it collected 95, and 12 at `n = 13` where
it collected 313 — the algebraic oracles gain 10–30× because they pay per
relation.  Enumeration was already cheap per trial, so it gains only
what the trial count gives (9× at `n = 23`).  And the pair table does
**not** win at `m = 2` on these instances: its one-time build is
`|F|(|F|+1)/2` point additions at ≈ 9 µs each on this curve arithmetic
(4.4 s at `|F| = 991`, 12 s at `|F| = 4005`), against a dozen trials of
`|F|` additions for enumeration.  It pays when trials are many or
`m ≥ 3`, where enumeration is `|F|²` per trial and the table is `|F|`
lookups — which is exactly the regime past the old cap below.  On the
cofactor-2 curves `K_1/2^11`, `K_1/2^17` and `K_1/2^23` three summands
are inadmissible (no cofactor-class cancellation), and every solver
reports that up front instead of searching.

### Past the old cap

`MAX_N` is 40 (was 24); `KoblitzCurve::new` uses the sparse irreducible
search, which a test pins equal to the exhaustive one for every `n ≤ 24`,
so no existing fixture changes.  Above 23 only four Koblitz curves have a
prime-order subgroup larger than their cofactor — `K_1/2^29`
(r = 42 457), `K_0/2^31` (r = 1 439 393, h = 1492), `K_0/2^37` and
`K_0/2^39` — and the constructor rejects the rest, as it did before.

On `K_0/2^31` (h = 1492) the single-factor family cannot be used at all
in `ic` (`ord_31(2) = 5` gives 63–65-point bases whose three-summand
coverage is 0 on every one of them), so this is the first instance the
search is *necessary* for rather than merely better.  `ic search
--degree 31 --curve-a 0 --summands 3 --family divisor --max-dimension 11
--targets 256 --validate-top 2 --holdout 2` scores the 42 divisor bases
of dimension 5, 6, 10 and 11 and their pruned variants — 72 candidates
on 256 sampled targets in 13 minutes, three quarters of it the
`|F|²` pair tables — then validates the top two by child runs:

| rank | base | points | columns | coverage | expected trials | validated |
|-----:|:-----|-------:|--------:|---------:|----------------:|:----------|
| 1 | divisor `[0, 1, 6]`, pruned to 35 orbits | 2 170 | 35 | 0.605 | 62.8 | 2/2 |
| 2 | divisor `[0, 1, 6]` | 2 421 | 39 | 0.648 | 64.8 | 2/2 |

The pruning step removed four orbits and 251 points at a cost of 4
points of coverage, which the column count more than repays.  The
selected recipe then solves fresh targets end to end:

| target | trials | relations (independent) | pair table | relations | linear algebra | wall |
|:-------|-------:|------------------------:|-----------:|----------:|---------------:|-----:|
| `[654009]G` | 56 | 33 (33) | 8.65 s | 0.69 s | 0.6 ms | 10.4 s |
| `[1153191]G` | 68 | 36 (34) | 9.07 s | 0.84 s | 0.7 ms | 10.9 s |

Relation collection is under a second for 2 170 points and 35 columns;
the run is the pair table.  The first validation runs took 68 s each,
and profiling that gap found the cofactor-admissibility walk paying
`|layer| · |classes| ≈ 2.2 · 10^6` point additions at `h = 1492`; the
walk now generates each sumset layer from one representative per signed
Frobenius orbit and closes it by squaring (the layers are `π`- and
negation-invariant because the classes are), which is what took the run
to 10 s.  The naive walk is kept as a test reference
(`orbit_generated_admissibility_walk_matches_the_naive_one`).  This is
a 31-bit prime-order subgroup solved by index calculus over a
materialised base; it is still a toy, and Pollard rho would take
milliseconds on it — the point is that the pipeline now reaches the
sizes the scaling note's step 0 asked for, with the factor base chosen
by measurement rather than by hand.

## The CADO-style split: factor-base logs, then descent — 2026-09-10

**Modules:** `koblitz_index_calculus::{solve_factor_base_logs,
individual_log, FactorBaseLogTable}`.
**Tool:** `ic logs` (precompute), `ic solve` (per-target descent).
**Docs:** `docs/ic/README.md`.

The driver above bakes the target `Q` into every relation
(`R = [a]G + [b]Q`) and rebuilds the whole relation matrix for each `Q`.
A number-field-sieve pipeline instead solves the factor-base logarithms
**once per curve** and then recovers each target with a single
individual-logarithm relation.  Both stages are now implemented.

- **Precompute.**  `solve_factor_base_logs` draws `R = [a]G` probes
  (`b = 0`), decomposes each over the factor base, and rewrites it as a
  row `Σ_o c_o x_o ≡ h·a (mod r)` over the projected columns.  Once the
  rows determine every column, the whole logarithm vector is read off in
  one dense modular solve.  Only the projected representation is used:
  its columns are the canonical cofactor projections `R_o ∈ ⟨G⟩`, so a
  column logarithm `x_o = log_G R_o` is a genuine discrete log that
  certifies itself — the table is accepted only when `[x_o]G == R_o`
  for **every** column, which depends on nothing but the curve and the
  table.
- **Descend.**  `individual_log` draws `R = [a]G + [b]Q` until one
  decomposes, giving `h·a + h·b·d ≡ Σ_o c_o x_o (mod r)`; with the
  column logs known, `d = log_G Q` is one modular inverse, re-checked as
  `[d]G == Q`.  No relation matrix — one decomposition and a lookup.

Measured on `K_0/2^31` (r = 1 439 393) over the search-selected
dimension-11 base (35 projected columns, `m = 3`):

| stage | work | wall |
|:------|:-----|-----:|
| `ic logs` precompute (once) | 129 probes to full column rank, all 35 logs certified | 15.4 s |
| `ic solve` descent (per target) | 1–2 relations | ~10 s* |

`*` each `solve` is a separate process that rebuilds the pair table
(≈9 s at this base size); in a single process the table is built once
and the descent itself is milliseconds.  The point is structural: the
per-target cost is one decomposition, not a 35-column matrix rebuild, so
`n` targets cost one precompute plus `n` cheap descents rather than `n`
full solves.  This is the first stage that makes precomputation reuse
and batched targets meaningful, and the shortest path to the individual
logarithm / descent stage a full pipeline needs.  Tested by
`factor_base_logs_precompute_and_descend_every_target` (cross-checked
against brute-forced logs) and `ic`'s
`factor_base_logarithm_database_precomputes_then_descends`.

## Fast binary-field hot path — 2026-09-23

**Module:** `koblitz_fast_arith.rs` (single-word field + batched affine
sums); hooks in `PairSumTable::{build, witnesses}` and trial sampling in
`koblitz_index_calculus_dlp_observed`.
**Tool:** unchanged (`ic run --solver pair-table`).

Three changes, all bit-exact with the general arithmetic (new tests pin
`add`, `double`, `scalar_mul`, batch-vs-single, and `pack_words ==
pack_point`; every decomposition is still re-added on the slow path
before it becomes a relation, and the final `[d]G == Q` check is
untouched):

1. **Single-word field.** `F2mElement` heap-allocates a `Vec<u64>` per
   operation even though every index-calculus curve has `n ≤ 63`.  The
   hot path now uses the carry-less (`pclmulqdq`/`pmull`) + table-reduction
   field from `semaev_decomp::Gf2`, with the identical affine group law.
2. **One inversion per batch (Montgomery).** Pair-table build was one
   Itoh-Tsujii inversion per pair sum; it is now one inversion per table
   row.  The `m = 3` witness step (`R − P_k` for all `k`) is one inversion
   per target.  Doubling/negation edge cases take the slow path one by one.
3. **Fast trial sampling.** `R = [a]G + [b]Q` per attempt runs double-and-add
   on words.  Same RNG stream, same targets, same relations — only faster.

Same machine, same commands, `ic run --solver pair-table --summands 2`:

| curve | stage | before | after |
|:------|:------|-------:|------:|
| `K_0/2^13` (\|F\| = 4005, 8.0M sums) | pair-table build | 2.93 s | 0.21–0.34 s (~11×) |
| `K_1/2^23` (\|F\| = 2071, 2.1M sums) | whole process (CPU) | 14.1 s | 0.20 s (~70×) |
| `K_1/2^23` | relation sampling (56 trials) | 76 ms | ~2 ms (~30×) |

`m = 3` at `K_0/2^13` verifies in ~0.13 s table + 1 ms collection;
`ic search --degree 15` scores 142 candidates on the full 210-target
subgroup in 6 s.  Gates green: 99 `koblitz_*` lib tests, 14
`ic_framework`, 4 `ic_progress`, 5 new `koblitz_fast_arith` — zero
disagreements, enumerate/groebner oracles still agree target by target.

What is still slow: the Gröbner/SAT oracles (algebra untouched), `m = 4`
per-target cost (fast singles, not batched), and `F2mElement` everywhere
outside this hot path.

## Fast S3 scan kernel + degenerate-pair soundness fix — 2026-09-23

**Module:** `s3_x_roots` in `koblitz_fast_arith.rs`; call sites in
`examples/koblitz_s5_sat_instance.rs` (`scan_exceptional_relative_states`,
`scan_regular_relative_states`, stats/equivariance, lazy-theory callbacks,
positive-installer twists).
**Ledger link:** priorities #1/#2 — the edge-free scan that must replace
explicit pair-table setup.

The relative-Frobenius scans decide `S_3(left, right, ·)` solvability on
`R²·n` states (7.49M at n=53) with general-field arithmetic (~6µs/state:
Itoh-Tsujii inversion + half-trace + two S3 asserts, every op allocating).
The single-word port costs one fast inversion per state (~100–200ns).
Same-machine `relative_pair_stats` backend, identical invocations:

| rung | states | scan before | scan after | exceptional before → after |
|:-----|-------:|------------:|-----------:|:---------------------------|
| n=23 (4 reps) | 368 | 135 ms | 4.4 ms (31×) | 4 → 0 |
| n=41 (73 reps) | 218,489 | 1496 ms | 39.6 ms (38×) | 73 → 0 |
| n=53 planted edge-free | 7,492,928 | 46,393 ms | 1,247 ms (37×) | 376 → 0 |

The exceptional column is the second finding: **every pre-fix exceptional
state at these rungs was a false positive.** The old kernel gave up
whenever its combined inversion hit zero — `left == right` (diagonal) or a
zero side — but `S_3` degenerates to linear there with a real root
(diagonal `x_3 = x² + b/x²`; zero side `x_3² = b/(left+right)²`; only
`(0, 0)` is genuinely rootless). Each false exceptional became a wrong
nogood forbidding valid doubling/2-torsion pairs — e.g. the "76 nogood
clauses" at n=19 were exactly the diagonal states. The fast kernel returns
these single roots doubled (`[r, r]`, which the positive installer's
`first == second` arm pins correctly).

Gates: lib tests prove full root-set equality on all 16,384 pairs of
`F_{2^7}`, slow-kernel equality on generic pairs at n ≤ 53, and S3
vanishing independently; every stats run keeps 0 equivariance
discrepancies. End to end, the n53 planted edge-free run reproduces
byte-for-byte (same planted indices, SAT, 1 group-valid relation, 0 bad
lifts, **identical 3631 conflicts**) with scan 46.4 s → 1.2 s and 39,856
false nogoods gone. Pair-support certificates bumped to schema 2.0 so
pre-fix certs fail closed instead of installing stale nogoods. No ledger
promotion claimed: digests move (fewer exceptionals), unrestricted arms
stay UNKNOWN.

## Fast factor-base construction + Itoh-Tsujii inversion — 2026-09-23

Two changes in the shared library (`Gf2::inv`, `finish_factor_base_domain`):

1. **Itoh-Tsujii `Gf2::inv`.** The naive Fermat loop spent `n − 2`
   multiplications per inversion; the addition chain spends
   `⌊log₂(n−1)⌋ + popcount(n−1) − 1` (7 vs 51 at n = 53) for the same
   power, so outputs are bit-identical. Pinned by exhaustive inverse-law
   at n = 7 plus independent square-and-multiply Fermat and `F2mElement`
   cross-checks at n ≤ 63. Op-count win only: per-op microbenchmarks sit
   below this host's noise floor (±2× round variance), so no wall-clock
   number is claimed for it alone — it multiplies through every
   single-add path (sampling, `m = 4` walk, label checks).
2. **`finish_factor_base_domain` on words.** Lifts via the lib-tested
   fast `points_with_x`, orbit walks with word Frobenius/negation over a
   `(u64, u64)` index (exactly as injective as `point_key` at n ≤ 63, no
   packing limit), points converted back once. Same sets, same orders,
   same labels.

A/B on 27 divisor bases at n = 31 (stashed-baseline binary, same machine):
0.37–0.42 s → 0.03 s (~13×) with identical point counts; `ic search
--degree 15` replays deterministically (identical validated selection).
Gates: 105 `koblitz_*` lib tests, 14 `ic_framework`, 4 `ic_progress`,
S5 stats digests unchanged, n13 planted SAT at identical 13 conflicts.

## Batched enumeration oracle — 2026-09-23

`enumerate_decompose` — the shared reference oracle (surveys, search
census filter, cross-oracle tests, `ic run --solver enumerate`) — ran the
textbook recursion with one slow group operation per step. It now
projects the base to words once, answers the `m == 1` base case from a
word-keyed index, and batches each recursion level (`target + (−P_i)`
for all `i`) through one Montgomery inversion. Visit order and first
witness are unchanged; the slow recursion stays as the cross-check.

Microbench, 56 deterministic targets at `K_1/2^23`, `m = 2`
(stashed-baseline binary, same machine): 480–482 ms → 5.1–5.3 ms
(**~93×**, identical 23/56 found). End-to-end `ic run --solver
enumerate` at n = 23: 216 ms → 81–95 ms relation collection (the rest is
run-loop bookkeeping shared by all solvers, now the dominant share).
Gates: new `fast_enumeration_matches_slow_recursion_exactly` test pins
witness-for-witness equality with the slow recursion on subgroup,
off-subgroup, and Infinity targets at `m = 2, 3`; all cross-oracle
suites green with zero disagreements.

## Dispatch integration, batched sampling, logs rank gate — 2026-09-23

Three integration gaps closed around the fast arithmetic above (all
bit-exact, so RNG streams, trial counts, and tables do not move):

1. **Dispatch through fast enumerate.** Both run-loop dispatches
   (`koblitz_index_calculus_dlp_observed`, `decompose_once`) called the
   slow textbook recursion directly, bypassing `enumerate_decompose`.
   Routing them through it: end-to-end `ic run --solver enumerate` at
   n = 23 goes 216 ms → ~2.5 ms relation collection (~80×; the remainder
   is shared run-loop scaffolding).
2. **Batched multi-scalar sampling** (`batch_double` +
   `batch_scalar_mul` in `koblitz_fast_arith.rs`). MSB-first
   double-and-add with one batched doubling and one batched addition per
   bit position — ~2 inversions per bit for the whole batch instead of
   ~1.5 per pair per bit. Pinned bit-exact against `scalar_mul`
   (mixed bases incl. infinity, mixed lengths incl. zero) and wired
   into trial sampling with the RNG draw order preserved; the
   loop-carried CADO probes/descent draws use fast singles.
3. **Logs-precompute rank gate** (`U64RankTracker` in
   `koblitz_relation_solver.rs`). The precompute re-ran a dense BigUint
   Gaussian elimination after nearly every new row, but a solve can only
   succeed at full column rank; tracking rank in `O(cols²)` u64 ops per
   row skips every provably-doomed attempt (skipped systems are
   underdetermined against random true logs, so could never verify).

`ic logs` at n = 23, pair-table (stashed-baseline binary, same machine):
1.8–2.8 s → 0.065–0.076 s (**~25×**) with identical trials (271),
relations (106), and a verifying `ic solve`. Gates: new tracker tests
(counts, cross-check vs incremental solver on identical streams) plus
the fixed-monic-pivot bug they caught (elimination by non-monic pivots
overcounts rank — found by test, fixed by normalising on insert);
115 lib + 14 + 4 suites green; S5 stats unchanged.

## Parallel search scoring — 2026-09-23

`search_with_progress` scored candidate factor bases serially while the
pair-table build inside each score already saturated all cores in
bursts: n31 divisor search (140 candidates × 64 targets) took 40 s wall
at ~130% cpu. Specs are now scored on the process-wide rayon pool;
scoring is pure given the fixed spec/target set (RNG lives only in the
serial spec/target generation), so results are collected in spec order
with progress replayed serially — rankings, selections, and report
ordering are identical to serial scoring, only wall time and `*_ms`
fields move. Nested parallelism is safe (table builds split further on
the same pool).

Same machine, same commands: n31 divisor census 39.8 s → 5.9 s wall at
~960% cpu (**6.7×**) with byte-identical top candidates (family, points,
columns, coverage, expected trials all equal); n15 full search with
validation 6.0 s → 0.9 s with the identical validated selection.
Cost: transient peak RSS 1.68 GB at n31 (concurrent big tables, freed
per spec; disclosed in the report's resource counters) vs ~0.15 GB
serial. If that ever bites, the fallback is wave-chunked scoring at a
small throughput cost — not implemented, no evidence it is needed.
Gates: 115 lib + 14 + 4 suites green (search selection tests pin the
deterministic rankings).

## Parallel validation and comparison children — 2026-09-23

`ic search` validation and `ic compare` spawned one child `ic run`
process per (candidate, fixture) pair serially — up to 136 children for
compare, each a full pipeline run. New `run_children_ordered` helper in
`src/bin/ic/experiment.rs`: scoped-thread waves with a bounded
concurrency cap (`cpus/3` clamped to 2..=6 — each child carries its own
rayon pool and transient 100MB+ tables), per-job recipe files (the old
shared rank-tagged path would have collided), results collected in job
order with the first error surfacing as serial `?` would. Rankings,
winner selection (strict min in rank order), eligibility, and medians
are computed serially from ordered reports, so JSON is deterministic;
only wall time, per-child timings, and display interleaving move.

Measured: n31 search + 2×2 validation 9.8 s → 8.0 s (validation slice
~2.4 s → ~0.7 s); n23 compare (8 children) 0.16 s wall; groebner n17
search validation overlaps 4 slow children at ~845% cpu with correct
inconclusive eligibility. Value scales with child cost (groebner/sat at
larger n). Degenerate inputs (0 samples/holdout) now yield inconclusive
instead of panicking in `median`. Gates: framework compare/search tests
unchanged and green; no temp-file leaks (per-job Drop verified).

## Batched cofactor projection and target sets — 2026-09-23

Two remaining slow-scalar-mul loops, both replaced by the tested
`batch_scalar_mul` with identical outputs:

1. **`projected_signed_orbit_map` cofactor projection.** Every
   factor-base point was projected by `[h]` one slow inversion-heavy
   scalar mul at a time — 305 ms for a 2049-point base at n = 31, paid
   per candidate in every search (the dominant census cost). One batch
   call now projects the whole base: **305 ms → 4 ms (~77×)**, same
   projected set. n31 divisor census CPU halves (56.6 s → 28.6 s) with
   byte-identical top candidates.
2. **`TargetSet::new` point construction.** All subgroup (or sampled)
   target points in one batch: n13 exhaustive 2002 targets 72–85 ms →
   3 ms (~25×); n23 sampled 256 targets 51–56 ms → 0.8 ms (~65×).
   Searches replay deterministically (identical validated selection).

The point helpers `to_fast_point`/`from_fast_point` are now `pub(crate)`
for the search module's use. Gates: 115 lib + 14 + 4 suites green.

## Word-based admissibility: classes and layer walk — 2026-09-23

The factor-base stage hid a 0.2 s admissibility check at n = 31 (the
whole `ic run` was 0.25 s): `distinct_cofactor_classes` spent one slow
`[r]` scalar mul per signed orbit, and `m_can_decompose` walked
reps × classes sumset layers with one slow addition per pair (~30K at
`h = 1492`). Both now run on words — one batched `[r]` over the orbit
reps, one Montgomery batch per layer grid replayed in identical nested
order, Frobenius/negation/packing via word ops — with slow fallbacks
for curves too wide for one word. Sets, keys, traversal, and verdicts
match exactly (packing is proven equal; class order was already
HashMap-nondeterministic and only the boolean is consumed).

Isolated probe on the n31 recipe base (stashed-baseline binary, same
machine): classes 11–15 ms → 0.5–0.8 ms (~20×); `m_can_decompose`
197–209 ms → 1.3–1.6 ms (**~140×**), same 1054 classes and verdict.
End-to-end n31 run: 0.25 s → 0.05 s. Gates: the naive-walk agreement
test and the pointwise class-set test pin the rewrite; 115 + 14 + 4
suites green.

## Word-level projection-map walks — 2026-09-23

A deliberately small one: with cofactor projection batched, the
remaining `projected_signed_orbit_map` cost was its orbit-merging walks
(slow squares, negations, `BigUint` keys per step). They now run on
words with `(x + 1, y)` pair keys, which order *identically* to
`point_key` (numerical order on values `< 2^64`, `(0, 0)` still
reserved), so canonical choices, sort order, and locations match
provably — plus a slow-path fallback and a new test asserting
representative-for-representative, location-for-location equality
against it. Measured ~1.3× on the map stage (4.0 ms → 2.8 ms at
2049 points); small absolute, but it removes the last slow-field
arithmetic from the projection path. Column order is unaffected for
consumers (relation columns are self-consistent per run; log databases
match columns by point identity). 116 + 14 + 4 suites green.

## Shared factor-scan memo and large-ord cap — 2026-09-24

Two related finds in `x^n − 1` factorisation, the `≈ 0.6 s` cold cost
behind every divisor-base build at rungs like `n = 41` (`d = 20`):

1. **Shared scan.** `factor_x_n_minus_1` ran its own `2^ell` scan per
   process alongside `all_factors_of_x_n_minus_1`'s scan.  For
   `ell ≤ 24` it is now exactly the degree-`ell` slice of the cached
   full list (same discovery procedure, same scan order, `x + 1`
   excluded by degree) — one cold scan per process instead of two.
   Pure function, process-wide `OnceLock` memo, so hits are exact.
2. **Large-ord cap.** Degrees above 24 return empty instead of scanning:
   the old `ell < 40` bound admitted a `2^36` scan at `n = 37` that hung
   for hours (`factor_x_n_minus_1(37)` never returned).  Anything it
   could define has `≥ 2^25` abscissae — past every materialisation cap
   (`MAX_ABSCISSAE`) — so empty is the honest answer, and every caller
   already handles it (`ic compare` at n37 now fails fast with the
   dimension-limit error instead of hanging).

Gates: existing `x⁷ − 1` / `x⁹ − 1` value tests pin slice equality; new
`large_ord_factor_list_is_empty_not_a_hang` asserts `factor_x_n_minus_1(37)`
is empty in milliseconds and cross-checks the slice against the full
list at every `n ≤ 31`.  Measured: factor probe completes all rungs
(n37 instant-empty, n41 single 635 ms cold scan); n41 search unchanged
(inconclusive, 36 candidates); 119 + 14 + 4 suites green.

## Parallel LSD radix sort for pair-table build — 2026-09-24

Profiling showed the pair-table sort at ~40–60% of build time (8M
entries at n13: ~70–145 ms of ~175–280 ms).  `par_sort_unstable` is now
a parallel 16-bit LSD radix sort: per-chunk histograms, one serial
prefix, disjoint-offset scatter (raw-pointer writes justified by the
counting invariant, joined before reuse), adaptive pass count from a
key-OR pre-pass (n13's sub-2¹⁵ keys sort in a single pass), few wide
chunks (count/scatter are bandwidth-bound; the serial prefix costs
`chunks × 64K`).  Stable, so within-key order is the build's row-major
order; no gate depends on order beyond key-grouping (validity,
sortedness, verdict agreement — all green on both paths).
`KOBLITZ_SORT=comparison` restores the old sort for A/B.

Gates: new `radix_sort_orders_entries_by_key_stably` (edge cases,
adversarial duplicates, exact equality with stable `sort_by_key`);
52 lib + 14 + 4 suites green on both sort paths.  Measured under heavy
load, so wall-clock is reported as ranges with CPU-time confirmation:
full build min-of-8 187 ms (radix) vs 122 ms (comparison) wall —
overlapping — but process CPU over identical probe runs is 2.65–2.77 s
vs 3.64–3.81 s across 4 interleaved rounds (~1.37×), and end-to-end
`ic run` at n13 is 0.29 s vs 0.56 s CPU (~1.9×, build-dominated).

## Open problems from the talk (unimplemented)

- Couveignes–Lercier invariant factor bases via isogenies between
  algebraic tori and elliptic curves — a second, denser family.
  **Closed for this ladder (2026-10-08):** elliptic periods for
  `F_{2^n}/F_2` need an auxiliary elliptic curve over `F_2` with a rational
  point of order `n`; Hasse bounds `#E'(F_2) ≤ 5`, so the construction is
  empty for every `n > 5`, in particular at `n = 131`.  See
  `docs/ic/RESEARCH_ISOGENY_CONDUCTOR_GAP_20261007.md` §1.
- Exploiting the block/homogeneous structure of the resulting
  polynomial systems in the Gröbner step.  The systems are now actually
  built (`koblitz_groebner`), so this is measurable rather than
  hypothetical: the `m ≥ 3` chained systems are where it would pay.
- Precise complexity estimates across characteristics.
- Multiple invariant factor bases: a single one need not yield `n`
  *independent* relations, so the paper takes several.  The module
  exposes `factor_index` to select which factor of `x^n − 1` builds the
  base, but the driver uses one at a time.

## Tests

```bash
cargo test --release --lib cryptanalysis::koblitz_index_calculus
```

- `irreducibility_test_matches_known_polynomials` — Rabin's test.
- `find_irreducible_returns_degree_n_irreducible` — `n = 3 … 16`.
- `x_n_minus_1_factors_have_degree_ord_2_mod_n` — `x⁷ − 1`, `x⁹ − 1`.
- `point_counts_match_direct_enumeration` — the trace recurrence
  `#E = 2^n + 1 − s_n` against a full sweep of `F_{2^n}`, both `a`,
  `n = 3 … 9`.
- `frobenius_acts_as_lambda_on_the_subgroup` — `π = [λ]` across the
  subgroup, and `λ² − tλ + 2 ≡ 0`.
- `factor_base_is_frobenius_invariant` — `π(F) = F`, `−F = F`, orbit
  bookkeeping, orbits partition `F`.
- `subspace_is_closed_under_squaring`.
- `relation_rows_are_consistent_with_the_orbit_logs` — a relation row
  checked against brute-forced logs.
- `solves_the_dlp_on_k0_over_f512`, `solves_the_dlp_on_k1_over_f512`
  (algebraic oracle), `solves_the_dlp_by_enumeration_too`.
- `the_two_algebraic_engines_agree` — matrix-F4 vs Buchberger.
- `all_three_oracles_answer_every_target_identically` — search, Gröbner
  and SAT on 39 targets.
- `sat_refutes_undecomposable_targets`, `sat_handles_the_cubic_chained_system`,
  `solves_the_dlp_with_the_sat_oracle`.
- In `semaev_sat`: `xor_parity_encoding_is_exact_at_every_width`
  (regression), `cubic_monomials_are_encoded`, `contradictory_system_is_unsat`.
- `undecomposable_targets_are_rejected_algebraically` — 28 refutations
  with no branch searched.
- `speedup_model_tracks_the_predicted_factors`.

In `koblitz_groebner`:

- `symbolic_mul_matches_field_mul`, `symbolic_square_matches_mul_by_self`,
  `structure_constants_match_the_curve_field` — the Weil-restriction
  arithmetic against the field's own.
- `symbolic_s3_matches_scalar_s3` — the symbolic `S₃` against the scalar
  one from `binary_semaev`.
- `known_decomposition_is_a_root_of_the_system`, `system_is_quadratic`.
- `splitting_solver_agrees_with_exhaustive_evaluation` — the solver's
  roots against every point of `{0,1}^12`.
- `infeasible_system_is_closed_by_the_basis`.

## References

- S. D. Galbraith, R. Granger, S.-P. Merz, C. Petit, *On index calculus
  algorithms for subfield curves*, SAC 2020 / ePrint 2020/1315.
- J.-M. Couveignes, R. Lercier, *Elliptic periods for finite fields*,
  Finite Fields Appl. 15 (2009).
- I. Semaev, *Summation polynomials and the discrete logarithm problem
  on elliptic curves*, ePrint 2004/031.
- R. Gallant, R. Lambert, S. Vanstone, *Improving the parallelized
  Pollard lambda search on anomalous binary curves*, Math. Comp. 69
  (2000).
- N. Koblitz, *CM-curves with good cryptographic properties*, CRYPTO
  1991.
- J.-C. Faugère, L. Perret, C. Petit, G. Renault, *Improving the
  complexity of index calculus algorithms in elliptic curves over binary
  fields*, EUROCRYPT 2012 — the Weil-restriction-to-Boolean-system
  analysis the decomposition oracle implements.
- J.-C. Faugère, *A new efficient algorithm for computing Gröbner bases
  (F4)*, J. Pure Appl. Algebra 139 (1999).
- G. Bard, *Algebraic Cryptanalysis*, Springer 2009, ch. 13 — the
  Boolean-ring representation.
