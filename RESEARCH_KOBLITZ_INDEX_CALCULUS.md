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

## The resumable workflow driver — 2026-09-11

**Tool:** `ic workflow --params wf.json --dir runs/…` (`src/bin/ic/workflow.rs`).
**Docs:** `docs/ic/README.md`.

The stages above now run as a number-field-sieve-style workflow: a
parameter file names the curve, the factor-base source (an explicit
recipe or the census search), the oracle and the targets; the driver
executes **select → logs → solve** with every output on disk in a run
directory and a `state.json` manifest carrying the parameter digest and
each stage's status.  A rerun reloads existing artifacts, re-verifies
them against the reconstructed curve — a stale or tampered
`logs.json` is an error, never trusted — and continues from the first
incomplete stage; the solve stage rewrites `solutions.json` after each
target, so an interrupted batch resumes at the first unsolved one.  A
parameter file whose digest differs is refused, so one directory never
mixes two experiments.  `individual_log_with_pair_table` lets the batch
build the `|F|²` pair table once instead of once per target.

Exercised end to end on `K_0/2^9` (stop after select, resume through
logs, finish solve, full reuse, tamper rejection, digest refusal) and on
`K_1/2^15` with `mode: search`, where the census picked the same
30-point one-column pruned union as the interactive search and all five
targets descended in one relation each; the interactive `ic` tests cover
the stage-by-stage resume.

## Relation filtering and block Wiedemann — 2026-09-11

**Module:** `koblitz_sparse_la` (`filter_relations`,
`block_wiedemann_kernel`, `solve_sparse_system`).
**Wiring:** `KoblitzIcOptions::linear_algebra`, `ic logs
--linear-algebra sparse|dense --block-size`, the workflow's
`linear_algebra` block.  **Docs:** `docs/ic/README.md`.

The precompute solved its relation matrix by dense big-integer
elimination, re-run on the whole matrix after every new relation.  Each
relation has at most `m` nonzero entries, so the matrix is as sparse as
a sieve matrix, and the stages a number-field-sieve driver puts between
relation collection and the logarithm solve now exist for the ECDLP
pipeline:

- **Filtering** (`filter_relations`), in CADO's order: duplicate rows
  out; singleton columns eliminated with their row and recovered by
  back-substitution; excess rows removed down to a target excess by the
  clique rule (rows are joined by weight-2 columns with union–find, the
  heaviest row of the largest component goes first, and each removal
  cascades through the singletons it creates); light columns merged by
  structured Gaussian elimination under a fill-in bound.  Every
  elimination is kept as `(pivot row, column)`; reconstruction runs the
  stack backwards, then propagates through the *original* rows while
  some row has one unknown, then solves whatever is left densely — so
  the sparse path determines exactly the columns the original system
  determines (`filtering_preserves_the_solution_and_reconstructs_every_column`
  checks determinacy against an independent rank computation on random
  systems).
- **Block Wiedemann** (`block_wiedemann_kernel`).  The core is made
  square by folding its excess rows into random earlier rows with random
  coefficients, homogenised to `M (x, 1)ᵀ = 0` (never materialised: the
  operator applies `A` and subtracts `t·b`), and a kernel vector is
  found from the `m × n` Krylov blocks `X Mⁱ Y`.  The matrix
  Berlekamp–Massey step is the iterative shifted minimal approximant
  basis of `[a(λ) | I_m]` with shift `(0ⁿ, 1ᵐ)`: every basis column
  `[f; g]` has `a·f + g ≡ 0 (mod λ^L)` with `deg g ≤ δ − 1`, which is
  exactly the block recurrence `Σ_k S_{s−k} f_k = 0` for `δ ≤ s < L`.
  Two details that the tests forced:  (i) Horner runs from the top
  power (`f_0` multiplies `M^e`), and (ii) `Y = M·Z` for random `Z`
  (Coppersmith) — with random `Y` the vector identities
  `Σ_k M^{δ−k} Y f_k = 0` are exact in the block case and yield
  nothing, while evaluating the same generator on `Z` gives a vector
  `M` kills within `δ − e + 1` steps.  Sparse matrix-times-block
  products run in parallel over rows; `u64` arithmetic throughout
  (`r < 2^63`).
- **Safety.**  The solve is skipped until every column occurs in some
  row (a cheap test the dense path never made — it eliminated the full
  matrix to discover the deficiency), a failed Wiedemann run is retried
  with a fresh fold, the solution is checked against every relation,
  and the table is still certified in the group.  The two paths certify
  identical databases (`sparse_and_dense_linear_algebra_certify_the_same_log_table`
  on `K_0/2^9` and `K_1/2^11`; `ic`'s
  `factor_base_logarithm_database_precomputes_then_descends` compares
  the written files).

Measured:

| system | filter | core | block Wiedemann | solve |
|:-------|:-------|:-----|:----------------|------:|
| `K_0/2^31`, 35 columns, 74 relations (`ic logs`, sparse) | 4 singletons, 7 excess rows, 20 merged | 11 × 11, 158 nonzeros | 4×4, 14 Krylov terms, 17 products | 0.4 ms for all 40 attempts (dense: 22 ms) |
| synthetic, 3 000 columns, 18 000 rows of weight 3 (unit test) | 1 411 singletons, 14 968 excess rows, 1 273 merged | 316 × 316, 4 644 nonzeros | 4×4, 168 terms, 248 products | 48 ms |

At 35 columns the linear algebra was never the cost (the run is the
`|F|²` pair table, 16 s either way), so the point of the measurement is
that the sparse path is exact and no slower; the second row is the
regime it is for.  The dense reference is `O(rows · cols²)` big-integer
operations per attempt and keeps the full `rows × cols` `BigUint`
matrix in memory, which is what stopped the precompute from being run
on bases with thousands of columns; the sparse path keeps `≤ m`
entries per relation and costs one Krylov sequence.  The matrix
Berlekamp–Massey step is the quadratic iterative one (`O(L² m² (n+m))`),
fine to `10^4`–`10^5` columns; Thomé's subquadratic variant is the next
step if the base sizes get there.

## Open problems from the talk (unimplemented)

- Couveignes–Lercier invariant factor bases via isogenies between
  algebraic tori and elliptic curves — a second, denser family.
- ~~Exploiting the block/homogeneous structure of the resulting
  polynomial systems in the Gröbner step.~~  **Measured; it does not
  pay as a Macaulay-shift restriction.**  See
  [below](#the-block-structure-measured--2026-09-12).
- Precise complexity estimates across characteristics.
- Multiple invariant factor bases: a single one need not yield `n`
  *independent* relations, so the paper takes several.  The module
  exposes `factor_index` to select which factor of `x^n − 1` builds the
  base, but the driver uses one at a time.

## The block structure, measured — 2026-09-12

**Module:** `koblitz_groebner::{blocks, block_degrees, poly_block_degrees,
matrix_f4_f2_blocked}`.
**Bench:** `cargo run --release --example blocked_macaulay_bench`.

The open problem above supposed the decomposition systems have a block
structure worth exploiting, and that `m ≥ 3` is where it would pay.
The structure is real and stronger than expected; the exploitation does
not work, for a reason that is clear once the structure is written down.

**The structure.**  Partition the variables the way
`build_decomposition_system` lays them out — one block of `ℓ` per
summand, one block of `n` per intermediate point of the chain.  The
systems are **multilinear** with respect to that partition: degree at
most *one* in every block, at every `m`, even though the total degree
is 2 for `m = 2` and 3 once the chain appears.  Squaring is `F_2`-linear
in characteristic 2 and `x₁x₂` is bilinear, which is where it comes
from.  `the_decomposition_systems_are_multilinear_in_their_blocks` pins
it for `m = 2, 3`.

**Why it cannot be used as a shift restriction.**  A Macaulay matrix
multiplies each equation by shift monomials.  If `p` has block degree
`≤ 1` everywhere and `m` is a shift, then `p·m` has block degree
`1 + deg_i(m)` in block `i`.  So the *only* shift that preserves
multilinearity is the trivial one, and a Macaulay matrix bounded to
multidegree `(1,1,…,1)` is exactly the original equations — no new
information at all.  Allowing a bound of `2` per block admits shifts
again, but over `k` blocks that reaches total degree `2k`, which is a
*larger* matrix than the comparable total-degree one.

**Measured.**  Random `x_R`, so the systems do not decompose and the
event measured is the infeasibility certificate — rejection is most of
an attack's work, and refuting is what the note says the algebra is for.
Cost in 64-bit word XORs, summed over the targets where both engines
refuted:

| `n` | `m` | `v` | blocks | plain refuted | plain deg / ops | blocked refuted | blocked bound / ops | speedup |
|---:|---:|---:|:--|---:|---:|---:|---:|---:|
| 9 | 2 | 12 | `[6, 6]` | 5/8 | 2 / 169 | 5/8 | 1 / 169 | 1.00× |
| 9 | 3 | 27 | `[6, 6, 6, 9]` | 0/8 (capped) | — | 0/8 (capped) | — | — |
| 13 | 2 | 24 | `[12, 12]` | 3/8 | 2 / 660 | 3/8 | 1 / 660 | 1.00× |
| 13 | 3 | 49 | `[12, 12, 12, 13]` | 0/8 (capped) | — | 0/8 (capped) | — | — |
| 15 | 2 | 8 | `[4, 4]` | 8/8 | 3 / 19,383 | 8/8 | 2 / 91,925 | 0.21× |

Reading it:

- **`1.00×` is not a coincidence.**  Where the minimum setting refutes,
  the two matrices are *the same matrix*: the deciding total degree
  equals the system's own degree, so the only shift is trivial, and the
  blocked bound of 1 admits exactly that shift too.
- **Where they differ, blocked is worse** (`0.21×` at `n = 15`), for the
  reason above: a per-block bound of 2 over 2 blocks reaches total
  degree 4 where the plain matrix needed 3.
- **`m ≥ 3` is out of reach, not negative.**  At `v = 27` and `v = 49`
  the Macaulay matrix exceeds the size caps before any setting refutes,
  in both engines.  The note predicted `m ≥ 3` is "where it would pay";
  this says nothing either way about that, because neither engine gets
  there.  `(capped)` marks those rows so they are not read as results.

Class: **engineering**, and a negative one — the ratio to the floor does
not move and neither does `S`.

**What the real technique would be.**  Multilinearity is exploitable,
just not by restricting shifts.  Faugère, Safey El Din and Spaenlehauer's
work on bilinear and multihomogeneous systems uses the *Jacobian*
syzygies the bilinear structure forces — rows that are predicted to
reduce to zero and can be left out — and computes in multidegree with a
bound derived from the block sizes rather than a uniform one.  That is a
different algorithm, not a parameter change to this one, and it is
unbuilt.  The uniform per-block bound tried here is the cheap version of
the idea, and the cheap version is the one that does not work.

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

## Distributed relation collection — 2026-09-11

**Modules:** `koblitz_index_calculus::{RelationCollector, RelationWorkUnit,
CollectedRelation, probe_scalar, verify_collected_relation,
solve_factor_base_logs_from_relations}`.
**Tool:** `ic workflow --collect-units …` (worker), `ic workflow` (driver).
**Docs:** `docs/ic/README.md`.

Relation collection is the embarrassingly parallel half of the pipeline,
and a number-field-sieve run splits it into work units that clients
sieve independently and a server merges.  The workflow now has the same
shape.  The probe scalar of trial `t` is a function of the parameter
seed and `t` alone (a generator keyed by the pair), so the probe
sequence is one fixed sequence and a **work unit** `k` is the slice
`[k · unit_trials, (k + 1) · unit_trials)` of it.  `RelationCollector`
builds the decomposition state once per process (index map, field
structure, the `|F|²` pair table) and runs the trials of a unit in
parallel over the cores, returning the relations in trial order; the
single-process precompute `solve_factor_base_logs` draws its probes from
the same sequence in parallel batches of 64, so a local run and a
distributed one see identical relations (pinned by
`work_units_partition_the_probe_sequence_exactly`: the union of any
partition of a range, in any order, equals the whole range).

A worker (`ic workflow --params … --dir … --collect-units 4-7`) writes
`relations/unit-NNNNN.json` per unit and stops.  A relation file carries
only the probe scalar and the factor-base point indices of each
relation, bound to the parameter digest, curve, base recipe and summand
count.  The driver merges every valid unit present, **re-verifies each
relation in the group** (`[a]G == Σ P_i`, exactly `m` indices in range)
and drops exact duplicates before the linear algebra, so a corrupt or
forged file costs a rejected relation and nothing else, and a file from
another run or base is ignored rather than merged.  If the accepted
relations do not determine every column, the driver collects further
units itself up to `collection.max_units`, solving again after each,
and otherwise fails closed with the counts.  The workflow test runs a
worker for one unit, forges a relation in its file, drops in a file
from another digest, lets the driver fill the missing unit, and checks
that the forgery is the one rejected relation, the foreign file the one
ignored unit, and the three targets still descend; the library test
solves the same database from two out-of-order units with a duplicate
and a forgery mixed in and compares it column for column with the
single-process table.

Measured on `K_0/2^31` over the search-selected 35-column base
(`m = 3`, `unit_trials = 64`): two workers run concurrently on the same
four cores collected units 0 and 1 in 23 s each (each builds its own
`|F|²` pair table, which is the whole cost; the 64 probes are
milliseconds), then the driver merged 78 relations from the two units,
re-verified them, filtered the 78 × 35 system to a 13-column core, ran
block Wiedemann and certified all 35 logarithms in 0.7 s for the logs
stage, and descended three targets in two relations each — 13.8 s
wall for the driver, again dominated by its pair table.

This is the stage that makes the pipeline's collection cost scale with
machines rather than cores: every unit is independent, the merge is one
scalar multiplication per relation, and the linear algebra of the
previous entry is what turns the merged relations into the database.

## Beyond Koblitz: subfield curves `E/GF(2^k)` over `GF(2^{ke})` — 2026-09-11

**Modules:** `koblitz_index_calculus::{KoblitzCurve::subfield, invariant_factors,
top_factor_indices, subspace_basis_for_factors, q_linearised_kernel_basis,
subfield_group_order, frobenius_eigenvalue_q}`; `binary_semaev::solve_artin_schreier`.
**Tool:** `ic run|search|logs|solve --subfield k --curve-a a --curve-b b`,
workflow `curve.subfield` / `curve.curve_b`.
**Docs:** `docs/ic/README.md`.

Everything the pipeline does with a Koblitz curve uses one fact: the
curve is defined over a subfield, so a power of Frobenius is an
endomorphism acting as a known scalar `λ` on the prime-order subgroup,
and the invariant factor bases are the `π`-stable subspaces.  The
pipeline now takes that fact at its generality (GGMP §4 with `q = 2^k`):

- `KoblitzCurve::subfield(k, n, a, b)` builds `y² + xy = x³ + ax² + b`
  with `a, b ∈ GF(2^k) ⊂ GF(2^n)`, `n = k·e`, `e` odd, the coefficients
  named by coordinates in the `F_2`-basis of the subfield (the kernel of
  `X^{2^k} + X`).  `#E(GF(2^k))` is counted by enumeration, the trace
  `t` of the `q`-power Frobenius follows, and `#E(GF(2^n))` comes from
  `s_i = t·s_{i−1} − q·s_{i−2}` — the Koblitz recurrence with `2`
  replaced by `q`; `λ` is the root of `λ² − tλ + q` with `π(G) = [λ]G`.
  `KoblitzCurve::new(a, n)` is `subfield(1, n, a, 1)` and a test pins
  every field of the two equal, so the Koblitz path is unchanged.
- The invariant subspaces are the kernels of `q`-linearised polynomials
  `Σ c_i X^{q^i}`, `c_i ∈ GF(q)`, classified by the irreducible factors
  of `x^e − 1` over `GF(q)`.  Those are found by Cantor–Zassenhaus over
  `GF(q)` run inside `GF(2^n)` (`F2mPoly` arithmetic with coefficients
  in the subfield: distinct-degree splitting by `gcd(x^{q^d} − x, ·)`,
  equal-degree splitting with the absolute-trace map of random
  `GF(q)[x]` elements), sorted canonically; for `k = 1` the list is the
  `F_2` bit-mask list in the same order, so recipe indices, documents
  and reports are unchanged.  A factor of degree `d` gives `2^{kd}`
  abscissae and orbits of length dividing `e`.  Every factor-base
  constructor, the search families, the pair table, the cofactor
  admissibility walk, relation rewriting, filtering, block Wiedemann,
  the log database and the descent then work as they are, with
  `kc.frobenius` the `q`-power map; the only oracle-side change is the
  Kosters–Yeo trace row's constant, `Tr(a)` for `a` in the subfield.
- The Artin–Schreier solver used to brute-force `2^m` candidates for
  even `m` (and refuse `m > 20`).  Even field degrees are the normal
  case for even `k`, so `t² + t = c` is now solved as the `F_2`-linear
  system it is (bit-packed elimination, `m ≤ 63`), which took a
  `GF(2^18)` point lookup from 29 ms to microseconds.

Two structural observations from the first instances:

- Point counts check out against full enumeration on `E_{0,2}/GF(4)`
  over `GF(2^6)`, `GF(2^10)` and `GF(2^14)` (76, 964 and 16 636 points);
  `E_{0,2}/GF(4)` over `GF(2^14)` has `r = 4159`, `h = 4`, and `x^7 − 1`
  factors over `GF(4)` as degrees `(1, 3, 3)`, giving invariant bases
  of 3, 43 and 71 points — the 71-point one solves the DLP and
  precomputes a certified log database with two summands (pinned by
  `subfield_factor_bases_are_invariant_and_solve_the_dlp`).
- When `x^e − 1` splits into binomials `x^d − c` over `GF(q)` (`e = 9`,
  `q = 4`: degrees `1, 1, 1, 3, 3`), each invariant subspace minus `0`
  is a multiplicative coset `{x : x^d = c}` and its inverses form the
  reciprocal factor's coset, which has absolute trace `0`.  With
  `Tr(a) = 1` and `b = 1` the Artin–Schreier condition
  `Tr(x + a + b/x²) = 0` then fails for every `x ≠ 0` in every such
  subspace: `E_{2,1}/GF(4)` over `GF(2^18)` has *no* factor-base points
  at all beyond `(0, √b)`, on any of its five invariant subspaces, while
  `b ∈ {ω, ω²}` on the same field restores them.  The search scores
  such bases at zero and a run over one reports a base with no usable
  columns; the choice of `b` is a genuine parameter of the family, not
  a cosmetic one.

## Scaling to the ledger rungs, and a ρ baseline that was lying — 2026-09-11

The pipeline now runs end to end at the boundary ledger's Koblitz rungs
— `K_0/F_{2^n}` for `n = 31, 37, 39, 41`, 32 synthetic known-answer
targets each, with the in-process signed-Frobenius ρ baseline. The
parameter files are `docs/ic/params/k0n{31,37,39,41}.json` and the
measurements `docs/ic/runs/koblitz-scaling-20260911.json` (two
consecutive series from one build on one host).

### Single-word arithmetic

Every hot loop used to run on `Vec<u64>`-backed field elements that
allocate on each operation, which for the `n ≤ 62` fields this pipeline
actually uses costs about two orders of magnitude more than the
arithmetic. `cryptanalysis::koblitz_fast` puts a point in two `u64`s on
top of the existing carry-less `Gf2`:

- `add`, `double`, `mul`, the `2^k` Frobenius, and the packed identity,
  all checked against `binary_ecc::curve` on random points including
  `O`, `P + P`, `P + (−P)` and points outside the prime-order subgroup;
- `add_many` (one summand, many addends) and `add_pairwise` (matched
  slices), each spending **one** field inversion for the whole slice by
  Montgomery's trick, which is what a point addition's inversion costs
  when enough of them happen together.

The pair-sum table builds its rows with `add_many`, and answers a lookup
through a bucket index plus a presence filter — one bit per hashed key,
four bits per entry — so the overwhelming majority of lookups, which
miss, never touch the hundreds of megabytes of entries. Measured at
`n = 41`: the descent went 242 → 157 ms per target from the filter
alone. The descent's target-independent setup (signed-orbit map, column
logarithms, oracle tables) moved into `IndividualLogSolver`, built once
per batch rather than once per target.

Degree 31, 32 targets, before → after the single-word rewrite: descent
715 → 1.6 ms per target, collection 12.4 → 1.1 s.

### The ρ baseline was failing, not losing

The first version of the scaling run reported a **charged crossover of
1.009 at `n = 41`**. It was an artifact, and the way it failed is worth
recording.

That baseline used Floyd cycle finding. On a negation-quotient walk the
canonical representative flips the sign of the coefficients, so a walk
can return to *its own previous state*: going round such a cycle adds a
jump and then subtracts it, and the accumulated coefficients come back
unchanged. Floyd's meeting then happens at that cycle — tortoise and
hare holding identical `(a, b)` — and the collision carries no
information. Every restart burned itself on one within a few hundred
steps: at `n = 41` the baseline recovered **0 of 32** logarithms while
charging 1.4M iterations, and a ρ that never finishes is trivially
"slower" than anything.

The replacement is the distinguished-point method:

- one trajectory per walk; a direct-mapped cache of recent points
  (4096 slots) catches short cycles *and* nearby collisions; a sparse
  table of stored points — about one in `2^6` of the expected walk
  length — catches the long-range collision;
- a cycle is escaped by doubling **the cycle's own smallest state**, so
  every path entering the same cycle leaves it at the same point and the
  walk stays a deterministic map; a repeat doubles further;
- walks are stepped in batches sharing one field inversion
  (`parallel_walks`, default 32), scaled down when the whole walk is
  shorter than the setup would cost — at `n = 31` a 171-step walk would
  otherwise spend more on 32 jump tables than on walking.

The general-arithmetic walk stays as `koblitz_signed_frobenius_rho_reference`
and is asserted equal step for step, on every charge counter. Two tests
hold the walk to theory: the mean step count against
`√(πr/2) / √(2n)`, and the degree-41 rung (ignored by default).

### What the rungs actually measure

| `n` | `r` bits | base points | columns | relations | core | collect | logs | trials/target | descent/target | ρ/target | ρ steps | charged ρ/IC |
|-----|---------|-------------|---------|-----------|------|---------|------|---------------|----------------|----------|---------|--------------|
| 31 | 21 | 2170 | 35 | 78 | 13 | 1.1 s | 0.6 s | 2 | 1.6 ms | 1.2 ms | 171 | 0.73 |
| 37 | 28 | 3851 | 26 | 39 | 6 | 1.5 s | 1.4 s | 62 | 22.9 ms | 4.8 ms | 2259 | 0.21 |
| 39 | 27 | 5151 | 33 | 49 | 7 | 12.2 s | 8.2 s | 125 | 52.8 ms | 2.7 ms | 1072 | 0.05 |
| 41 | 40 | 4759 | 29 | 52 | 10 | 2.0 s | 5.3 s | 387 | 141 ms | 90 ms | 126425 | 0.64 |

All 32 targets solve and all 32 ρ runs verify at every rung. **No rung
crosses**: the charged ratio is below 1 everywhere, and the closest
(`n = 41`, ≈ 0.64–0.73) is the rung where ρ has the most room left,
not the least.

A third series was run later on a verifiably idle host, from the
committed parameter files, to check the first two for measurement load.
Per-target times came out 3% to 20% faster at `n = 37` and `n = 39`, so
the first two series did carry some background load — but the ρ step
counts are identical and the charged ratios agree across all three
(0.72/0.74/0.75, 0.20/0.21/0.21, 0.067/0.051/0.060, 0.73/0.64/0.66),
because load moves both sides of a same-process comparison together.
That is the reason to quote the ratio rather than the wall time: the
ratio is what survives an imperfectly quiet machine.

Three things the table says plainly:

1. **The descent is decomposition-bound.** Trials per target track
   `1/coverage` — 2 at `n = 31`, 387 at `n = 41` — and each trial is a
   pair-table search over the whole base. Everything else (linear
   algebra at 0.2 ms, relation collection, the log database) is noise
   by comparison. A faster descent means a base with better coverage
   per abscissa, not faster field arithmetic.
2. **The linear algebra has stopped being the problem.** Filtering
   reduces 26–35 columns to cores of 6–13 and block Wiedemann solves
   them in under 0.3 ms; `n = 39`'s 5151-point base costs 12 s to
   *collect* relations for and 8 ms to solve.
3. **`n = 39` is the wrong kind of rung.** Its cofactor is 8012 and its
   subgroup only 27 bits, so ρ finishes in a thousand steps while the
   base is the largest of the four. Degree 41, with `h = 4` and a 40-bit
   subgroup, is the honest one.

### Where the remaining headroom is

- **Coverage per abscissa.** At `n = 41` the census measures coverage
  `2^-8` on a 4759-point base: 387 trials per target. The union and
  2-torsion-saturated families at seed dimension 6 are the largest
  bases inside the 4096-abscissa cap; past that the pair table is
  quadratic in abscissae, so the next gain has to come from bases whose
  *witness* count per point is higher, not from more points.
- **ρ still has room**, which is the honest reading of any ratio near 1:
  the walk is a single core with batched inversions, while the published
  records for this size use many cores and tighter inner loops. Until a
  rung beats a ρ that has had the same attention, a ratio near 1 is a
  statement about two implementations, not about two algorithms.

## The factor base was the wrong shape — 2026-09-12

Every factor-base family in this pipeline was linear: an invariant
subspace, a union of Frobenius translates of one, or a subset of those.
That is not a design choice, it is an inheritance. Semaev's polynomials
and the Weil descent are *written over a subspace*, so an algebraic
oracle needs one. The pair-table oracle does not — it is a
meet-in-the-middle search that never looks at the base's structure — and
the structure turns out to be expensive.

### What the structure costs

At `n = 41`, over 3000 random subgroup targets with `m = 3`, measuring
how often a target is a sum of three base points:

| base | points | columns | decomposes | vs `|F|³/(3!·r)` |
|------|--------|---------|------------|------------------|
| 2-torsion-saturated union, seed dim 6 | 4759 | 60 | 1 in 273 | 8.9× worse |
| Frobenius-closed, random abscissae | 5494 | 67 | 1 in 67 | 3.4× worse |
| **drawn from the prime-order subgroup** | 4838 | 59 | **1 in 28** | **matches** |

Two separate effects, and the second is the larger:

1. **Linearity.** The sums of a subspace union concentrate: the sumset of
   a structured set is far from uniform, so its sums hit a given target
   less often than three random points would.
2. **The cofactor.** A relation and a descent both ask for a target *of
   `⟨G⟩`* to be a sum of base points. A point outside `⟨G⟩` can only sum
   into it when the cofactor parts of the summands cancel, so a base
   drawn from the whole curve spends most of its sums in the wrong
   coset. Restricting the base to `⟨G⟩` removes that loss entirely — the
   measured rate lands exactly on the random-base expectation.

The column count is unchanged (`≈ |F|/2n` either way), so the linear
algebra does not pay for the improvement.

### Selecting a subgroup base without knowing any logarithm

`build_subgroup_orbit_factor_base(kc, seed, points)` draws a random
abscissa, lifts it to a point `P`, and takes `[h]P` for the cofactor
`h`. That lands in `⟨G⟩` for **every** `P`, so no draw is wasted, and
knowing `P` says nothing about `log [h]P`.

The first version tested `[r]P = O` and rejected what failed. That is
also logarithm-free, but at `n = 39` the cofactor is 8012, so it threw
away 8011 of every 8012 draws and the select stage took **607 seconds**.
Multiplying by the cofactor instead made it immeasurable. When a
rejection test and a construction agree, take the construction.

### What it does to the rungs

Same curves, same 32 targets, same ρ baseline, same process; only the
factor base changes:

| `n` | descent/target, union → subgroup | ρ/target | charged ρ/IC | crossovers (charged / amortised / whole-process) |
|-----|----------------------------------|----------|--------------|--------------------------------------------------|
| 31 | 1.5 ms → 1.3 ms | 1.2 ms | 0.95 | no / no / no |
| 37 | 18.3 ms → 2.9 ms | 5.8 ms | 2.01 | **yes** / no / no |
| 39 | 44.1 ms → 2.7 ms | 3.6 ms | 1.36 | **yes** / no / no |
| 41 | 152.9 ms → 14.4 ms | 119.5 ms | **8.33** | **yes** / **yes** / **yes** |

All 32 targets solved and all 32 ρ runs recovered and verified at every
rung. Evidence: `docs/ic/runs/koblitz-subgroup-bases-20260912.json`.

**Read the degree-41 row carefully.** The whole-process verdict compares
one precompute plus 32 descents (3.4 s) against 32 ρ walks (3.8 s). It is
a *bulk* statement. For a **single** target the precompute is paid whole,
so IC costs 3.0 s against ρ's 0.12 s and ρ wins by about 24×. The
crossover is in amortising a database over many targets, which is what
the CADO-style split was built for and what the ledger's charged class
was defined to measure.

Two more things this is not. It is not asymptotic: the work per target is
still `Θ(r/|F|²)` with `|F|` capped by the pair table's memory, so this
moves the constant and not the exponent. And ρ has room left — the
baseline here is one core with batched inversions at about 1.2 µs a step,
where a tuned implementation would do better — so the honest reading of
a ratio near 1 is still "two implementations", not "two algorithms". At
8.3× the degree-41 charged ratio has more margin than that objection can
absorb, but the amortised and whole-process verdicts there (1.11 and a
0.4 s gap) do not.

### A driver that pays its setup once

Sizing collection to the base exposed a waste: when the first units do
not determine the columns, the driver collected another unit and called
the solver again — and the solver rebuilt the signed-orbit map of the
base and re-verified *every* relation each time, quadratically.
`FactorBaseLogSolver` now holds that setup and takes relations
incrementally, verifying each exactly once. At `n = 39` the logs stage
went from 16.6 s to 3.6 s over four extension rounds, at `n = 31` from
8.0 s to 1.1 s over nine.

### How large should the base be? — 2026-09-12

With a subgroup base the decomposition rate is the one a random base
gives, so the cost model is finally clean and can be used to *choose* a
base rather than to explain one. A target decomposes with probability
`|F|³/(3!·r)`, and an `m = 3` trial costs `|F|` table lookups, so the
work per target is

```text
    (3!·r / |F|³) · Θ(|F|)  =  Θ(r / |F|²)
```

— four times cheaper for every doubling of the base, until the rate
saturates at one witness per trial and growth only makes each trial
dearer. The pair table is quadratic in `|F|`, so the precompute rises
four times per doubling at the same moment.

Measured at `n = 41` (`docs/ic/runs/koblitz-base-size-20260912.json`):

| points | abscissae | table | build | trials/target | ms/target |
|--------|-----------|-------|-------|---------------|-----------|
| 1312 | 656 | 14 MB | 0.1 s | 400 | 47.8 |
| 2624 | 1312 | 55 MB | 0.4 s | 133 | 29.2 |
| 5248 | 2624 | 220 MB | 1.7 s | 29 | 13.7 |
| 10496 | 5248 | 881 MB | 6.6 s | 3 | 3.5 |
| 20992 | 10496 | 3.5 GB | 32.3 s | 1 | 2.3 |

The `1/|F|²` law holds exactly to 10496 points, where the rate saturates
and the last doubling buys 1.5× for four times the memory.

**So the base is a function of the target count, not of the curve.**
Total cost is `precompute(|F|) + T·descent(|F|)`, and end to end over
256 targets at `n = 41`:

| base | precompute | descent/target | total | ρ total |
|------|-----------|----------------|-------|---------|
| 5248 points | 3.0 s | 16.2 ms | **7.2 s** | 31.5 s |
| 10496 points | 9.0 s | 7.2 ms | 10.9 s | 31.1 s |

256 of 256 solved, 256 of 256 ρ walks recovered and verified, both runs.
The two lines cross at about **665 targets**: below that the smaller base
wins on total cost, above it the larger. Against ρ's 121 ms per target
they cross at roughly 29 and 85 targets respectively — which is the
honest way to state a whole-process win, as a *number of targets* rather
than a verdict.

The abscissa cap was `2^12`, which put the measured optimum out of
reach; it is now `2^13`. What the cap was really protecting is the pair
table, so `PairSumTable::build` now refuses past a stated byte budget
(4 GiB by default, about 16000 points) instead of attempting an
allocation the machine cannot meet — a base one doubling too large now
fails with a number rather than an OOM.

### Collection and the descent need not agree on `m`

They share the factor base and its pair table. They do not share the
number of summands, and they should not: the two have different cost
shapes.

```text
    three summands:  3!·r/|F|³ probes  ×  |F| lookups each
    two summands:    2r/|F|²   probes  ×  1 lookup each
```

Collection wants **few** probes, because each one costs two scalar
multiplications to build — so it takes three summands and does a long
scan per probe. The descent can make a probe nearly free: a probe is any
`[a]G + [b]Q`, so step one by `+G` instead of drawing a new pair, and
step 64 of them together so a single field inversion serves the lot.
Once a probe costs less than the lookup after it, two summands win.

The walks have to start a stride apart on the same line rather than from
64 independent draws — 3 scalar multiplications and 63 additions instead
of 128 scalar multiplications. On a walk of only a few hundred rounds
that setup was most of the descent: it cost 9.4 ms a target before the
change and 4.65 ms after.

Measured at `n = 41`, collection at three summands throughout:

| base | descent `m` | ms/target | probes/target | charged ρ/IC |
|------|-------------|-----------|---------------|--------------|
| 5248 | 3 | 14.35 | 21 | 8.3 |
| 5248 | **2** | **8.63** | 43529 | **13.6** |
| 10496 | 3 | 7.20 | 2 | 16.5 |
| 10496 | **2** | **4.65** | 9427 | **24.8** |

All 32 targets solved in every configuration. `descent_summands` in a
workflow parameter file selects it; it defaults to `summands`.

A floor is now visible underneath: a solved target ends with one
`[d]G = Q` check in the general arithmetic, which costs 1.16 ms at this
degree against 0.025 ms in the single-word representation. It stays in
the general arithmetic on purpose — a check that shares no code with the
arithmetic that did the search is worth more than the milliseconds — and
ρ pays the same on its own candidate. At 8.6 ms a target it is 13% of
the descent, and at the wide base a quarter of it.

## Degree 53, and where this stops winning — 2026-09-12

`K_0/F_2^53` carries the largest subgroup in this family: 44 bits,
`r = 21 044 858 204 113`, 38 times the degree-41 rung. It runs end to
end, with a 15264-point subgroup base, collection at three summands and
the descent walking two:

| stage | |
|-------|--|
| base | 15264 points, 7632 abscissae, 144 columns |
| relations | 458 from 15888 probes |
| precompute | 35.9 s (collect 17.7, logs 9.2, pair table the rest) |
| descent | **49.2 ms** a target, 243863 probes |
| ρ | **1.216 s** a target, 996618 steps |
| verdict | charged **24.7**, amortised 1.04, whole-process **yes** |

32 of 32 targets solved, 32 of 32 ρ walks recovered and verified.

The first attempt lost one target of the 32 — not to the mathematics but
to `max_trials`, which was capped at a million. That ceiling was set when
a probe cost two scalar multiplications and a million of them was an
hour; a walked probe costs 0.2 µs, so a million is a fifth of a second
and the cap was cutting off targets that happened to need five times the
mean. It is now a hundred million.

### Where it stops

At a fixed memory budget the factor base is fixed, and then everything is
determined:

```text
    descent        2r/|F|²      probes        — linear in r
    collection     9r/(n|F|)    lookups       — linear in r
    ρ              √(πr/2)/√(2n) steps        — square root of r
```

Linear loses to a square root. The only question is where, and the
constants this run measured — 0.20 µs a probe, 1.22 µs a ρ step, a
15264-point base — answer it:

| subgroup bits | descent | ρ | charged | precompute | break-even targets |
|---------------|---------|---|---------|------------|--------------------|
| 40 | 0.002 s | 0.156 s | 82 | 2 s | 16 |
| 45 | 0.061 s | 0.881 s | 14 | 79 s | 96 |
| 48 | 0.488 s | 2.49 s | 5.1 | 632 s | 315 |
| 50 | 1.95 s | 4.98 s | 2.6 | 2528 s | 834 |
| 52 | 7.80 s | 9.96 s | 1.3 | 10112 s | 4679 |
| 53 | 15.6 s | 14.1 s | 0.9 | 20224 s | never |
| 56 | 125 s | 39.9 s | 0.3 | — | never |

**The charged advantage survives to about a 52-bit subgroup and is gone
by 53.** The whole-process advantage dies sooner in practice: it needs a
batch of targets that grows with `r` — 16 at 40 bits, 96 at 45, 315 at
48, 834 at 50 — because the precompute grows linearly while ρ grows as a
square root.

So **relation collection is what stops this, not the descent**. The
descent is the part that has been optimised hardest and it is not the
binding cost past 48 bits; collecting the relations that build the
database is.

More memory moves the boundary by about half a bit per doubling of `|F|`
— `|F|` enters the descent squared, so a 4× table buys one bit. That is
the shape of the whole result: every improvement in this note moved a
constant, and the exponent has not moved at all.

A last caveat on the table: it is arithmetic on two measured constants,
not a series of runs. Checked against the degree-53 run it under-predicts
the win by about 1.7×, mostly because that run's ρ walks took 1.8× their
expected step count.

## Collection, taken to its floor

The previous section ended by naming relation collection as the binding
cost past 48 bits. It was, and it was spending three times what it had
to.

A three-summand probe asks whether a random `R = [a]G` is a sum of three
factor-base points, and it asks by scanning: for each `k` it looks up
`R − P_k` in the pair-sum table. If `R = P_i + P_j + P_k` the scan finds
that triple three times — once with each of its indices standing as the
third summand — and the code kept exactly one of the three, the one whose
witness came out sorted (`j ≤ k`). The other two chances were computed
and thrown away.

They needn't be. Scanning a **window** of `w` summands instead of all of
them, and dropping the sorting condition, keeps all three chances at
`w/|F|` of the cost: a triple is caught whenever *any* of its three
indices lands in the window, which happens with probability
`1 − (1 − w/|F|)³ ≈ 3w/|F|`. Three times the relations per lookup, paid
for with more targets rather than more scanning — so the targets have to
be cheap, which is why the probes are now walked (runs of 64 trials
sharing one scalar multiplication and stepping by a seed-fixed stride)
rather than multiplied one at a time.

There is a floor to this, and the window reaches it. A single lookup
succeeds with probability `|F|²/2r`, so no oracle of this shape can spend
fewer than

```text
    2r/|F|²   lookups per relation
```

Measured against that bound, the full scan and a window of `|F|/32`:

| | lookups per relation | ideal | ratio |
|---|---|---|---|
| degree 41, full scan | 119 273 | 39 922 | 2.99 |
| degree 41, window `|F|/32` | 41 323 | 39 922 | **1.04** |
| degree 53, full scan | 529 507 | 180 651 | 2.93 |
| degree 53, window `|F|/32` | 192 973 | 180 651 | **1.07** |

The waste was exactly the factor of three the argument predicts, and
what is left is within 7% of the bound. There is nothing further to win
here without a different oracle.

Note what the limit says: `2r/|F|²` lookups per relation is the same
count as the walking descent's `2r/|F|²` probes. As the window narrows
to a single summand, a collection probe *is* a two-summand descent probe
— one target, one lookup. Collection and the descent were never two
algorithms; the full scan was the only thing making them look different.

### What it costs in practice

At `K_0/F_2^53`, 32 targets, against a full-scan control run in the same
session on the same binary:

| | full scan | window `|F|/32` |
|---|---|---|
| trials | 15 888 | 174 768 |
| relations | 458 | 432 |
| summands scanned | 242.5 M | 83.4 M |
| collect stage | 14.49 s | 9.78 s |
| precompute | 32.08 s | 27.29 s |
| amortised ratio vs ρ | 1.155 | **1.331** |

Both solved 32 of 32 and had all 32 confirmed by the ρ walk. The collect
stage carries the one-time pair-table build — about 6.5 s at this base
size, the same for both — and net of it the scanning costs 17.5 ms a
relation against 7.6 ms, about 2.3 times. That is short of the 2.7 times
fewer lookups because a windowed lookup costs some 20% more in wall
time: its prefetch pipeline restarts every window instead of running the
length of the base. Batching several windows into one scan to fix that
was tried and measured; it was a wash, and the simpler code stands.

### A cost that was not algorithmic at all

The same measurement turned up something else. The collector's
point-index map is keyed by big integers, and it was being rebuilt for
every work unit. It costs about as much to build as a short unit costs to
run, so the price of a run depended on how finely its trials had been
partitioned — precisely the parameter a distributed collection is
supposed to be free to choose:

| | 1 unit | 4 units | 8 units |
|---|---|---|---|
| degree 41 collect, before | 0.92 s | 1.25 s | 1.70 s |
| degree 41 collect, after | 0.92 s | 0.97 s | 0.98 s |

Building it once per stage removed it. At degree 53 the same fix took the
full-scan collect stage from 25.20 s to 14.49 s — larger than the window
saves, on a base of 15264 points. Together the two changes take the
degree-53 collect stage from 25.20 s to 9.78 s.

It is worth saying plainly that the second of these was not a discovery
about index calculus. It was a cost that had been sitting in the
distributed collector since it was written, invisible because nothing had
ever compared a run against itself at two partition sizes.

## What was left once collecting stopped being the cost

Flooring collection had a side effect worth reporting on its own. With
the scan no longer dominating, the degree-53 precompute turned out to be
three fixed setup costs of about eight seconds each:

```text
    select the factor base          8.4 s
    projected signed-orbit map      8.2 s   (built twice per run)
    pair-sum table                  7.5 s
```

Two of the three were accidental.

**The orbit map.** It projects every base point by the cofactor and then
walks each one's whole signed Frobenius orbit — about 1.6 million point
operations on a 15264-point base — and it did all of it in big-integer
arithmetic while the rest of the pipeline had long since moved to single
words. Moving it costs an argument rather than a leap of faith: the
general map canonicalises and sorts by `point_key = (x + 1, y)`
lexicographically, and `FastPoint::pack` is `((x + 1) << 1 | sign)` where
`sign` marks which of `y` and `x ^ y` is the larger — so within one
abscissa it separates `P` from `−P` in exactly the order `y` does. The
two orders agree, the two maps choose the same representatives, and a
test compares them point for point at four degrees. **8.2 s → 0.047 s**,
twice per run.

**Selecting the base.** It drew abscissae in batches of eight and rebuilt
the entire base after each batch — 954 rebuilds over a growing list to
reach 15264 points, quadratic in the abscissae. It need not rebuild at
every boundary: an abscissa carries at most two points, so `2·|abscissae|`
bounds what a representative set can yield, and while that bound is below
the target a build could only have come back short and the loop would
have gone round again. Skipping those builds changes nothing — same
15264 points, same 144 columns. **8.4 s → 0.93 s**.

The pair table's 7.5 s is real work: `|F|²/2` point additions and a
1.9 GB sort. It is now the dominant fixed cost, and unlike the scanning
it does not grow with `r`.

### The precompute over this line of work

Measured at `K_0/F_{2^53}` on a 15264-point base, each cost in isolation:

| | before | after |
|---|---|---|
| select the factor base | 8.60 s | **0.94 s** |
| projected signed-orbit map | 8.16 s | **0.048 s** |
| pair-sum table | 6.67 s | 6.67 s (untouched) |

End to end on the parameters the repository ships, 32 targets:

| | precompute | amortised ratio vs ρ |
|---|---|---|
| before any of this | 35.9 s | 1.04 |
| with the orbit map alone | ~24.1 s | — |
| with both | **16.5 s** | **2.16** |

A collection window of `|F|/32` takes it further still, to 14.8 s and a
ratio of 2.41, but a window is chosen on its own selection run rather
than carried by a rung parameter file, so the table above is the
like-for-like reading. The window and these fixes are independent and
their effects compose.

32 of 32 solved and 32 of 32 confirmed by ρ at every step. The
whole-process advantage at this degree went from marginal to a factor of
two and a half, and none of it came from a better algorithm — the cost
laws in the previous sections are unchanged. What changed is that the
implementation now spends its time where those laws say it should.

It is worth being blunt about the pattern. Both of these had been sitting
in the pipeline since the code was written, and neither was visible until
something else stopped being the bottleneck. The window made collection
cheap enough that eight seconds of setup became the thing to look at;
before that, they were noise inside a larger number. Profiling only finds
what is currently largest, which means a real speedup is usually followed
by another one.

## The memory wall, and what a pair costs to store

Every section above moved a constant inside a cost law. This one moves
the quantity the laws divide by.

The descent needs `2r/|F|²` probes and collection `2r/(n|F|)` lookups, so
`|F|` is what the reach is made of — and `|F|` is capped by the pair
table, which holds `|F|²/2` sums. At sixteen bytes a sum, four gibibytes
buys a base of about 23000 points and not one more. The relation is
worth writing down: at a fixed budget `B`,

```text
    |F| = √(2B / bytes per pair)
```

so the descent's `2r/|F|²` is directly proportional to the width of a
stored pair. **Halving the bytes halves the descent.**

Sixteen bytes is `(packed sum, i, j)`. The sum is what a lookup asks
about; `i` and `j` are what it answers with. Four bytes hold a bucketed
hash of the sum — never a false negative, and a false positive about one
time in `2²⁸`. What that gives up is the answer: a hit knows a
decomposition may exist but not of what.

A hash rather than the key's own low bits, which would be exact if the
bucket index were wide enough to leave 32 of them over. "Wide enough" is
`key_bits − 32`, and a packed sum is `n + 2` bits: at `n = 53` that asks
for 23 bucket bits while run length already wants 26, so it costs
nothing, but at `n = 61` it asks for 31, and an index of `2³¹` words is
four gibibytes spent on nothing but making a rest fit. A 36000-point base
needs 2.97 GiB hashed at any degree, against 6.72 at `n = 61`. Neither
form can give a wrong answer — build and lookup narrow a key the same
way — but only one of them fits.

The answer is recoverable, and recovering it is what makes the
approximation safe. `target − P_i` is a base point exactly when `i` is a
summand, so one `|F|`-long scan finds the pair — and the scan asks the
group, not the table, so a false positive costs a scan that comes back
empty and never a wrong answer. Hits are rare by construction — that is
what a decomposition oracle *is* — so the scan is paid about once per
relation rather than once per probe. Storage falls from sixteen bytes to
about four and a half, counting the bucket index and the presence filter.

At four gibibytes that is a base of 42302 points instead of 23169, and
3.33 times fewer descent probes.

### What it costs, measured

`K_0/F_{2^53}`, 32 targets, the same four-gibibyte budget:

| | narrow (15264 points) | wide (36464 points) |
|---|---|---|
| pair table | full, 16 B/pair | compact, 4.9 B/pair |
| precompute | 16.5 s | 87.8 s |
| descent probes/target | 243 863 | **32 877** |
| descent | 50.4 ms | **16.3 ms** |
| charged ρ/IC | 25.0 | **75.2** |
| amortised ρ/IC | 2.17 | 0.44 |

Both solved 32 of 32 and had all 32 confirmed by ρ. The charged ratio
triples, which is about three and a half bits of reach since the ratio
falls as `√r`. The precompute quintuples, because the table is quadratic
in a base that grew and the base has 344 columns to determine instead of
144.

So this is a trade, not a free win, and the crossing point is the target
count: 71 extra seconds of precompute against 32 ms saved per target is
**about 2200 targets** before the wider base is the cheaper one. For a
single logarithm the narrow base is better; for a campaign against a
whole curve the wide one is, and it is the only one of the two whose
charged advantage survives past 53 bits.

### Building it

The compact table is built by a counting sort in two parallel passes and
no comparison sort at all — the bucket of a key is known before any key
is compared, so counting the buckets and scattering into them puts every
rest in its run, and a run is sixteen entries, one cache line, scanned
rather than searched. The pair sums are recomputed in the second pass
rather than held between the two, because holding them would cost the
eight bytes a pair the representation exists to avoid.

It costs 66.3 s for 664.8 million pairs against 6.7 s for 116.5 million —
100 ns a pair against 57. That factor is the second pass, and it is now
the largest single item in a wide precompute.

## Degree 61: a 48-bit subgroup, and the boundary moves

`K_0/F_{2^61}` has a 48-bit prime-order subgroup,
`r = 162 888 033 982 417` — 7.7 times the degree-53 rung and the largest
this family offers inside single-word arithmetic. On a 36112-point
compact-table base, 32 targets:

| | |
|---|---|
| factor base | 36112 points / 296 columns |
| precompute | 79.0 s |
| descent | 53.1 ms/target, 213 108 probes |
| signed-Frobenius ρ | 3.248 s/target, 2 263 934 steps |
| charged ρ/IC | **61.2** |
| amortised ρ/IC | **1.29** |

32 of 32 solved, 32 of 32 confirmed by the ρ walk, and **all three
verdicts true** — charged, amortised and whole-process — at 32 targets.

That last point is worth dwelling on, because the earlier projection said
otherwise. The reach table two sections up put the charged ratio at 5.1
at 48 bits and said the whole-process class would need about 315 targets
there. Measured: **61.2, and 32 targets sufficed.** The projection
under-predicted by a factor of twelve.

It was not wrong when it was written. It used the constants of a
15264-point base with a full scan, and since then every one of those
constants has moved: collection reached its `2r/|F|²` floor, two setup
costs left the precompute, and the compact table more than doubled the
base. A projection is only as good as the implementation it was measured
on, and this one was overtaken three times in an afternoon.

### Where the boundary is now

The charged ratio falls as `√r`, so from 61.2 at 48 bits the charged
advantage runs out around

```text
    48 + 2·log₂(61.2) ≈ 60 bits
```

against the 52 the earlier table gave. That is an extrapolation from one
point along a theoretical slope, and it should be read as one — the two
wide-base rungs actually measured (75.2 at 45 bits, 61.2 at 48) fall more
slowly than `√r` predicts, which would put the boundary further out
still. The slope is used because it is the one the cost laws predict, not
because two points confirmed it.

None of this touches the exponent. The descent still costs `2r/|F|²`
probes and ρ still costs `√(πr/2)/√(2n)` steps; linear still loses to a
square root, and the boundary is still a boundary. What moved is where it
sits, and it moved because a constant that looked fixed turned out not to
be — four times in a row.

## The whole thing in one law, and what measuring it cost

Four sections of this note moved constants, and one moved the quantity
the constants divide. It is worth collapsing all of it into the statement
it adds up to — and then testing that statement, which turns out to be
the part that matters.

Two facts do the work. The descent needs `2r/|F|²` probes. The factor
base is bounded by memory, `|F| = √(2M/c)` for `M` bytes at `c` bytes a
stored pair. Substitute, set equal to ρ's `√(πr/2)/√(2n)` steps, and

```text
    r_max  ∝  M²
```

— two bits of reach per **doubling** of memory. That is the clean answer,
and it is wrong, in a way that only a measurement finds.

### Testing it

The law rests on `2r/|F|²` probes at a probe cost that does not depend on
`|F|`. Both halves are checkable: run one degree at three base widths and
fit. At `K_0/F_{2^61}`, 32 targets each:

| `\|F\|` | probes/target | `2r/\|F\|²` | measured/predicted | µs a probe | charged |
|---|---|---|---|---|---|
| 9 760 | 3 068 736 | 3 419 948 | 0.90 | 0.148 | 7.20 |
| 18 544 | 545 117 | 947 354 | 0.58 | 0.214 | 27.79 |
| 36 112 | 213 109 | 249 814 | 0.85 | 0.249 | 61.16 |

Fitting across the three:

```text
    probe count      ∝ |F|^−2.03      the model says −2
    seconds a probe  ∝ |F|^+0.40      the model says 0
    descent seconds  ∝ |F|^−1.64
```

**The probe-count law is confirmed to within 1.5%.** What is not
confirmed is the part nobody wrote down: a probe is not a fixed cost. It
is a random access into a table that grows quadratically with the base,
and it gets slower as the table grows — 0.148 µs into a table of 0.2 GB,
0.249 µs into one of 3 GB. Sixty per cent dearer, across a range where
the table is comfortably in DRAM the whole time.

Carry that through and the law becomes

```text
    r_max  ∝  M^1.64
```

which is 1.64 bits a doubling, not 2. The difference compounds: at 64
bits it asks for 0.67 TiB instead of 0.25, at 72 bits for 20 TiB instead
of 4, and at 80 bits for **576 TiB instead of 64** — nine times the
memory, for the same reach.

### What the probe cost is

The `|F|^0.40` is not a property of this code. It is this machine's
memory latency curve, and once that is measured the power law turns out
to be the wrong shape entirely.

A probe is one random read of the presence filter, and the filter is
quadratic in the base: 32 MB at 9760 points, 128 MB at 18544, 512 MB at
36112. Measuring the machine directly — a pointer chase, every load
waiting on the last, over working sets from 1 MB to 4 GB:

| working set | dependent | independent |
|---|---|---|
| 32 MB | 83 ns | 8.4 ns |
| 128 MB | 115 ns | 14.7 ns |
| 512 MB | 149 ns | 25.1 ns |
| 4 GB | 236 ns | 37.0 ns |

Dependent latency is **linear in the logarithm** of the working set —
`−31.9 + 21.1·log₂(MB)`, fitting eight sizes with `r² = 0.949`. And the
descent's probe cost, fitted against its own three filter sizes, is
`26.9 + 25.2·log₂(MB)`.

**25.2 ns per doubling against 21.1 ns per doubling.** Two independent
measurements — one of the pipeline, one of the bare machine — agreeing on
the slope. That is a mechanism rather than a curve fit, and it says the
per-probe cost grows as a *logarithm* of the table, not as a power of
`|F|`. The `0.40` was a local power-law fit to a logarithm across a range
of only 16×.

### Which makes the law better than the power-law fit said

Carrying a logarithm through instead of an exponent:

```text
    descent  ∝  (r/|F|²)·(K + 2β·log₂|F|)
    r_max    ∝  M² / (C + β·log₂M)²
```

Quadratic in memory with a log-squared penalty — not `M^1.64`. A
logarithm is nearly a constant, so this sits much closer to the clean
`M²` than the power-law fit suggested:

| subgroup | if a probe were free of the table | power-law fit `\|F\|^0.40` | measured log law |
|---|---|---|---|
| 48 bits | 1.0 GiB | 0.8 GiB | 1.4 GiB |
| 64 bits | 0.25 TiB | 0.67 TiB | 0.61 TiB |
| 72 bits | 4.0 TiB | 19.6 TiB | 11.7 TiB |
| 80 bits | 64 TiB | 576 TiB | **220 TiB** |

So the honest figure at eighty bits is about 220 TiB rather than the 576
the power-law fit gave — the fit over-charged the far end by a factor of
2.6, because it extrapolated an exponent that was only ever a local
approximation to a logarithm.

### Where this stops being true

All of it holds while the table is DRAM-backed. The 21 ns a doubling is a
DRAM curve measured on DRAM; it is not a law of nature and it does not
continue across the storage boundary. A 220 TiB table is not DRAM on any
machine one buys, and the step from memory to flash is a discontinuity of
two to three orders of magnitude, not another 21 ns. What the log law
legitimately covers is the range where the table still fits in memory —
which on a large server today is a few terabytes, so around 68 to 70
bits.

Beyond that the right model is not this one, and this note does not have
it. The 80-bit row is arithmetic, not a prediction, and it is the
arithmetic of a machine that does not exist.

### What the probe cost is not

The obvious suspect for that `|F|^0.40` is the TLB. A probe is one random
access into the presence filter, the filter is quadratic in the base —
326 MB at 36112 points — and four-kilobyte pages make that eighty
thousand page-table entries for a structure walked at random. Two-megabyte
pages would cut the entries by five hundred.

So: `madvise(MADV_HUGEPAGE)` on the filter and the rests before either is
filled, and the three widths re-run. No change — 1.01, 0.94 and 1.02
times the previous cost a probe, the fitted exponents unmoved.

That is not a result, though, and it is worth being clear why.
`AnonHugePages` in `/proc/meminfo` stayed at zero for the life of a
process holding a three-gigabyte table: this kernel granted no huge pages
at all, whatever the advice, and the sysfs knob reading `madvise` did not
mean the sandbox would deliver. The measurement says the mechanism was
unavailable, not that the hypothesis was wrong. The change was reverted
rather than kept unverified, since it wanted an `unsafe` call to buy
something nothing here could demonstrate.

So the cause of the `0.40` is still open — TLB misses, cache pressure,
memory bandwidth, or some mixture. It is worth settling by whoever has a
host that grants huge pages, because the exponent is exactly what turns
64 TiB into 576 at eighty bits.

### The exponent is an upper bound, not an estimate

The 0.40 was measured across tables of 0.2 to 3 GB. Every one of them fit
in this machine's memory and was served by DRAM. The whole point of the
`M` law is to make the table bigger, and past DRAM a lookup becomes an
SSD or a network access — two to four orders of magnitude slower, not
sixty per cent. So 1.64 is what the exponent looks like *before* the
memory hierarchy has its real say, and the true large-scale figure is
lower. The 576 TiB row is a floor on the cost, not an estimate of it.

I had written `M²` into this note and into a pull request before running
the three-width test, on the strength of the algebra alone. The algebra
was right about the probes and silent about their cost, and it was the
silence that mattered. The honest form of the result is the shape, which
survives both versions:

**Reach grows with memory, sub-quadratically and by a factor the memory
hierarchy sets; precompute grows linearly in `r`.** This is a method for
many logarithms on one curve, and never for one — 17 targets to pay for
itself at 48 bits, and hundreds to thousands beyond that.

## The table was 2n times larger than the group allows — 2026-09-13

Everything above treats the pair table's size as fixed by the base:
`|F|(|F|+1)/2` sums, and the only question was how few bytes a sum could
be stored in. Sixteen went to four and a half, and the reach law that
came out of it — `r_max ∝ M²/(C + β log M)²` — took the count of pairs
as given.

It is not given. The factor base is closed under the Frobenius `π` and
under negation; that is what makes the relation columns signed orbits
rather than points, and it has been true of every base this pipeline has
built. But if `F` is closed under both, so is the set of its pair sums:

```
π(P_i + P_j) = π(P_i) + π(P_j)        −(P_i + P_j) = (−P_i) + (−P_j)
```

both again sums of two base points, because `π` is an endomorphism and
`F` contains the images. So `S = {P_i + P_j}` is a union of orbits of
`G = ⟨π, −1⟩`, a group of order `2n`, and a table that answers *is `rest`
in `S`* needs one key per orbit and not one per pair. The table had been
storing every orbit `2n` times over.

### What a canonical key costs

The orbit's key has to be computable from a point alone, and the obvious
one is the least element. It is cheaper than it looks. Negation does not
move the abscissa at all — `−(x, y) = (x, x + y)` — so the sign half of
the fold is free, and `π` acts on the abscissa as a squaring. The key is

```
canon(P) = 1 + min_k x^{2^k}
```

`n − 1` squarings and a running minimum, no memory touched. That the
packed point already encoded the sign in its low bit is why the negation
costs nothing: the two halves of a `±` pair differ in exactly that bit
and in nothing else.

### Building it without canonicalising a quadratic number of pairs

Folding is worthless if the build has to form all `|F|²/2` sums to find
out which are canonical. It does not. Walk `i` over one representative of
each signed orbit and `j` over the whole base:

> Given any pair `(P, Q)`, let `g ∈ G` carry `P` to the representative
> `R` of its own signed orbit. Then `(R, gQ)` lies in the same `G`-orbit
> of pairs and *is* enumerated.

So the enumeration is onto the orbits whatever the stabilisers are —
there are no false negatives — and it forms `|F|²/2n` sums rather than
`|F|²/2`. What it does not remove, at first, is the unordered `i ≤ j` symmetry, so
each orbit is stored twice and the fold is worth `n` rather than `2n` —
61.0 times fewer stored pairs at `n = 61`, against the 61 that argument
predicts.

That second factor is recoverable, and cheaply. A sum orbit is
enumerated once from each of its two summands' orbits, and one of the
two will do: keep the entry whose row is the **smaller** of them. The
rule is *row `α` keeps `j` only when `orbit(j) ≥ α`*, and covering
survives it — given `rest = P_a + P_b` with `orbit(a) = α ≤ β =
orbit(b)`, the `g` carrying `P_a` to `rep(α)` puts `g · rest` in row `α`
with its second summand in orbit `β ≥ α`, so `canon(rest)` is still
produced, and row `β` skips the mirror image. Ordering the base by orbit
makes "orbit at least `α`" a contiguous suffix, so a row is a slice and
no addend list is ever built.

With it the fold is worth the full `2n`: **121.9× measured**, and the
base a 4 GiB budget affords goes 330376 → **467128**.

### What it buys, and what it costs

At a 4 GiB budget and `n = 61`:

| representation | bytes a pair | base it affords | `2r/|F|²` probes |
|---|---|---|---|
| full, with summands | 16 | 23169 | 6.07e5 |
| compact | ~4.5 | 42302 | 1.82e5 |
| folded | ~4.5, `2n` times fewer | **467128** | **1.49e3** |

A base 11.0 times wider, and 122 times fewer probes a target. The
temptation is to call that a 61-fold speedup. It is not, and the run says
so plainly. Two things eat it:

- **A probe got 3.1 times dearer** — 111 ns to 340 ns. That is with the
  canonicalisation done as a rotation (below); computed as `n − 1`
  squarings it was 1125 ns and 9.8 times, and the squaring chain is what
  most of the gap was.
- **Recovery grew with the base.** The compact table does not store
  summands; it recovers them by an `|F|`-long scan on a hit. That was
  0.96 ms at `|F| = 16592` and is 12.61 ms at `|F| = 177632`, and a scan
  pays it once per witness it *finds*, not once per witness it keeps.

Measured end to end at equal memory — seconds per decomposed target,
which is the only figure immune to the fact that a scan stops at its
first witness:

**0.455 s → 0.0027 s, a factor of 168** — but only after four further
changes, and on the code as it stood it was 2.0×. What happened in
between is the rest of this section, and it is the more useful half of
the result.

### Two things the measurement had to be rescued from

The first attempt compared a folded table against a compact one small
enough to sit in cache, which flattered the baseline; the second divided
elapsed time by an assumed probe count, which is wrong in the one
direction that mattered — the folded base decomposes 99.2% of targets
against the compact base's 0.4%, so its scans stop early and never spend
the probes they were charged for. Seconds per decomposition is the
honest denominator, and it is the one the table above uses.

The early stopping also found a real inefficiency. The scan prepared its
whole row — `|F|` subtractions, and on a folded table `|F|`
canonicalisations — before probing any of it, which is linear work paid
in full however early the stop comes. On a base wide enough to decompose
most targets that setup *is* the cost. Both now run a block at a time:
long enough to amortise one field inversion and keep the squaring chains
interleaved, short enough that an early exit throws away at most a block.

### The squaring chain was not the only way to name an orbit

`π` is a squaring only in a *polynomial* basis. In a normal basis
`{β, β², β⁴, …}` it is a one-bit cyclic rotation, because

```
(Σ c_k β^{2^k})² = Σ c_k β^{2^{k+1}}
```

So the Frobenius orbit of `x` is the set of rotations of its coordinate
word, and the least rotation names it: one change of basis — an `F_2`
linear map, applied as eight byte-table lookups — and `n` rotations.
Nothing squared and nothing reduced.

The name it gives an orbit is *not* the least element of that orbit; it
is a different function of the point. That does not matter. A key has
only to be constant on orbits and distinct across them, which a bijective
linear map followed by a rotation-invariant minimum is — and a table is
built and queried with the same one. Both properties are tested directly
at degrees 13 through 61, against the squaring chain it replaces.

| | squaring chain | rotation |
|---|---|---|
| a probe | 1125 ns | **340 ns** |
| the build, 2.59e8 pairs | 116.7 s | **37.4 s** |
| the fold, end to end | 2.0× | **2.8×** |

The build gained the same factor as the probe, which it should: it
canonicalises every pair it stores.

### The recovery, and the search that threw it away

With the canonicalisation cheap, the cost moved to **summand recovery**.
The compact representation does not store which two base points made a
sum; it recovers them by a scan over the whole base, `O(|F|)` — 0.94 ms
at `|F| = 16592` and 12.88 ms at `|F| = 177632`. It grows with the very
base the fold exists to widen.

Two things were wrong with it. The smaller: `recover_pair` looked each
difference up in a `std::collections::HashMap`, whose default hasher is
SipHash — a strong hash bought for keys that are already the output of a
packing, and paid `|F|` times. An open-addressed table using the hash the
presence filter already computes, and a scan blocked so its scratch stays
in cache rather than four megabytes of it going out and coming back, took
12.88 ms to **8.66 ms**.

The larger was not in the recovery at all. A three-summand search paid it
once per witness it *found*, not once per witness it kept. The condition
`j ≤ k` is a sorted-witness condition: it is how a triple found three
times over — once for each of its indices playing the role of `k` — gets
counted once, and enumerating every sorted witness is exactly what exact
yield needs. But a descent does not want an enumeration. It wants one
decomposition, and it re-checks the sum in the group before returning it.
So it was rejecting two witnesses in three *after* paying the `O(|F|)`
scan that produced them, to arrive at the same answer later.

`decompose_fast` now takes any witness and sorts on the way out, so what
a caller sees is unchanged. Recoveries per decomposition fell from about
eleven to about one, and the scan reaches a witness about three times
sooner.

### Recovery was still `O(|F|)`, and that cancelled the whole point

At 18 ms a decomposition, about 6.7 ms was the probes and about 8.7 ms
the single recovery. That is not merely unbalanced — it *cancels the
lever the fold exists for*. Probes fall as `1/|F|²` and an `O(|F|)`
recovery rises as `|F|`, so widening the base by `√2`:

| | probes | recovery | total |
|---|---|---|---|
| `\|F\| = 177632` | 6.67 ms | 8.66 ms | 15.3 ms |
| `\|F\| × √2` | 3.34 ms | 12.2 ms | 15.6 ms |

A wash. The fold buys a base eight times wider and the recovery hands it
straight back.

### The orbit of a summand survives the fold

It need not. A folded entry was built as `canon(P_{rep(r)} + P_j)`, so a
target whose key matches it satisfies `g · target = P_{rep(r)} + P_j` for
some `g ∈ G`, and therefore

```
target = g⁻¹P_{rep(r)} + g⁻¹P_j
```

with the first summand somewhere in the **same signed orbit `r`** —
orbits being exactly what `G` preserves. So recovery never had to walk
the base. It walks the `2n` points of orbit `r`: 122 against 177632.

And the tag costs no memory, because it is spent out of the rest rather
than added to it: the stored word becomes `(orbit << 16) | rest16`. The
rest keeps sixteen bits, so a probe gets through wrongly about one time
in four thousand rather than one in `2²⁸` — and a wrong admission now
costs an orbit walk of a few microseconds instead of a scan of the base.
Over the probes a decomposition spends that is some four false positives
and tens of microseconds, against the milliseconds the shortened scan
saves. Still no false negatives; the group has the last word either way.

**Recovery: 8.66 ms → 0.037 ms.**

### A bug the tag exposed

Folded tables had been bucketing on the canonical key's high bits. A
canonical key is a *minimum over `n` rotations*, which sits far below the
middle of its range, so the whole table piled into the first few buckets.
The presence filter hid it completely — a probe that misses never reaches
a run — but every *hit* then walked an enormous one. It only became
visible when recovery stopped dominating and a single lookup on the
degenerate key took seconds.

Folded tables now bucket on the hash the filter and the rest are already
made of, which is flat. Probe 374 → 285 ns, and far more than that on
hits.

The degenerate key is worth naming: `O` is the sum of every `±` pair, so
every orbit representative stores an entry under its key, and a lookup
there is offered *every* orbit at once. Two membership tests that were
linear-per-push went quadratic on it; they are sorts now, and it costs
14 ms for the 88816 pairs it legitimately returns.

### A measurement that was not comparing what it said

Before the stage table, a correction to it. The folded base in those runs
was chosen by scaling the compact point count by `√(2n)` — the right
factor only once each sum orbit is stored *once*, which before the row
filter it was not. So the folded table held about **1.8 times the bytes**
of the compact one it was being compared against, and "at equal memory"
was not.

The example now picks the folded base by bisecting the sizing law against
the compact table's byte count, so the two match by construction: 0.64
GiB against 0.64 GiB, 137655528 stored pairs against 138074720. Every
figure in the final row below is at matched bytes; the earlier rows are
kept for the shape of the progression, and each of them flattered the
fold by that 1.8.

### What the fold is actually worth

| | s per decomposed target | the fold |
|---|---|---|
| the code as it stood | 0.229 | 2.0× |
| key as a rotation | 0.168 | 2.8× |
| faster recovery | 0.150 | 3.2× |
| any witness, not the sorted one | 0.018 | 27.8× |
| the orbit tag | 0.0027 | 168× |
| the row filter, at matched bytes | **0.0025** | **197×** |

The fold was worth 2.0× on the code as it stood. The rest is four costs
that only a base eight times wider makes visible, and three of the four
were not about the fold at all — the compact table had been paying them
quietly at every width.

What matters more than the 168 is that the probes are now 98.6% of a
decomposition. The base-widening lever the fold exists for is connected
again: a `√2` widening should be worth about 1.9× rather than nothing.

### Which is the test worth running

That last sentence was written before the run, and it is the only claim
in this section that was a prediction rather than a description. Widening
the folded base by `√2`, from 177632 to 250832:

| | `\|F\|` | s a decomposition | recovery |
|---|---|---|---|
| before | 177632 | 0.0027 | 0.037 ms |
| `× √2` | 250832 | **0.0014 / 0.0013** | 0.042 / 0.044 ms |

Two runs of the same configuration, so **1.93× and 2.04×** — call it 1.9
to 2.0 against a prediction of about 1.9. The recovery grew by the `√2`
the model says it should, and the hit rate reached 1.000: every one of
the 21779 and 22619 targets decomposed at the first scan.

(The spread is why both are here. A thirty-second window at a
millisecond-scale cost is a few thousand samples of a geometric variable,
and this note has already once mistaken two points of that noise for a
law.)

That is the result, more than the 168 is. A speedup is a number; a
restored lever is a direction. Before the orbit tag, more memory bought
nothing at all — the `O(|F|)` recovery exactly ate the `1/|F|²` fall in
probes. The `M²/(log M)²` reach law had quietly stopped applying to this
implementation, and nothing in the timings said so until the two halves
were measured apart.

The `M²/(log M)²` reach law is unchanged in shape by any of this. What
the fold moves is the constant, by putting `n` times more base behind the
same byte.

The lesson is the one this note keeps relearning. The algebra said `2n`
and it was right about `2n`; what it could not say was that at `n` times
the base, four costs nobody had been watching — a hash chosen for
strength, a sortedness condition kept for an enumeration nothing was
enumerating, a linear recovery, and a bucket index taken off bits that
were not uniform — would be the whole of the bill.

Though the algebra had the last word after all: that the orbit of a
summand survives the fold is a one-line consequence of `G` preserving
orbits, and it is what turned the linear recovery into a constant one.

## End to end at a width only the fold can hold — 2026-09-13

Everything above is a microbenchmark: a decomposition oracle timed on its
own. The pipeline has now been run whole at a width only the folded table
can hold — `K_0/F_{2^61}`, a 48-bit subgroup, **300608 base points in
2464 signed orbits** — through select, collect, logs and solve, with a
signed-Frobenius ρ baseline on the same targets.

**32 of 32 targets verified.** 73495 relations from 80000 probes; the
sparse solve filtered `73495 × 2464` to a core of 236 and finished in
0.088 s.

### The A/B that was not planned

The first attempt ran on a binary built two commits earlier, before the
orbit tag and the row filter. That was a mistake, and it turned into the
best measurement in this note: the same parameter file and the same base,
run twice on binaries differing only by those two changes.

The relation counts came back **identical, unit for unit** — 18347,
18363, 18394, 18391 — which is a stronger check that the tagged,
row-filtered table answers exactly as the untagged one than any unit
test. A whole pipeline agreeing to the relation.

| | before | after | |
|---|---|---|---|
| collect, the four units | 779.2 s | **24.8 s** | 31.4× |
| collect stage | 948.1 s | 123.2 s | 7.7× |
| solve stage | 162.2 s | **2.5 s** | 64.9× |
| precompute | 1061.6 s | 235.7 s | 4.5× |
| peak resident | 4.68 GiB | **2.40 GiB** | the row filter, visible |
| descent | 0.0381 s/target | **0.0100 s/target** | 3.8× |
| charged, against ρ | 86.7× | **330.7×** | |

The *before* column is explained by one number: collection spent 11.85 ms
a probe, which at this width is about one `O(|F|)` recovery scan per
probe — precisely the cost the orbit tag removes.

### Against the best this note had before

| | 36112 points | 300608 points |
|---|---|---|
| columns to solve | 296 | 2464 |
| precompute | 79.0 s | 235.7 s |
| descent | 0.0531 s/target | **0.0100 s/target** |
| charged, against ρ | 61.2× | **330.7×** |
| amortised over 32 targets | **1.288** | 0.45 |

Charged, the folded width is 5.4 times better. Amortised over 32 targets
it is **worse**, and that is not a footnote — the base has 8.3 times as
many points and therefore 8.3 times as many orbits to find logarithms
for, so the precompute triples. A wider base buys a cheaper descent and
charges for it once, up front.

### Most of that precompute was not needed

The run collected **73495 relations to certify 2464 columns** — a
thirtyfold oversupply — and then spent 87.6 s of the logs stage verifying
all of them, against 0.088 s for the solve itself. The unit counts had
been inherited from the 36112-point file, where a probe yielded far less.

Sized from the covering bound instead: each relation touches three
orbits, so covering every column needs about `cols·ln(cols)/3 = 6414`
relations before any margin for independence, and at the measured 0.919
relations a probe, 14000 probes should give about 12900. It gave 12843,
and **all 2464 columns certified with none rejected** — the margin is
sized rather than lucky.

| | oversupplied | sized |
|---|---|---|
| relations | 73495 from 80000 probes | 12843 from 14000 |
| logs stage | 87.6 s | **24.7 s** |
| precompute | 235.7 s | **156.9 s** |
| amortised over 32 targets | 0.45 | **0.67** |
| charged, against ρ | 330.7× | 326.0× |

The charged ratio is unchanged within noise, as it must be: sizing the
collection does not touch the descent.

*(A first attempt at this measurement gave 137.9 s for collect. The test
suite was running during that stage's table build. It was re-run with
nothing else on the machine, and those are the numbers above.)*

### Where the precompute actually goes

| | seconds |
|---|---|
| select | 25.1 |
| pair-table build | ~102 |
| probing | 4.5 |
| verifying relations | 24.7 |
| **the linear algebra** | **0.029** |

The sparse solve is two hundredths of a percent of precompute. The block
Wiedemann machinery earlier in this note — the filtering, the Krylov
sequence, all of it — solves a problem that is no longer anywhere near
the cost. What precompute is made of at this width is the pair-table
build and verifying relations in the general group arithmetic.

### The crossing

```
(156.9 − 79.0) s  ÷  (0.0531 − 0.0101) s/target  ≈  1800 targets
```

Below about 1800 targets the narrower base wins; above it the folded one
does. It was 3600 before the collection was sized. That is the cleanest
statement this work has produced of what it actually is: **a method for
many logarithms on one curve, and never for one.** The reach law says how
far memory can take you; this says how many targets you must have before
taking it is worth anything.
