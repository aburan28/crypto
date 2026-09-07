# SAT-encoded binary Semaev systems

**Modules:** `src/cryptanalysis/semaev_sat.rs` (encoders),
`src/cryptanalysis/binary_semaev_s4.rs` (symmetrised `S₄` + descent),
`src/cryptanalysis/semaev_corpus.rs` (reference corpus),
`src/cryptanalysis/sat.rs` (CDCL + native XOR reasoning)
**Bench:** `cargo run --release --example semaev_sat_bench`
**Provenance:** The "Hybrid SAT / algebraic" entry in the ECDLP
research-direction map.  Rebuilt against the prior art reviewed in
[`RESEARCH_TRIMOSKA_BENCHMARKS.md`](./RESEARCH_TRIMOSKA_BENCHMARKS.md).

## Where this got to

The pipeline used to stall at `n = 4`–`5`.  It now decides every
satisfiable `n = 19, l = 6` instance in the reference corpus — eleven
of them, in 50 s total — verifying each decoded decomposition against
the original polynomial over `F_{2¹⁹}`.

Three modelling changes made it possible; five further changes made it
about 50× faster (see [Making it fast](#making-it-fast)).

### 1. Confine the unknowns to a factor base

`weil_descend_s3` gave `2n` bit-variables because `X₁, X₂` ranged over
all of `F_{2ⁿ}` — which is not index calculus, it is "solve Semaev over
the whole field".  `weil_descend_s3_subspace(n, l, …)` confines each
`Xᵢ` to the `l`-dimensional subspace `⟨1, z, …, z^{l−1}⟩`, the standard
Gaudry/Diem factor base.

The variable count matters less than the *shape*: at `l ≈ n/m` the
system is `n` equations in `m·l ≈ n` unknowns — determined, and
something unit propagation can bite on — instead of `n` equations in
`2n` unknowns with `2ⁿ` solutions scattered through a dense quadratic
system.

### 2. Native parity constraints

`Solver::add_xor` installs a parity row outside the CNF;
`Solver::propagate` alternates watched-literal propagation with
Gauss-Jordan elimination over those rows until neither derives
anything new.  A row reduced to one unknown propagates; an inconsistent
row is a conflict.

The subtlety is the reason clause.  Each row carries the *full*
variable set of its combination alongside its original right-hand side,
so the value implied by the current trail is
`rhs ⊕ parity(mask ∩ assigned-true)`.  Both parts XOR when rows
combine, so that relation survives elimination and the reason clause
can be read straight off a reduced row: every assigned variable appears
negated-as-assigned, leaving the clause false except for the implied
literal.  `analyze()` and `backjump()` needed no changes at all.

Why it matters: refuting a dense parity constraint by resolution takes
exponentially many steps (Urquhart 1987).  Without this the solver has
to rediscover linear algebra, one conflict at a time.

### 3. Symmetrised `S₄`

`m = 2` has no asymptotic payoff over Pollard ρ; `m = 3` is where index
calculus starts to mean something.  `binary_semaev_s4` builds the
Faugère–Gaudry–Huot–Renault symmetrisation over `e₁, e₂, e₃`, which
drops the descended system from cubic to **quadratic** — the cubic part
moves into a correspondence system tying each `eᵢ` back to a symmetric
function of the `Xᵢ`.

Squaring there is pure relocation, not multiplication: in
characteristic 2 there are no cross terms, and each coefficient is an
ANF over Boolean variables where `a² = a`.  Of the twelve terms of the
descended `f₃` — including `e₁⁴, e₂⁴, e₃⁴, e₃³` — only four need a real
multiplication.

## Making it fast

The first working version decided one `n = 19` instance in **231 s**.
The whole eleven-instance `n19l6` family now takes **22.8 s**, median
**2.0 s** — about **115× per instance**.  Six changes, in the order
they mattered:

### 1. Branch only on the free variables

The single largest win, and not a solver micro-optimisation but a
modelling one.  Of the 767 variables in an `n = 19` instance, only the
**18 x-bits are free**: every monomial auxiliary is an AND of others,
and every `e`-variable is fixed by a correspondence row.  Deciding one
of those is case-splitting on something propagation already knew.

`Solver::set_branch_priority` orders the branching heap by
`(priority class, activity)`, so the free variables are exhausted
before any other is ever chosen.  The search tree goes from `2^767` to
`2^18`.  It is an *order*, not a restriction — name an incomplete set
and you get a slower solve, not a wrong answer.

### 2. Learnt-clause minimization

Ordinary CDCL hygiene that turns out to matter far more here than
usual.  A parity row's reason clause names *every* assigned variable of
its combined mask, and Gauss-Jordan makes those masks dense — so
unminimized learnt clauses ran to hundreds of literals and then had to
be walked on every propagation.  Local (self-subsuming) minimization
cut `n15l5` from 1.1 s to 95 ms on the instance it was first measured
on.

### 3. Symmetry breaking

Nothing forbade permutations of `(X₁, X₂, X₃)`, so the search
rediscovered every solution `3! = 6` times.  `S4Options::break_symmetry`
adds lexicographic `X₁ ≤ X₂ ≤ X₃` constraints — `2l + 1` clauses per
adjacent pair, using an "equal so far" auxiliary per bit position.
Worth a steady ~4.5× across all three families.

### 4. No more clause-database growth from parity reasoning

Every parity implication used to push a reason clause onto the clause
database and never reclaim it — millions of clauses of hundreds of
literals over a long solve.  Reasons now live in a per-variable buffer
rewritten in place (`Reason::XorPropagated`), and learnt clauses are
periodically forgotten (`Solver::reduce_db`, detaching from the watch
lists so every stored index stays valid).

### 5. An incremental Gauss-Jordan matrix

The elimination pass used to be redone from scratch at every
propagation fixpoint, `O(rows² × words)` a time.  It is now carried
across propagations, and only rows whose pivot has since been assigned
are re-pivoted.

What makes per-row reasoning safe is one invariant:

> every row has a **pivot** variable that is unassigned and occurs in
> no other row.

Restricted to the unassigned columns the matrix is then `[I | B]` after
permutation, so a sum of `k` rows still contains all `k` of their
distinct pivots.  A combination can be unit only for `k = 1`, and
inconsistent only if a single row is already fully assigned — exactly
what a per-row check finds.  Nothing is lost by not re-eliminating.

Two consequences fall out.  Row operations are algebraically valid
whatever the trail says, so they are never undone; and backjumping only
*unassigns* variables, which cannot break the invariant, so it costs no
work at all.  The matrix is rebuilt from the original rows at restarts,
which bounds the fill-in that row operations accumulate.

Worth **2.4×** at `n = 19` and growing with the instance — the removed
cost is quadratic in the row count.  `incremental_matrix_matches_full_elimination`
checks the completeness claim directly against a from-scratch
elimination oracle, over 300 random systems under random partial
assignments, rather than trusting the argument.

### 6. Data-structure work

A position-tracked activity heap instead of a linear scan over every
variable per decision; `analyze` scratch reused rather than reallocated
per conflict; a bitset mirror of the assignment so a parity row's
read-off is `mask & !assigned` plus a popcount — `O(words)` instead of
`O(set bits)` with a two-byte lookup each.

## Measured

`cargo run --release --example semaev_sat_bench`.  Aggregates over
**every satisfiable instance** of each corpus family: per-instance
times swing by an order of magnitude on search-trajectory luck, so a
single instance measures nothing.  (An early single-instance reading
suggested symmetry breaking *hurt*; across ten instances it is a clean
3× win.  Worth remembering before drawing a conclusion from one
number.)

| family | symmetry | instances | total | median | conflicts |
|---|---|---:|---:|---:|---:|
| `n15l5` | off | 10 | 5.5 s | 582 ms | 93 593 |
| `n15l5` | on | 10 | 1.8 s | 192 ms | 33 974 |
| `n17l6` | off | 10 | 67.0 s | 6.0 s | 467 663 |
| `n17l6` | on | 10 | 13.9 s | 1.4 s | 156 612 |
| `n19l6` | off | 11 | 209.8 s | 8.5 s | 1 012 122 |
| `n19l6` | on | 11 | 22.8 s | 2.0 s | 253 496 |

`n19l6` has eleven satisfiable instances rather than ten because one
upstream `-U` instance is misannotated; the corpus records that, so it
falls into the SAT set automatically.

### Where the time went

Median solve on an `n = 19` instance, one change at a time:

| after | median | vs. previous |
|---|---:|---:|
| first working version | 231 s | — |
| free-variable branching, heaps, reused scratch | ~21 s | 11× |
| learnt-clause minimization | — | folded in below |
| symmetry breaking | 4.7 s | 4.5× |
| incremental Gauss matrix | **2.0 s** | 2.4× |

Encoding size, `n = 19, l = 6`, symmetry breaking off:

| encoding | vars | clauses | parity rows |
|---|---:|---:|---:|
| native | 767 | 2 364 | 52 |
| Tseitin CNF | 2 892 | 25 074 | 0 |

Tseitin expansion costs **10.6×** the clauses.

### Independent agreement on instance size

The native row — 767 variables, 2 364 ordinary clauses, 52 parity rows
— is *exactly* what the unrelated C generator in
`mtrimoska/EC-Index-Calculus-Benchmarks` emits for the same family
(`p cnf 767 2416` with 52 `x`-lines; 2416 − 52 = 2364).  Two
implementations written from the algebra rather than from each other,
landing on the same instance.  Locked in by
`s4_encoding_size_matches_upstream_generator`.

## Reference corpus

`semaev_corpus.rs` carries the parameters of all 60 upstream instances
— field, target x-coordinate, and the planted decomposition — but not
their generated files: upstream is GPL-3.0 and this crate declares no
licence, so regenerating from parameters avoids a licensing decision
and exercises more of our own pipeline besides.

All 30 planted decompositions satisfy our symmetrised `S₄`.  Note that
`n19l6-19-U` is **satisfiable** despite its label: upstream's `-U`
instances come from a random target and the label is an expected
outcome, not a certificate.  The table records `labelled_sat` and
`truly_sat` separately.

## Tests

```bash
cargo test --release --lib cryptanalysis::sat            # solver + XOR engine
cargo test --release --lib cryptanalysis::binary_semaev_s4
cargo test --release --lib cryptanalysis::semaev_sat
cargo test --release --lib cryptanalysis::semaev_corpus
# still slow: n=19 solve, and the exhaustive corpus audit
cargo test --release --lib cryptanalysis::semaev -- --ignored
```

The `n = 15` round trip and the corpus round trip used to take minutes
and were `#[ignore]`d; the whole default suite for these modules now
runs in under half a second.

Worth singling out:

- `xor_engine_agrees_with_brute_force` — 200 random dense parity
  systems, decided both by the solver and exhaustively; every verdict
  and every model checked.
- `aggressive_clause_reduction_preserves_answers` — 60 mixed CNF+parity
  instances solved twice, once with a forgetting budget of 2 and once
  with a budget that never forgets; the verdicts must match and the
  model from the forgetting run must still satisfy everything.
- `lex_le_accepts_exactly_the_ordered_pairs` — the symmetry-breaking
  encoder over all 256 four-bit pairs.
- `symmetry_breaking_preserves_satisfiability` — the constraint may
  only remove duplicate solutions, never the last one.
- `subspace_descent_matches_field_evaluation` — the descended
  equations, evaluated at subspace points, against `S₃` computed
  directly over `F_{2ⁿ}`.
- `symmetrised_s4_vanishes_on_corpus_planted_solution` — our
  twelve-term transcription against decompositions an independent
  generator planted.
- `s4_encoding_size_matches_upstream_generator` — the clause-exact
  agreement described above.
- `incremental_matrix_matches_full_elimination` — 300 random parity
  systems under random partial assignments; every implication the
  incremental matrix derives is compared against an oracle that
  re-eliminates from scratch.  This is the test that guards the
  completeness argument above, and it is worth having: a gap there
  would surface as a missed propagation, not a crash.
- `corpus_labels_match_exhaustive_search` (ignored) — re-derives the
  mislabelling audit in Rust.

## What's still open

- **Re-pivoting is still eager.**  When a pivot is assigned the row
  immediately looks for a replacement, where CryptoMiniSat watches two
  unassigned variables per row and only acts when *both* are gone.
  Watching would cut the remaining per-pass `O(rows)` scan to work
  proportional to the rows actually touched.
- **Fill-in is only bounded at restarts.**  Row operations make the
  matrix denser over time, which lengthens the reason clauses read off
  it.  Choosing pivots to minimise fill-in (the Markowitz rule from
  sparse LU) would keep rows shorter than picking the lowest-numbered
  unassigned variable does.
- **`S₄` is specialised to `b = 1`** — the Koblitz curve
  `y² + xy = x³ + x² + 1`.  `S₄` does not depend on `a₂`, but it does
  depend on `b`, and the `b`-powers are folded into the twelve
  constants.  General `b` needs re-deriving from
  `Res_X(S₃(X₁,X₂,X), S₃(X₃,x_R,X))` followed by re-symmetrisation;
  `weil_descend_s4` panics rather than returning a wrong system.
- **Bit-sliced coefficients.**  `AnfPoly` is a `BTreeSet` of monomials.
  A bit-packed representation over a ranked monomial space (as upstream
  uses) would be far denser and make relocation a shift.  Deferred
  deliberately: it buys descent time, and descent is 11 ms against a
  4-minute solve.
- **UNSAT is still out of reach.**  Every result above is a
  *satisfiable* instance, where the search stops at the first witness.
  Refuting a `-U` instance means exhausting the space, which the solver
  does not finish in reasonable time at any corpus size — so
  `corpus_round_trips_through_the_solver` asserts only the SAT
  direction.  Symmetry breaking is now in place, which was one of the
  two prerequisites; the incremental Gauss matrix is the other.

## References

- **M. Soos, K. Nohl, C. Castelluccia**, *Extending SAT solvers to
  cryptographic problems*, SAT 2009 — native XOR reasoning.
- **A. Urquhart**, *Hard examples for resolution*, JACM 1987 — why
  parity constraints must not go through CNF.
- **J.-C. Faugère, P. Gaudry, L. Huot, G. Renault**, *Using symmetries
  in the index calculus for elliptic curves discrete logarithm*,
  J. Cryptology 2014.
- **P. Gaudry**, *Index calculus for abelian varieties of small
  dimension and the elliptic curve discrete logarithm problem*, 2009.
- **G. Tseitin**, *On the complexity of derivation in propositional
  calculus*, 1968.
