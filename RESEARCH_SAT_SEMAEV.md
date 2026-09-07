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

The pipeline used to stall at `n = 4`–`5`.  It now decides the
`n = 19, l = 6` symmetrised-`S₄` instance from the reference corpus,
end to end, and verifies the decoded decomposition against the
original polynomial over `F_{2¹⁹}`.

Three changes did it, and they are independent of each other.

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

## Measured

`cargo run --release --example semaev_sat_bench`:

| system | params | parity | vars | clauses | rows | solve |
|---|---|---|---:|---:|---:|---|
| `S₃` | `n=5, l=5` | Tseitin | 66 | 452 | 0 | 299 µs |
| `S₃` | `n=5, l=5` | native | 32 | 66 | 5 | 82 µs |
| `S₃` | `n=7, l=7` | Tseitin | 163 | 1 246 | 0 | 4.2 ms |
| `S₃` | `n=7, l=7` | native | 63 | 147 | 7 | 451 µs |
| `S₃` | `n=12, l=6` | Tseitin | 128 | 1 001 | 0 | 619 µs |
| `S₃` | `n=12, l=6` | native | 48 | 108 | 12 | 116 µs |
| `S₄` | `n=19, l=6` | Tseitin | 2 892 | 25 074 | 0 | not reached |
| `S₄` | `n=19, l=6` | native | 767 | 2 364 | 52 | 231 s |

Tseitin expansion costs **10.6×** the clauses at `n = 19`.

### Independent agreement on instance size

The `n = 19` native row — 767 variables, 2 364 ordinary clauses, 52
parity rows — is *exactly* what the unrelated C generator in
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
cargo test --release --lib cryptanalysis::sat            # XOR engine
cargo test --release --lib cryptanalysis::binary_semaev_s4
cargo test --release --lib cryptanalysis::semaev_sat
cargo test --release --lib cryptanalysis::semaev_corpus
# the slow ones: n=15 and n=19 solves, and the exhaustive corpus audit
cargo test --release --lib cryptanalysis::semaev -- --ignored
```

Worth singling out:

- `xor_engine_agrees_with_brute_force` — 200 random dense parity
  systems, decided both by the solver and exhaustively; every verdict
  and every model checked.
- `subspace_descent_matches_field_evaluation` — the descended
  equations, evaluated at subspace points, against `S₃` computed
  directly over `F_{2ⁿ}`.
- `symmetrised_s4_vanishes_on_corpus_planted_solution` — our twelve-term
  transcription against decompositions an independent generator
  planted.
- `corpus_labels_match_exhaustive_search` (ignored) — re-derives the
  mislabelling audit in Rust.

## What's still open

- **Solve time is now the bottleneck, not encoding size.**  `n = 19`
  encodes in 11 ms and solves in ~4 minutes.  The Gauss-Jordan pass is
  rebuilt from scratch at every propagation fixpoint; an incremental
  matrix with watched pivots (CryptoMiniSat's actual design) is the
  obvious next step and should be worth a large constant factor.
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
- **Symmetry breaking.**  Nothing yet forbids permutations of
  `(X₁, X₂, X₃)`, so every solution is found `3!` times over.

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
