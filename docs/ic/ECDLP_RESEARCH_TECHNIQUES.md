# ECDLP index-calculus research technique matrix

This is the canonical inventory for the elliptic-curve index-calculus techniques
implemented, exercised, or queued in this repository. It maps a research name to
a concrete native implementation and states the boundary of that implementation.

This page is an inventory, not a performance claim. A component or stage result
does not establish an ECDLP speedup. Any comparison with Pollard rho still has to
follow [the repository's full-cost measurement rules](README.md), use the same
one-target workload and resource envelope, recover and verify the scalar, and
price every stage.

## Status vocabulary

| Status | Meaning |
|:--|:--|
| **Integrated** | Used by a native end-to-end or resumable ECDLP index-calculus path and covered by correctness checks. |
| **Experimental** | A working bounded research implementation or harness, but not a general production path. |
| **Component** | Reusable native code exists, but it is not connected to the main ECDLP pipeline for this role. |
| **Backlog** | No complete native implementation currently exists for this role. |

“Integrated” says that the path exists; it does not say it beats rho or scales to
a deployed curve.

## Stage summary

| Stage | Integrated | Experimental or component | Explicit backlog |
|:--|:--|:--|:--|
| Factor base | Subspace bases; Frobenius-invariant bases; signed-orbit quotienting; orbit representatives | Trace-zero tests and subspaces | General trace-zero-variety factor-base pipeline |
| Relation collection | Semaev systems; point decomposition; Weil restriction; F4 and hybrid F4/F5 | XL; symmetrised systems; single/double-large-prime experiments | Mainline symmetrised oracle; reusable multi-large-prime collector |
| Filtering | Deduplication; singleton peeling; clique removal; structured merging; Markowitz-style pivoting | Partial-relation combination in bounded IC experiments | General large-prime graph/hypergraph filter |
| Linear algebra | Parallel block Wiedemann; black-box sparse Krylov products | Sequential Wiedemann component | Block Lanczos in the ECDLP pipeline; distributed/GPU sparse matrix-vector backend |
| Individual logarithm | Randomized one-relation target decomposition; signed-Frobenius-orbit lookup and certified recovery | Recursive decomposition and large-prime descent experiments | General recursive descent scheduler; independently valid Frobenius-orbit recursive descent |

## 1. Factor-base construction

### Subspace factor bases — Integrated

[`FactorBaseSpec::StandardSubspace`](../../src/cryptanalysis/koblitz_factor_base_search.rs)
and the divisor-family search provide reproducible linear-subspace factor
bases. The Diem demonstration independently provides the base-field
construction for extension-field curves in
[`diem_descent.rs`](../../src/cryptanalysis/diem_descent.rs).

### Frobenius-invariant bases — Integrated

[`koblitz_index_calculus.rs`](../../src/cryptanalysis/koblitz_index_calculus.rs)
constructs kernels of linearised polynomials obtained from factors of
\(x^n-1\), as well as Frobenius closures of seed spaces.
[`koblitz_factor_base_search.rs`](../../src/cryptanalysis/koblitz_factor_base_search.rs)
searches the divisor lattice and nonlinear Frobenius unions using measured
coverage and column counts.

### Quotient by \(P\mapsto -P\) and orbit representatives — Integrated

The Koblitz pipeline stores one unknown per signed Frobenius orbit. A point
\(\pi^k(Q)\) contributes \(\lambda^k\) to the representative column, while
negation changes the coefficient sign. The factor-base recipes preserve the
canonical orbit representatives needed to replay the same column assignment.

This is a representation quotient, not permission to drop relation witnesses:
every accepted decomposition is still re-added and verified in the curve group.

### Trace-zero constructions — Component

The factor-base search detects when an invariant abscissa subspace is contained
in the trace kernel and records the resulting yield property. The GHS code in
[`ghs_descent.rs`](../../src/cryptanalysis/ghs_descent.rs) also provides an
end-to-end \(m=1\) trace descent and structural \(m=2\) machinery.

The repository does **not** yet contain a general trace-zero-variety ECDLP
pipeline. In particular, higher-genus GHS descent still needs the explicit
smooth model before its Jacobian path is end to end.

## 2. Relation collection

### Semaev summation polynomials — Integrated

- Prime-field \(S_3\) and recursively constructed \(S_4\):
  [`ec_index_calculus.rs`](../../src/cryptanalysis/ec_index_calculus.rs).
- Higher and binary-field systems:
  [`semaev_higher.rs`](../../src/cryptanalysis/semaev_higher.rs),
  [`binary_semaev.rs`](../../src/cryptanalysis/binary_semaev.rs), and
  [`binary_semaev_s4.rs`](../../src/cryptanalysis/binary_semaev_s4.rs).
- Pair-and-root \(S_4\) decomposition over binary subspaces:
  [`semaev_decomp.rs`](../../src/cryptanalysis/semaev_decomp.rs).

### Point-decomposition algorithms — Integrated

The Koblitz pipeline exposes enumeration, meet-in-the-middle pair tables,
Gröbner, SAT, WDSat, and MQ-FES decomposition strategies. Returned witnesses are
never trusted directly: the points are added again and checked against the
target before the relation is accepted.

### Weil restriction — Integrated

The Koblitz Gröbner and SAT paths descend extension-field Semaev equations to
Boolean systems over the selected subspace. The symbolic Petit–Quisquater
construction is in
[`pq_descent_symbolic.rs`](../../src/cryptanalysis/pq_descent_symbolic.rs);
the bounded degree and cost measurements are in
[`ic_descent_degrees.rs`](../../src/cryptanalysis/ic_descent_degrees.rs).

### Gröbner F4/F5, XL, and hybrid guessing

- **Integrated:** native F4 and the Koblitz hybrid matrix-F4,
  matrix-F5, and inherited-F4 engines are registered through
  [`ic_framework/solvers.rs`](../../src/cryptanalysis/ic_framework/solvers.rs).
- **Component:** Boolean XL is implemented in
  [`pq_xl.rs`](../../src/cryptanalysis/pq_xl.rs) and registered as a bounded
  reference solver. It is deliberately capped where its matrix becomes
  impractical.
- **Integrated as hybrid solvers:** the matrix engines combine bounded
  Macaulay reduction, propagation, and splitting/guessing. Budget exhaustion is
  reported as inconclusive, never as “no decomposition.”

### Symmetrised polynomial systems — Component

[`symmetrized_semaev.rs`](../../src/cryptanalysis/symmetrized_semaev.rs)
implements elementary-symmetric coordinates and symmetrisation for the available
Semaev polynomials. It is useful for formula and density experiments, but it is
not yet a selectable mainline decomposition oracle. Integrating it requires a
bijection-preserving witness lift back to factor-base points and an end-to-end
cost comparison.

### Large-prime variants — Experimental

Three distinct research paths exist:

- residual collisions analogous to the large-prime variation:
  [`residual_walk.rs`](../../src/cryptanalysis/residual_walk.rs);
- bounded single-large-prime Jacobian experiments:
  [`hyperelliptic_index_calculus.rs`](../../src/cryptanalysis/hyperelliptic_index_calculus.rs)
  and [`hyperelliptic_ic_large_prime.rs`](../../examples/hyperelliptic_ic_large_prime.rs);
- bounded Gaudry cubic runs with configurable large-prime count and merge level:
  [`gaudry_cubic.rs`](../../src/cryptanalysis/gaudry_cubic.rs).

These are not silently treated as one general ECDLP large-prime collector.

## 3. Relation filtering

### Duplicate removal and singleton peeling — Integrated

[`koblitz_sparse_la.rs`](../../src/cryptanalysis/koblitz_sparse_la.rs)
deduplicates rows, removes singleton columns with their rows, and records
back-substitution data so every eliminated logarithm can be reconstructed.

### Cycle/clique selection and peeling — Integrated for the main sparse matrix

Surplus rows are removed with a clique rule designed to cascade through
weight-two columns. This is matrix filtering after complete relations have been
formed. It should not be confused with a general graph or hypergraph engine for
combining arbitrary multi-large-prime partial relations.

### Structured merging and Markowitz-style elimination — Integrated

The sparse logarithm solver merges light columns under a fill-in bound.
The generic research framework also provides a lightest-column,
Markowitz-style incremental solver in
[`ic_framework/linalg.rs`](../../src/cryptanalysis/ic_framework/linalg.rs).

### Large-prime graphs/hypergraphs — Backlog as a reusable layer

Some experimental collectors combine their own partial relations, but the
repository has no common graph/hypergraph filter that accepts single, double,
and multi-large-prime relations from every ECDLP collector. A complete version
must preserve edge provenance, reject trivial cycles, emit independently
verified complete relations, and report peeling/merge work separately.

## 4. Sparse linear algebra

### Parallel block Wiedemann and black-box Krylov — Integrated

[`koblitz_sparse_la.rs`](../../src/cryptanalysis/koblitz_sparse_la.rs)
implements block Wiedemann over the filtered relation core. It forms a block
Krylov sequence using only sparse matrix-times-block products, computes a
shifted minimal approximant basis, reconstructs eliminated variables, checks
every original relation, and leaves the final group certification to the
caller. Sparse products run in parallel over rows.

### Sequential Wiedemann — Component

[`pq_wiedemann.rs`](../../src/cryptanalysis/pq_wiedemann.rs) provides the
scalar black-box algorithm and sparse matrix-vector primitives. The main
Koblitz path uses the block implementation instead.

### Block Lanczos — Backlog for ECDLP

Block Lanczos is documented as an alternative but is not registered as an
ECDLP relation-matrix backend. An implementation must avoid unsound plain
normal-equation reasoning in characteristic two, expose deterministic counted
work, and cross-check its kernel or solution against the existing dense and
block-Wiedemann references.

### Distributed/GPU sparse matrix-vector multiplication — Backlog

The current block-Wiedemann products are shared-memory CPU code. GPU algebra in
other repository modules does not constitute an ECDLP sparse-matrix backend.
A future backend needs a stable sparse format, deterministic modular
accumulation, sharded checkpoint/restart, exact replay certificates, and
CPU-reference cross-checks before performance measurement.

## 5. Individual logarithm

### Randomized target decomposition — Integrated

The `ic solve` path samples
\(R=[a]G+[b]Q\) until one decomposition is found, substitutes certified
factor-base logarithms, recovers \(d=\log_G Q\), and accepts it only after
checking \([d]G=Q\). See the factor-base-log and descent section of
[the IC framework guide](README.md).

### Signed-Frobenius-orbit recovery — Integrated

For Frobenius-invariant bases, descent relations are rewritten into the signed
orbit columns used by the precomputed logarithm database. This is an orbit-aware
one-relation descent, not a claim that every target decomposition has an
independent Frobenius symmetry.

### Recursive point decomposition — Experimental

The repository contains degree, coordinate, and algebraic descent experiments,
including
[`coordinate_descent.rs`](../../src/cryptanalysis/coordinate_descent.rs),
[`descent_algebraic.rs`](../../src/cryptanalysis/descent_algebraic.rs), and
[`diem_descent.rs`](../../src/cryptanalysis/diem_descent.rs). They do not yet
form one general recursive scheduler with a shared cost ledger and termination
certificate.

### Large-prime descent — Experimental

Large-prime combination exists in bounded Gaudry, residual-walk, and Jacobian
experiments, but it is not a selectable target-descent policy in the main
Koblitz workflow.

### Frobenius-orbit recursive descent — Backlog unless separately proved

Frobenius quotienting is valid for factor-base columns because the action and
its scalar \(\lambda\) are known. A per-target symmetry of a point-decomposition
system does not follow automatically. Any recursive orbit descent must state
the applicable curve and subspace hypotheses, prove witness transport, and
verify each lifted relation before it can be integrated.

## Integration order for the remaining items

1. Add a reusable partial-relation graph/hypergraph layer with provenance and
   complete-relation verification.
2. Add a recursive target-descent scheduler over the existing decomposition
   oracle interface, with explicit node and total-work budgets.
3. Wire symmetrised systems into that interface with a certified witness lift.
4. Add block Lanczos as a relation-matrix backend and cross-check it on the
   frozen dense/block-Wiedemann corpus.
5. Separate sparse matrix-times-block behind a backend trait, then add
   distributed CPU and GPU implementations without changing the transcript.
6. Attempt general trace-zero or higher-genus GHS integration only after the
   explicit smooth-model map is available and independently checked.

Each item is a separate research iteration. Before measuring, it requires a
versioned protocol with frozen inputs, a matched baseline, a derived boundary,
success and stop conditions, exclusive cost accounting, and verified
one-target recovery. Stage-only results remain labelled stage diagnostics.

## References already represented by the code

- Semaev, summation polynomials (2004).
- Faugère–Perret–Petit–Renault and Petit–Quisquater, Weil-descent systems
  and algebraic solving (2012).
- Faugère–Gaudry–Huot–Renault, symmetry in elliptic-curve index calculus
  (2014).
- Gaudry and Diem, extension-field and abelian-variety index calculus.
- Galbraith–Granger–Merz–Petit, Frobenius-invariant factor bases for
  subfield curves (2020).
- Wiedemann, Coppersmith, Kaltofen, and Montgomery, sparse black-box
  linear algebra.
