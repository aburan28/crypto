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
| Factor base | Subspace bases; Frobenius-invariant bases; signed-orbit quotienting; orbit representatives; trace-zero subspace policy | Higher-genus GHS/trace-zero experiments | General trace-zero-variety factor bases beyond the Koblitz subspace policy |
| Relation collection | Semaev systems; point decomposition; Weil restriction; F4/F5 and hybrid solving; symmetrised systems | XL; single/double/multi-large-prime experiments | Collector-specific adapters to the common multi-large-prime layer |
| Filtering | Deduplication; singleton peeling; clique removal; structured merging; Markowitz pivoting; large-prime graph/hypergraph filtering | Collector-specific partial-relation combination | None at the reusable filtering layer |
| Linear algebra | Parallel block Wiedemann; block Lanczos; black-box Krylov methods; serial/Rayon/sharded/worker SpMV | CUDA SpMV worker (hardware validation required) | Performance evidence for distributed/GPU deployments |
| Individual logarithm | Randomized target decomposition; signed-Frobenius recovery; budgeted recursive scheduler with large-prime children and orbit canonicalisation | Curve-specific recursive-descent adapters | Higher-genus GHS target descent |

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

### Trace-zero constructions — Integrated for Koblitz subspace bases

The `koblitz-trace-zero` framework policy constructs a Frobenius-invariant
abscissa subspace, proves that every element is in the absolute-trace kernel,
and refuses a divisor containing a trace-one element before collection starts.
The factor-base search also detects this property and records its yield. See
[`ic_framework/plugins.rs`](../../src/cryptanalysis/ic_framework/plugins.rs)
and [`koblitz_factor_base_search.rs`](../../src/cryptanalysis/koblitz_factor_base_search.rs).
The GHS code in
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

### Symmetrised polynomial systems — Integrated

[`symmetrized_semaev.rs`](../../src/cryptanalysis/symmetrized_semaev.rs)
implements elementary-symmetric coordinates. The Koblitz implementation in
[`koblitz_symmetrised.rs`](../../src/cryptanalysis/koblitz_symmetrised.rs)
constructs the torsion-closed factor base, solves the symmetrised systems, lifts
roots back to points, and verifies their group sum. `koblitz-symmetrised` plus
the `symmetrised` oracle is selectable through the IC framework and has an
end-to-end logarithm-recovery test.

### Large-prime variants — Experimental

Three distinct research paths exist:

- residual collisions analogous to the large-prime variation:
  [`residual_walk.rs`](../../src/cryptanalysis/residual_walk.rs);
- bounded single-large-prime Jacobian experiments:
  [`hyperelliptic_index_calculus.rs`](../../src/cryptanalysis/hyperelliptic_index_calculus.rs)
  and [`hyperelliptic_ic_large_prime.rs`](../../examples/hyperelliptic_ic_large_prime.rs);
- bounded Gaudry cubic runs with configurable large-prime count and merge level:
  [`gaudry_cubic.rs`](../../src/cryptanalysis/gaudry_cubic.rs).

The common [`large_prime_filter.rs`](../../src/cryptanalysis/large_prime_filter.rs)
layer accepts all three shapes. It preserves source coefficients for replay,
peels the large-prime hypergraph, closes graph cycles, and performs bounded
Markowitz elimination on arbitrary hyperedges. Collector adapters remain
explicit so different point encodings cannot be mixed silently.

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

### Large-prime graphs/hypergraphs — Integrated reusable layer

[`large_prime_filter.rs`](../../src/cryptanalysis/large_prime_filter.rs) accepts
single, double, and arbitrary multi-large-prime relations with coefficients
modulo the subgroup order. It deduplicates, recursively peels degree-one
vertices, closes cycles, and uses fill-bounded Markowitz pivots on the
hypergraph two-core. Every complete row carries its exact sparse combination of
source relation ids, and the replay helper reconstructs the row before the
collector independently verifies its curve-group witness.

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

### Block Lanczos — Integrated

[`koblitz_sparse_la.rs`](../../src/cryptanalysis/koblitz_sparse_la.rs) provides
a selectable finite-field block-Lanczos/conjugate-direction backend. It applies
the symmetric congruence `A^T D A` only over the odd prime subgroup order,
retries deterministic random diagonals on breakdown, and accepts a result only
after checking every original equation `A x = b`. Tests cross-check it against
the dense reference and exercise selection through the full sparse solver.

### Distributed/GPU sparse matrix-vector multiplication — Integrated backend

The ECDLP CSR operator now selects serial, Rayon, deterministic row-sharded, or
external-worker products. The stable `SPMV1` contract and native host/CUDA
workers live in [`gpu/spmv`](../../gpu/spmv). External results are recomputed
and compared with the portable CPU product before use, so an accelerator can
only cause a fallback, never a false logarithm. A cluster launcher may split
the same independent row shards. No GPU or distributed speed claim is made
until matched hardware measurements and replay evidence exist.

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

### Recursive point decomposition — Integrated scheduler, experimental adapters

[`recursive_descent.rs`](../../src/cryptanalysis/recursive_descent.rs) implements
the common budgeted scheduler: memoisation, depth/node/relation caps, randomized
attempts, cycle rejection, verified relations, and a replay transcript. The
repository also contains degree, coordinate, and algebraic adapters and
experiments, including
[`coordinate_descent.rs`](../../src/cryptanalysis/coordinate_descent.rs),
[`descent_algebraic.rs`](../../src/cryptanalysis/descent_algebraic.rs), and
[`diem_descent.rs`](../../src/cryptanalysis/diem_descent.rs). Curve-specific
adapters decide which decomposition oracle supplies each recursion level.

### Large-prime descent — Integrated scheduler, experimental collectors

Large primes are ordinary recursive children in the scheduler. Their relations
can be reduced by the common hypergraph filter, while the bounded Gaudry,
residual-walk, and Jacobian collectors provide concrete experimental sources.

### Frobenius-orbit recursive descent — Integrated proof boundary

The recursive scheduler exposes canonicalisation as an explicit proof boundary:
an adapter returns a representative and the multiplier transporting its
logarithm. Signed Frobenius adapters use `±lambda^k`; adapters without a proven
action must return the identity. Every lifted relation and final logarithm is
verified independently, so an invalid per-target symmetry is rejected rather
than assumed.

## Integration order for the remaining items

1. Add collector-specific converters for the common large-prime relation type,
   retaining each collector's independent curve-group verifier.
2. Add curve-specific recursive-descent adapters and frozen termination
   transcripts for the Koblitz, Diem, and Gaudry paths.
3. Validate block Lanczos against the larger frozen block-Wiedemann corpus.
4. Validate the CUDA worker and a cluster launcher on named hardware without
   changing the `SPMV1` transcript.
5. Attempt general trace-zero or higher-genus GHS integration only after the
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
