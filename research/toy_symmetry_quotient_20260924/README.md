# Fixed toy symmetry-quotient pilot

This is the first executable slice of the symmetry experiment: exact point
decompositions, permutation quotients, symmetric-coordinate lifting, and
Frobenius transport on `y^2+xy=x^3+1` over GF(128). It has no third-party
dependencies. The field, curve and decomposition length are fixed in code.

**Status: accounting / complete-enumeration stage diagnostic.** This is not a
Gröbner implementation, an F4/F5 benchmark, or a complete index-calculus solver.
Observed solving degree, matrix dimensions, neighbor results and full-DLP
metrics are explicitly null or `NOT_RUN`. No speedup claim follows from this PR.

## Run

From the repository root:

```sh
python3 -m unittest discover -s research/toy_symmetry_quotient_20260924 -p 'test_*.py' -v
python3 research/toy_symmetry_quotient_20260924/experiment.py --output /tmp/toy-symmetry-new-run --repetitions 3
```

The output directory must not exist; saved evidence is never overwritten.
`instances.json` freezes the exact point encoding, supports and target domain.
`certificates.json` retains every target's complete point-multiset solution set,
orbit-count certificates and a concrete incorrect-rotation counterexample.
`results.json` records raw repetitions, source and input hashes, operation
counters, exclusive stage timers, coverage and comparison ratios. The source
hash identifies the actual working-tree implementation; `source_base_commit`
identifies the checkout it was developed against, not a claim that the new
implementation already existed in that commit.

The workflow `Toy symmetry quotient` runs both the tests and a one-repetition
pilot on PRs that change this experiment. It uploads fresh evidence as an
Actions artifact without overwriting the committed three-repetition baseline.

## Boundary and matched comparison

The reference visits `B^3` ordered triples. Every explicit enumeration of all
point multisets must visit at least `C(B+2,3)` distinct multisets. That count is
the boundary for this enumeration workload only; it is not a universal PDP
or algebraic-solver lower bound. Each summary row reports tuple sum checks and
ratios to both counts. Repeated summands are retained.

The variants use identical supports and cover all 116 rational targets:

| Variant | Implementation | Interpretation |
|---|---|---|
| `ordered` | Visit every ordered triple and deduplicate verified point multisets | Reference oracle; recomputed independently of candidates |
| `permutation` | Visit each point multiset directly | Exact permutation quotient |
| `invariant_lift` | Enumerate x-multisets, form their monic root polynomials, reconstruct roots with multiplicity, then enumerate admissible y lifts | Invariant-representation and lifting correctness; **not** polynomial-system solving |

In the invariant variant, coefficients are constructed from enumerated roots.
They are not found by a Gröbner solver. Returning only x-roots would lose point
sign information, so every admissible y choice is retained and every lifted
point sum is checked. Multiple y assignments for repeated x-coordinates can
produce the same point multiset; those extra tuple checks remain in the counts.

Two fixed seeds and two holdout seeds each provide a Frobenius-closed signed
orbit support and a random point support of the same cardinality. Both
three-dimensional invariant coordinate subspaces are separate controls. This
is a one-field pilot, not independent holdout curves or a scaling experiment.
`contract.json` records why this is an equivalent bounded enumeration suite
rather than a run of the WDSat solver corpus.

## Frobenius and permutation certificates

The Frobenius audit verifies all seven simultaneous rotations for every target,
including targets with no solutions. It compares the **entire** transported
solution set with the target at `tau^k(R)` and verifies the inverse rotation.
Only powers fixing R are counted as fixed-target stabilizers. It also checks
Burnside's identity on target/multiset pairs after permutation quotienting.

There is no Frobenius solve-time row: the audit receives the completed ground
truth and is correctness infrastructure. Counting its 20 target orbits as a
reduction from 116 independent solver calls would be an unsupported claim.
No relation-rank benefit is inferred from orbit images.

Permutation certificates independently count fixed ordered tuples under each
of the six slot permutations and reconstruct ordered multiplicities from each
multiset. This detects missing repeated-point cases.

## Accounting and limits

Primary counts are **tuple sum checks**, not calibrated common operations.
Nested field and curve-operation counters are retained for inspection, but
must not be added together. `field_multiply` includes calls inside squaring and
inversion; `ec_add` includes identity calls and doublings.

Each row includes enumeration and independent point-verification time.
Enumeration includes support validation, candidate generation, invariant
construction, root recovery, y lifting, target sums and canonicalization.
Shared curve/support setup is reported once for the suite, outside per-row
medians. Ground-truth and orbit audits are test infrastructure and are excluded
from candidate times. Consequently per-row timing is a diagnostic, not a cold
complete-method metric. Uninstrumented hardware/memory and polynomial-solver
costs are not inferred from these counts.

Any coverage mismatch fails the run. The plan's 20% solver improvement gate
requires actual solver runs at two field sizes; this pilot cannot satisfy it.
The next reviewable milestone is to connect these frozen fixtures to a bounded
polynomial-system backend and measure processed degrees, matrices and complete
lifting cost. A certified isogeny-neighbor comparison comes after that.

## Evidence

The committed `results/pilot-001/` directory contains the measured reference
and both candidates, three repetitions each. See `RESULTS.md` for the single
comparison table and interpretation. The canonical scoreboard links the same
evidence and keeps full-DLP metrics unmeasured.

The research design and references are retained in `SPECIFICATION.md`.
