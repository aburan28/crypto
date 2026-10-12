# Symmetry quotients for elliptic-curve point decomposition

Implementation note: this file preserves the original proposed study. The
initial fixed-field enumeration/correctness slice is now implemented; see
README.md and RESULTS.md. Polynomial-system and neighbor experiments below
remain proposed.

Experiment specification • Isogeny Labs • 24 September 2026

**Status:** designed, not benchmarked. The small combinatorial counts and invariant-subspace dimensions below were checked locally. No Gröbner run or neighbor construction was performed for this specification.

**Question:** Does identifying symmetry-equivalent decompositions reduce the work needed to solve a point-decomposition problem, and does a certified isogenous neighbor change that benefit?

**Plain English:** the solver may see several descriptions of the same answer. We want to keep one description without deleting genuinely different answers. Then we measure whether that makes the equations easier, including the cost of recovering the original points.

The study uses deliberately small, fully enumerable fields. Its outputs are exact solution-set comparisons and algebraic solver measurements.

## 1. Evidence and hypotheses

Three established research threads motivate the experiment:

- Faugère–Gaudry–Huot–Renault study symmetries and quasi-homogeneous structure in decomposition systems for particular curve models [1].
- Faugère–Huot–Joux–Renault–Vitse use small-order torsion to obtain more compact symmetrized summation polynomials [2].
- Galbraith–Granger–Merz–Petit study Frobenius-invariant factor bases and distinguish reducing the number of systems from reducing the cost of each system [3].

These results justify testing symmetry. They do not establish that descending an isogeny volcano lowers Gröbner regularity.

The current project report, *Koblitz-Volcano-n53-n83-Structural-Results.pdf*, was read before preparing this specification. Its small F_103 control reports equal d_reg_top across the two curves in 240 paired cases, matching transported decomposition counts in 504 cases, and a reduction from 3 to 2 after equation preprocessing for four-point bases on both curves. It reports no F4/F5 timing for that control. Those are reported prior results, not newly reproduced measurements.

Pre-register four hypotheses:

| ID | Hypothesis | Evidence that would support it |
|---|---|---|
| H1 | Permutation symmetry can reduce solver work for a fixed target. | Same exact solution orbits; lower complete solve cost. |
| H2 | An invariant-coordinate encoding can improve algebraic work beyond output deduplication. | Smaller matrices or lower observed solving degree, with lifting included. |
| H3 | Frobenius bookkeeping can avoid repeated work across equivalent targets. | Correct target transport and measured reuse in the declared target distribution. |
| H4 | A certified neighbor changes the computational benefit of symmetry. | A reproducible curve-by-encoding interaction on transported instances. |

Failure to lower degree does not refute H1: smaller matrices at the same degree may still help. Conversely, a lower degree with slower lifting is not a performance win.

## 2. The mathematical contract

Fix a finite group of rational curve points H, an allowed point set B contained in H, a target R in H, and a tuple length m. Define

    D_R(B) = {(P_1,...,P_m) in B^m : P_1 + ... + P_m = R}.

The identity point is excluded from B in the primary experiment. Repeated points are allowed; a distinct-points-only slice is separately labeled. All baseline and candidate variants use the same convention.

**Permutation symmetry.** If every slot uses the same B, any permutation of the m slots preserves D_R(B). A canonical sorted tuple represents the resulting multiset. If slots have different supports, only permutations preserving those supports are allowed.

**Point-map symmetry.** If g is a group automorphism of H with g(B)=B, then

    g(D_R(B)) = D_{g(R)}(B).

This is a symmetry of the same target only when g(R)=R. The stabilizer of R is the subgroup of maps satisfying that condition.

**Frobenius.** For a model defined over F_2, tau(x,y)=(x^2,y^2) permutes its F_(2^n)-rational points. Applying tau to every summand also applies tau to R. Therefore:

- For a fixed target, quotient only by Frobenius powers that fix R.
- For a target-orbit experiment, transform the entire pair (R, decomposition), retain the transformation, and lift the answer back to the requested R.
- Never rotate each summand independently while retaining the original target. Those are generally different sum equations.

Negation has the same issue: simultaneous negation changes R to -R unless R=-R. An x-coordinate encoding also hides point signs; it does not make independent sign choices harmless.

For a verified finite action Gamma on a finite solution set D, the audit checks

    |D/Gamma| = (1/|Gamma|) sum_(g in Gamma) |Fix_D(g)|.

This is an integer counting identity. It is not permission to divide by |Gamma| inside a field of characteristic dividing |Gamma|.

**Simple counting control:** for eight possible points and three slots, there are 512 ordered triples but only 120 multisets. Restricting to distinct points gives 336 ordered triples and 56 sets, exactly a factor of six. These counts concern the whole candidate space, not the number of solutions to a particular R. Repeated points and other stabilizers change orbit sizes, so an automatic m! speedup claim is invalid.

## 3. Toy instances and support controls

Use the two small binary models

    K_a: y^2 + xy = x^3 + a*x^2 + 1,  a in {0,1}.

They are separate families; do not treat K_0 and K_1 at the same extension degree as an isogenous pair.

Primary field sizes: n in {5,7,9,11}. Start with n=7, a=0 and m=3. Expand to m=2 and m=4 after correctness passes. Cap the actual point support at 24 points, so even the largest m=4 candidate space has only 331,776 ordered tuples.

Compare three support families:

| Support | Purpose | Required record |
|---|---|---|
| Ordinary coordinate-subspace support | Familiar algebraic control | Basis, defining equations, actual point count, closure result |
| Union of complete Frobenius point orbits | Guaranteed point-set closure | All orbit lengths, actual point count, membership-encoding cost |
| Frobenius-invariant linear-subspace support, when small enough | Closure with a potentially simple algebraic description | Exact subspace, dimension, point count, admissibility certificate |

Desired support sizes are 8, 16 and 24; these are targets, not assumed attainable values. Select the closest admissible complete-orbit size at or below each cap, omit duplicates, and make the random-set control have exactly that actual cardinality. Never trim an orbit to force a size match. A dimension-matched subspace comparison and a point-count-matched set comparison are different experiments; label them separately.

Closure alone can be expensive: replacing a small support by its full orbit closure may increase both membership cost and the number of admissible points. Charge those changes explicitly.

There is also a structural obstruction to choosing arbitrary invariant subspace dimensions. In a normal basis, Frobenius is a cyclic shift. For these odd n, X^n-1 is squarefree, and its irreducible factor degrees determine the possible dimensions of invariant F_2-linear subspaces. The following small cases were verified by exact polynomial trial division over F_2:

| n | Irreducible factor degrees of X^n-1 | Possible invariant dimensions |
|---|---|---|
| 5 | 1, 4 | 0, 1, 4, 5 |
| 7 | 1, 3, 3 | 0, 1, 3, 4, 6, 7 |
| 9 | 1, 2, 6 | 0, 1, 2, 3, 6, 7, 8, 9 |
| 11 | 1, 10 | 0, 1, 10, 11 |

Thus n=7 permits a three-dimensional invariant subspace, while n=11 does not. Record the missing n=11 case as structurally unavailable. Do not call an arbitrary three-dimensional subspace invariant or expand a support past the size cap to fill the table. The broader construction issue is discussed in [3].

For every selected B, enumerate all admissible tuples once to build an independent table of exact targets, multiplicities and permutation orbits. Keep this ground-truth table out of solver inputs and caches.

## 4. Variants to compare

| Variant | Change from the matched baseline | What it tests |
|---|---|---|
| A: ordinary encoding | Existing polynomial encoding and preprocessing | Reference cost |
| B: permutation representatives | Keep one representative of each permitted slot-permutation orbit | Whether symmetry breaking pays for itself |
| C: invariant coordinates | Express symmetric information through suitable invariants; retain complete support and lifting conditions | Whether the algebra itself becomes cheaper |
| D: Frobenius target bookkeeping | Canonicalize whole target instances and retain inverse transformations | Cross-target reuse, separately from fixed-target cost |
| E: combined | Combine only changes that passed B–D independently | Interaction and total overhead |

Output-only sorting/deduplication is a negative control for solving time: it can reduce output storage while leaving the expensive solve unchanged. Do not report that as a solver optimization.

For three scalar coordinates, the elementary symmetric quantities are

    e1 = x1+x2+x3
    e2 = x1*x2+x1*x3+x2*x3
    e3 = x1*x2*x3.

All permutations have the same values. This replaces three variables by three invariant variables; it does not automatically reduce dimension or variable count. The associated root polynomial must split into admissible coordinates in the declared field, with the correct multiplicities. Reconstruct actual curve points and verify their sum. Non-splitting roots, excluded supports, wrong signs, denominator zeros and identity-point exceptions must be handled explicitly.

Record the number and degree of invariant generators, relations between them, auxiliary variables, and lifting constraints. A general invariant ring need not be a free polynomial ring. In characteristic two, do not use an averaging projector that requires division by an even group order.

Small-order torsion symmetries are a separate optional replication of [2]. Include them only after verifying that the tuple transformation preserves the target, the support and the chosen point domain. Translation by a nonzero torsion point need not preserve a prime-order subgroup. Keep this arm out of the first minimum experiment.

## 5. Connect it to the isogeny-volcano question

Start this phase only after the single-curve symmetry experiment is correct. Select toy ordinary isogeny classes with constructible low-degree edges and certify their endomorphism orders. Count a pair as descending only after its conductor change is certified. An arbitrary adjacent curve is not automatically on a lower stratum. Standard volcano theory supplies the structural framework [4].

For a constructed F_(2^n)-isogeny phi:E->E', choose a point subgroup H on which phi is injective and verify the map and its dual/composition identity. Use two separate comparisons:

1. **Transported support:** B'=phi(B), R'=phi(R). Exact decomposition tuples correspond bijectively. Any measured difference comes from representation and computation, not a larger mathematical solution set.
2. **Native support:** choose a support by a codomain coordinate rule. Match actual point count where possible, but measure yield again: this is a different decomposition problem.

Use a small factorial comparison: source versus certified neighbor, ordinary versus symmetry-aware encoding, transported versus native support. Include several isomorphic coordinate presentations on each curve. Otherwise a favorable coordinate change can be mistaken for an effect of the endomorphism conductor.

**Cheap Frobenius is not automatically preserved.** A neighbor over F_(2^n) need not have its displayed coefficients in F_2. Squaring then maps it to its Frobenius-conjugate curve, not necessarily to itself. The 2^n-power Frobenius is the identity on its rational points and does not supply an n-fold orbit.

On phi(H), a transported permutation phi*tau*phi^(-1) exists if tau preserves H. Its existence as a permutation does not give a cheap coordinate formula or a degree-two endomorphism. Descending increases the conductor and moves to a smaller endomorphism order, so additional cheap endomorphisms must be demonstrated, not presumed. Measure the actual action cost and record which maps survive on each model.

The useful test is whether the neighbor offers a better overall equation representation after these costs. There is no predetermined direction for that outcome.

## 6. Correctness and workload design

Freeze field polynomials, bases, point serialization, support points, target list, curve models, map certificates, software versions, monomial orders, preprocessing and resource caps in a manifest before timings.

For each cell, use four independently seeded supports when available. Use all targets when |H| is small; otherwise select 64 uniformly without replacement using a fixed seed. Include unsatisfiable targets in the main workload. Also keep an explicitly labeled challenge set containing repeated summands, short point orbits, fixed targets, zero-yield targets and exceptional coordinates.

Run two distinct tasks:

- **Exact validation:** enumerate and compare the full solution sets modulo the declared action. Verify no solution orbit is omitted, merged incorrectly, or introduced. Check orbit sizes and the counting identity in Section 2.
- **Performance:** find one verified decomposition or certify that none exists for each target. Use identical stopping rules. Include root extraction and validation in the timer. Complete solution enumeration is a separate timing if desired.

Frobenius gets two workload views. The primary view uses the pre-registered uniform target sample. A positive-control batch deliberately contains complete target orbits and checks that reuse/lifting works. Report that artificial repetition explicitly; do not use it to predict savings for fresh targets in a large group. Also report one-target-per-orbit timings, where target caching has no repeat work to remove.

A deliberately non-invariant support must fail the closure check. Incorrect per-summand Frobenius normalization and incorrect target-preserving claims must fail against the exhaustive oracle. These negative controls establish that the experiment can detect its main failure modes.

## 7. Measurements and decision rules

Use the same solver backend, thread count and machine for paired timings. Randomize execution order. Repeat timings five times; use medians per instance and report the complete distributions. Statistical resampling should cluster by support seed and, for Frobenius experiments, by target orbit; five repeats of one instance are not five independent mathematical samples.

| Measurement | Meaning |
|---|---|
| Exact ordered solutions and solution orbits | Correctness and actual redundancy |
| Closure, stabilizer and orbit-length distribution | Whether a claimed symmetry is valid and how large its orbits are |
| Variables, equations, monomials and input degrees | Size of the representation presented to the solver |
| Maximum degree processed by the backend | Observed solving-degree diagnostic |
| Peak matrix rows, columns, nonzeros and memory | Size of the expensive linear algebra |
| Setup, solve, canonicalization, lifting and verification time | Complete cost, including overhead |
| Solved/unsatisfiable/timed-out/error counts | Coverage and failures, without survivor bias |
| Distinct targets and target orbits served | Denominator for any reuse claim |

Do not label the final Gröbner-basis degree, maximum matrix degree, first-fall degree or a homogeneous regularity proxy interchangeably as d_reg. Record each separately with its exact definition. In particular, reducing modulo Boolean field equations or changing homogenization can alter the reported degree statistic without making the solve faster.

Primary timing for a fixed batch is

    T_batch = T_support/setup + sum(T_encoding + T_solve + T_lift + T_verify + T_bookkeeping).

Report cold cost and amortized cost separately. The precomputed validation oracle is test infrastructure and is excluded from the solver timing; any preprocessing used by the candidate itself is included. Charge map construction and support transport to the neighbor comparison.

Report total batch time and coverage first. A secondary cost-per-verified-decomposition metric includes failed attempts in its numerator and counts at most one success per requested target; if there are no successes it is undefined/infinite. Orbit images do not count as independently solved systems. Additional images are also not automatically independent algebraic relations.

Use the following pre-registered outcomes:

- **Correctness failure:** any mismatch with exhaustive ground truth; reject the variant regardless of speed.
- **Inconclusive timing:** unequal completion coverage, substantial censoring, or too few independent instances. Record the limit instead of deleting timeouts.
- **Promising toy performance:** at least 20% lower median complete paired batch cost on two field sizes, equal completion coverage, no correctness failures, and a paired uncertainty interval supporting improvement. This is an empirical progression criterion, not a complexity claim.
- **Algebraic improvement:** lower observed solving degree and/or lower matrix work under the fixed definitions. Report separately from wall-clock improvement.
- **Neighbor effect:** replicated curve-by-encoding interaction on transported supports after charging map costs and checking coordinate controls. A single favorable pair is a lead, not evidence that deeper strata are generally easier.

Use a 60-second and 2-GiB limit per solver instance, plus a 30-minute CPU cap for the first pilot. If that pilot is censored, report it; do not silently increase the resource budget or field sizes.

## 8. First executable milestone

The first implementation should do only this bounded experiment:

1. Enumerate K_0(F_(2^7)) and construct one small Frobenius-invariant support plus a matched ordinary support.
2. Build exact m=3 solution tables, including repeated points and unsatisfiable targets.
3. Compare the ordinary encoding, permutation representatives and an invariant-coordinate encoding, each with complete lifting.
4. Audit Frobenius transport separately on full target orbits and on one representative per target orbit.
5. Publish exact correctness counts, observed solving degrees, matrix sizes and complete timings, including negative results.

Only then expand the toy matrix and compare certified neighbors. Expected deliverables are a frozen instance manifest, ground-truth certificates, solver traces, a row-level results table and a short interpretation. No solver speedup or smaller regularity is claimed by this specification.

## References

[1] J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using Symmetries in the Index Calculus for Elliptic Curves Discrete Logarithm*, Journal of Cryptology 27, 595–635 (2014; DOI registered 2013). [Institutional record](https://researchportal.ip-paris.fr/en/publications/using-symmetries-in-the-index-calculus-for-elliptic-curves-discre/). The abstract and publication record were checked; the full-text HAL endpoint was unavailable during this review.

[2] J.-C. Faugère, L. Huot, A. Joux, G. Renault, V. Vitse, *Symmetrized Summation Polynomials: Using Small Order Torsion Points to Speed Up Elliptic Curve Index Calculus*, EUROCRYPT 2014, pp. 40–57. [Publisher abstract and bibliographic record](https://link.springer.com/chapter/10.1007/978-3-642-55220-5_3). This specification relies on the abstract-level result, not an independently reproduced implementation.

[3] S. D. Galbraith, R. Granger, S.-P. Merz, C. Petit, *On Index Calculus Algorithms for Subfield Curves*, SAC 2020, revised proceedings 2021. [ePrint 2020/1315](https://eprint.iacr.org/2020/1315); [public conference manuscript](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf). Relevant parts: Frobenius symmetry, invariant supports, and the distinction between per-system cost and system count. These are research precedents, not measurements of the proposed toy experiment.

[4] A. V. Sutherland, *Isogeny volcanoes* (2013). [Author manuscript](https://arxiv.org/html/1208.5370v3). Used for the ordinary-volcano/endomorphism-order framework, not as evidence for a degree-of-regularity improvement. Binary-field maps and endomorphism orders still require their own certificates.
