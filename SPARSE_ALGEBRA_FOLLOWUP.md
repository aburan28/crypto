# Sparsity-aware algebra: follow-up and applicability review

Status: applicability review and completed small synthetic Boolean-closure audit.
Date: 2026-09-28. No ECDLP speedup or blocked-workload result is reported here.

## What has been established

The existing report `degree_bounded_sparse_storage.md` (2026-09-28,
5,874 bytes; reviewed from its current saved contents) describes a standalone
Boolean/GF(2) component: on-demand monomial indexing, packed sparse rows,
incremental echelon rank and explicit fill-in budgets. It reports small
synthetic elimination checks and a separate 39,001-row storage-only check.
Those are prior reported observations, not reruns in this follow-up. The
39,001-row check did not perform elimination or reproduce the blocked workload.
No sparse-versus-dense speedup can be inferred from those figures alone.

The follow-up question is whether exact sparse polynomial computations retain
enough sparsity during elimination to reduce total memory or work. Compact
input storage, small input degree, and a low first-fall degree do not settle it.

## Literature applicability

1. Bender, Faugere and Tsigaridas, *Groebner Basis over Semigroup Algebras:
   Algorithms and Applications for Sparse Polynomial Systems* (2019):
   <https://arxiv.org/abs/1902.00208>.
   Provenance: retrieved; abstract and introduction of the primary manuscript
   reviewed in this session, together with section 2's coefficient-field
   assumption. Section 2 explicitly assumes characteristic zero. The work
   uses support/Newton-polytope structure,
   treats mixed supports, and gives results subject to regularity assumptions.
   Its solving formulation considers the algebraic torus and computes a
   saturation with respect to the product of the variables.

2. Bender, *Solving sparse polynomial systems using Groebner bases and
   resultants* (ISSAC 2022):
   <https://arxiv.org/abs/2205.09888>.
   Provenance: retrieved; primary abstract reviewed. This is a survey/tutorial
   connecting sparsity, toric geometry, resultants and Groebner computations,
   not evidence of a speedup for the present Boolean workload.

**Applicability finding (mathematical deduction).** A torus-only solution set
excludes zero coordinates. Boolean variables range over {0,1}; requiring every
coordinate to be nonzero retains only the all-ones assignment. Therefore a
torus-only solver cannot be substituted for a general Boolean solver while
claiming preservation of all solutions. Any proposed adaptation must account
for the missing coordinate-hyperplane solutions and establish the relevant
coefficient-field and regularity assumptions. The cited bounds are not yet
verified for this Boolean setting.

The existing component also uses the quotient relations x_i^2=x_i. Ordinary
polynomial monomials, squarefree Boolean monomials, and Laurent monomials have
different multiplication semantics. A representation change must preserve the
declared algebra; discarding terms beyond a degree cap changes the system in
general and must never be labelled a completed exact computation.

## General-algebra follow-up

This study concerns standalone synthetic polynomial and linear-algebra data.

- [x] Review the existing storage report and distinguish storage from elimination.
- [x] Check the primary sparse-algebra literature's stated domain.
- [x] Record the torus/Boolean solution-set mismatch.
- [x] Recover and pin the standalone component source and its original receipt
      before reporting a reproduced result.
- [ ] Define frozen synthetic inputs with sparse, block-structured and
      fill-in-heavy cases; include zero-coordinate solutions, dependent rows,
      duplicate cancellation and explicit budget-exceeded cases.
- [ ] Compare identical exact computations using the existing sparse component
      and an independent dense reference where feasible. First compare output,
      rank and algebra semantics; then compare costs.
- [ ] Record input degree distribution, input/stored/peak nonzero terms,
      intermediate fill-in, peak process RSS, rank, operation counts and elapsed
      time. Count conversion and index-construction costs. Label Python heap
      measurements separately from process RSS.
- [ ] Keep truncation, budget exhaustion and incomplete computations explicit.
      Never count a smaller or incomplete result as an improvement.
- [ ] Only after exactness and applicability checks, decide whether a
      support-aware representation has earned further investigation.

Before any performance experiment, freeze its inputs, seeds, source hashes,
reference, metrics, success/stop conditions and repetition policy in a
versioned protocol as required by AGENTS.md. The checklist above is a research
plan, not an executed experiment or an admission-ready benchmark protocol.
The earlier plan above is superseded in part by the bounded synthetic audit
below. Process RSS, blocked-workload performance and cryptographic-improvement
results remain unmeasured.

## Field-archetype arithmetic checked in this follow-up

Integer checks enumerated divisors d with 1<d<m and the least k>0 with
2^k congruent to 1 modulo m:

| m | Proper intermediate subfield degrees d | ord_m(2) |
| --- | --- | --- |
| 31 | None | 5 |
| 51 | 3, 17 | 8 |
| 53 | None | 52 |
| 83 | None | 82 |
| 131 | None | 130 |

These are exact integer checks, not curve experiments. The subfield statement
uses the finite-field subfield theorem. For prime m, the degree of each
irreducible factor of Phi_m over GF(2) is ord_m(2); this explains why the
degree-31 cyclotomic structure differs despite having no intermediate subfield.
The benchmark roles and transfer qualifications are in AGENTS.md section 8b.

## Additional structural questions raised in the follow-up

These are diagnostic hypotheses for general polynomial algebra. No corresponding
measurements on the blocked system have been obtained in this review.

| Priority | Structure | Question and limitation |
| --- | --- | --- |
| 1 | Variable-interaction graph | Do equations couple small variable groups through narrow separators, or do shared variables create a large coupled core? Few monomials do not imply small elimination width. |
| 2 | Elimination fill-in | Does low input support survive exact reductions? Record intermediate support, not only initial/final storage; a memory improvement on input rows need not survive elimination. |
| 3 | Algebraic redundancy | Are many generated rows dependent, and are low-degree dependencies already handled by the current solver? F5-style criteria are existing prior art, not a newly discovered mechanism. |
| 4 | Structural fidelity across field degrees | Keep the degree-31 split cyclotomic block distinct from the irreducible blocks at 53, 83 and 131. An observed gain must name which structure it needs. |

The strongest additional literature lead is **chordal elimination / chordal
networks**, which exploit variable interactions rather than merely counting
monomials. Cifuentes and Parrilo's papers establish benefits for suitable
structured systems; their favourable complexity statements have hypotheses and
are not universal consequences of a sparse input.

- *Exploiting chordal structure in polynomial ideals: a Groebner bases
  approach* (SIAM J. Discrete Math., 2016):
  <https://arxiv.org/abs/1411.1745>.
- *Chordal networks of polynomial ideals* (SIAM J. Applied Algebra and
  Geometry, 2017): <https://arxiv.org/abs/1604.02618>.

Provenance: both primary abstracts retrieved and reviewed on 2026-09-28.
Their application and completeness conditions still need a full-paper review
before adopting a solver. A first general-algebra study should establish the
interaction graph and observed elimination width on frozen synthetic examples.
A width from one heuristic ordering is an upper bound, not a proof of optimal
treewidth. High width or rapid fill-in is a reason to lower this lead's priority,
not proof that all sparse-algebra approaches fail.

Matrix row rank is a linear diagnostic; it does not alone certify ideal
equality or a complete polynomial solve. Keep those correctness obligations
separate when interpreting redundancy.

## Executed candidate: exact Boolean closure audit

The recovered source archive and its original four-entry SHA-256 manifest were
verified. A standalone wrapper (1 <= n <= 10) replaces lossy degree truncation
with degree-prioritized exact work and records each accepted row's origin. It
uses the full Boolean quotient, retains all degrees, and allows zero coordinates.
An independent dense-bitset verifier checks provenance, input containment and
closure under every variable. These checks certify ideal equality; budget
exhaustion is explicitly incomplete. This is a correctness oracle for small
synthetic algebra, not a scalable solver or an application integration.

The declared control suite passed nine unittest methods, including 60 seeded
random systems and deliberately broken certificates. Eight measured systems
at n=6 and n=8 produced the same independently verified ideals as an exhaustive
monomial-multiple reference using the same sparse reducer. The candidate was
slower and used more peak Python allocation memory in all eight cases: measured
median time ratios were 1.45-2.66 and allocation-peak ratios 1.18-2.42. Three
untraced timing repeats and one separate traced allocation run were used.
The reference is elementary; these tiny timings are descriptive, not robust
hardware conclusions or a comparison against a state-of-the-art solver.

The monomial cases submitted fewer rows, but queue overhead outweighed that
saving. Every case reached full degree n. Full ideal closure can require 2^n
independent rows, so the implementation does not resolve the scaling problem.
Keep the verifier and negative control; do not promote degree scheduling alone
as an established optimization. No mathematical novelty has been demonstrated.

Primary prior art reviewed: PolyBoRi's Boolean quotient arithmetic and shared
ZDD representation (<https://polybori.sourceforge.net/features.html>), and
Borderbasix's established prolongation/reduction/completion framework
(<https://www-sop.inria.fr/teams/galaad/software/bbx/borderbasis.en.html>).
This wrapper does not implement those complete systems or transfer the cited
characteristic-zero semigroup bounds.

Evidence: `boolean_sparse_adaptation.zip`, saved as a durable artifact for the
requesting user, SHA-256
`f177e883dfa8d4fd2604ebf1e0b759a48c57701a9b4c54dc6d003fa855bab52f`.
It contains source, unchanged baseline, frozen protocol, all generated inputs,
raw timings/counters, tests, an example certificate, and a file hash manifest.
The evidence archive is not checked into this repository. Recorded protocol
deviation: planted measurement inputs use seed 131+n instead of literal 131;
no seed was changed after inspecting results. Source hashes were recorded after
the run, not independently timestamped before it. Process RSS was not measured.
Accordingly this is not an admission-ready repository benchmark and supports no
cryptanalytic claim.

## Executed follow-up: sparse versus packed row representation

2026-09-29 UTC; source/input freeze commit
`aa56cbb51ff53d8df64d65bdebd8aaefeb594d51`. The prior audit now lives unchanged
in `boolean_representation_audit/prior/`; all new code, inputs, protocol,
raw compressed receipts and audit output are checked into the repository.
See [the complete result](boolean_representation_audit/README.md).

Twelve frozen synthetic systems at 6/8 polynomial variables (not field degrees)
were evaluated using sparse/packed backends under exhaustive/frontier schedules.
Certificates and deterministic work counters match between corresponding
backends. The evidence audit reconciles 504 worker samples and 1,416 verified
computations, and reproduces all 48 case/backend/schedule outputs.

The numerical all-holdout criterion passes for exhaustive scheduling and fails
for frontier scheduling. Under frontier scheduling, packed/sparse compute-time
ratios on the new chain, cycle, sparse-quadratic and dense-quadratic cases are
0.974, 0.994, 0.291 and 0.203, respectively. These are provisional stage timings:
A/A controls expose substantial virtualized-host noise, memory binding failed,
and process isolation was not established. No robust hardware rate, process
memory improvement or application-level gain is claimed.

The sparse-quadratic holdout grows from at most six input terms to 50 terms in
an intermediate reduced row. Equation-level interaction widths are 1, 2, 7, 7
for those four cases. This supports measuring fill-in and variable interactions;
it does not demonstrate a new sparse solving algorithm. The next bounded
question is whether a compact support-sharing representation can preserve
structure on a fresh small synthetic suite. Keep the existing exactness
verifier and packed baseline. m=83 and challenge-specific evidence remain
unperformed; end-to-end cryptographic costs and speedup are null.
