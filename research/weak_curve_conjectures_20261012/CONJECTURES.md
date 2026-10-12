# Candidate families and unexplained properties

This ledger proposes six mathematical conjectures and one statistical
research hypothesis. It also preserves a refuted screening hypothesis.
The proposed contributions concern exact class support, intersections and
primitive higher-degree geometry. Published covering constructions keep
their original attribution; an exact statement absent from the inspected
sources has **novelty unresolved**, not an established priority claim.

`WCQ1` identifies a hypothesis, not an exact curve, measured IC candidate or
verified weakness. [conjectures.json](conjectures.json) is the machine-readable
ledger. [METHODOLOGY.md](METHODOLOGY.md) specifies the testing and evidence
contract; [LITERATURE.md](LITERATURE.md) records the bounded source audit.

## WCQ1:boundary-strip — rejected finite screen

**Statement tested:** every ordinary signed trace `t=2 mod 4` absent from both
current hyperelliptic branches, in the seven complete existing censuses,
satisfies `2*p^3-|t| <= 64*p`.

**Result:** refuted. At `p=7,t=±22`, both branch labels are zero and the gap is
664, exceeding 448. Across 87,486 eligible rows, 1,184 union-zero rows violate
the strip. All violations are retained in [FALSIFIERS.tsv](FALSIFIERS.tsv).
Interior gaps must remain visible in every later explanation.

## WCQ1:gap-count — an asymptotic class-support law

**Domain:** primes `p>3`, `Q=p^6`, all ordinary signed traces `t=2 mod 4`,
and precisely the union of the cubic and nonsplit-quadratic hyperelliptic
branches in the existing theorem. Let `Z(p)` count classes with neither
branch and `N(p)=p^2*(p-1)` count the eligible classes.

**Conjecture:** for every real `epsilon>0`, there exist constants `C_epsilon`
and `p_epsilon` such that `Z(p)<=C_epsilon*p^(1+epsilon)` for all eligible
primes above `p_epsilon`. This would give conditional class density tending
to one, with zero fraction `O(p^(-2+epsilon))`.

**Evidence:** the observed zero counts are 46, 138, 182, 270, 290, 494 and 654
at p=7,11,13,17,19,29,37. They are descriptive support, not an asymptotic
bound. Published model-density results use a different denominator.

**First test:** use p=53,59,61 as retrospective extensions because the
filename audit found existing cubic census artifacts there. Their labels
were not read in this investigation. The prospective proposed primes are
p=71,73,79, for which no matching artifact was found in the inspected ISO-1
archive. Freeze their complete evaluation censuses and the separate finite
prediction `Z(p)<=24*p` before opening labels. Do not relocate a test prime
after observing its labels. A violation rejects this
finite prediction, not the asymptotic conjecture. Refuting the latter needs
a superlinear lower bound on an unbounded sequence, with a fixed positive
exponent margin.

**Proof target:** an upper bound on the exact torus/CM intersections' empty
trace classes. The failed boundary-strip assertion cannot supply that bound.

## WCQ1:cm-mass — global structure in the missing classes

**Hypothesis type:** statistical occupancy model, not an exact classifier.
For a branch B and a fixed declared trace stratum, define
`H_B(t)=sum h(O_c)` over its selected CM orders: `c|f_pi/4` for the cubic
branch, or `c|f_pi` with `v2(c)=v2(f_pi)` for the quadratic branch. Restrict
this pilot to `D_K<-4`; exceptional units need their own stratum. Let `M_B`
be the total branch parameter-orbit occurrences assigned to the same
stratum, and `mu_B(t)=M_B*H_B(t)/sum_u H_B(u)`.

**Prediction to investigate:** the observed zero fraction in fixed `mu_B`
bins approximately follows `exp(-mu_B)`, with a separately measured residual
dependence on marked CM/Frobenius orbit structure. Parameter occurrences,
distinct invariants, signed models and class support require different
counts; choose one occurrence law before computing `M_B`.

**Distinctive property:** classes with the same simple local signature can
have different support. A global occupancy/obstruction explanation could
organize those differences while retaining exact support certificates.

**Falsification protocol:** fit no coefficient using test labels. Freeze bins
`[0,1),[1,2),[2,4),[4,8),[8,infinity)` using the training corpus; pool bins
with fewer than 30 test classes only under a preregistered rule. The finite
calibration forecast requires the absolute observed-versus-predicted zero
fraction difference to be at most 1/5 in each unpooled bin with at least 30
classes. Report every discrepancy, prime-level summaries and the
conductor-only baseline. A violating bin rejects this finite forecast.
Cross-prime uncertainty must retain dependencies within a census; an
independent-binomial interpretation is not assumed. Calibration failure
does not disprove the exact CM intersection theorem.

**Novelty boundary:** CM mass and class numbers are established objects.
This branch-selected, class-labelled prediction is the proposed new study;
its exact relationship to prior distribution results remains to be audited.

## WCQ1:nonhyper-only — additional genus-3 class support

**Domain:** the published Type I/II `(2,2)` coverings over `F_(p^2)` of
ordinary elliptic models over `F_(p^6)`, with all distinct-root, twist and
isogeny-condition hypotheses enforced. The destination must actually be a
smooth non-hyperelliptic genus-3 curve.

**Finite conjecture:** at least one of the 46 exact p7 union-zero classes has
a Type I or Type II non-hyperelliptic covering representative.

**Unique property sought:** a class with no representative in either current
hyperelliptic branch, but with an independently verified non-hyperelliptic
covering representative. The known existence of non-hyperelliptic models
alone does not establish this class difference.

**First test:** use the frozen 46 signed trace classes as the evaluation
population. A positive requires the actual covering, smoothness,
non-hyperellipticity and the Weil-restriction isogeny condition. One positive
witness proves the finite existential statement. Its refutation requires a
complete declared Type I/II oracle over all 46 classes; a capped search
remains unresolved.

The Type I/II geometries are already classified by
[Momose–Chao](https://eprint.iacr.org/2009/236). The proposed contribution is
their exact support beyond the current two-branch class union.

## WCQ1:two-towers — primitive curves with two different base fields

**Domain:** ordinary `E/F_(p^6)` with j of exact absolute degree six and an
isogeny class that does not descend to any proper subfield. Seek a genus-2
cover over `F_(p^3)` for the degree-two tower, and a genus-3 cover over
`F_(p^2)` for the degree-three tower, both satisfying their actual
Weil-restriction isogeny condition.

**Finite conjecture:** such a jointly supported model exists at p=7.

**Distinctive property:** two different descent geometries for the same
primitive model and class. The existing *nonsplit quadratic branch* uses a
quadratic root factor inside a **degree-three** extension; it is different
from the degree-two tower in this conjecture.

**First test:** intersect the two published geometric conditions on a fixed
synthetic p7 family population. Record j degree, isomorphic-model degree,
and minimal isogeny-class descent degree separately. Reject apparent
witnesses obtained by base extension or by an isogeny to a base extension.
Verify both actual coverings independently. A complete empty intersection
refutes the finite conjecture; bounded absence leaves it unresolved.

Quadratic-tower covering constructions are already part of
[Scholten-form work](https://eprint.iacr.org/2005/277). The simultaneous
primitive intersection is the proposed property; novelty remains unresolved.

## WCQ1:smaller-genus — primitive higher-degree special loci

**Domain:** relative prime degree n=5 or 7 over `F_q`, q odd, and the existing
compositum covers of genus `G_n=1+2^(n-2)*(n-2)`, hence 25 or 161. Require
an ordinary source with primitive j and no proper-subfield class descent.

**Conjecture:** for at least one of n=5,7 there is a primitive special locus
with an `F_q`-defined quotient curve D of genus `n<=g(D)<G_n`, whose Jacobian
maps onto the full n-dimensional Weil restriction of the source.

**Distinctive property:** reduced destination genus while retaining the
complete source abelian variety. A small quotient carrying only one
geometric elliptic factor does not meet this statement. The lower bound n
is required by dimension, as in the published isogeny-condition framework.

**First test:** analyze automorphisms and quotient factors in prescribed
synthetic n=5,7 subfamilies before any large-field sampling. Keep the base
field of every automorphism and quotient explicit. The elementary choice
`alpha^q=c*alpha` with c in the base field forces the displayed Legendre
parameter into the base field; reject that route as a primitive-j witness.
Non-normalizing symmetries and different cover presentations remain
separate possible mechanisms.

**Falsification:** complete negative results reject only the frozen
subfamily. Refuting the existential conjecture requires a theorem excluding
all primitive loci, not a finite search cap. Audit candidate loci against
[the general covering classification](https://eprint.iacr.org/2006/347) and
[Iijima–Momose–Chao](https://eprint.iacr.org/2013/487) before calling them new.

## WCQ1:primitive-j — degree rigidity of quadratic torus lifts

**Domain:** p>3, q=p^2, odd n>=3, and a nonidentity lambda of exact absolute
degree 4n in the torus of order `(q^n+1)/(q+1)` inside `F_(p^(4n))`. Put
`j=256*(lambda^2-lambda+1)^3/(lambda^2*(1-lambda)^2)`.

**Conjecture:** j has exact absolute degree 2n. This extends the recorded
degree-six result to primitive higher lifts and prevents a mistaken
interpretation of those lifts as proper-field models.

**First test:** certify the complete Frobenius period of lambda, then every
proper candidate period of j, on frozen synthetic n=5,7 controls. A single
lambda meeting every hypothesis with smaller j degree refutes the universal
statement. Parameters of smaller degree are outside this domain.

**Proof route:** the anharmonic S3 action describes the Legendre j fibre;
its possible cyclic Frobenius orbit lengths constrain the lambda/j degree
ratio. Prove the actual finite-field fibre and exceptional cases in Lean,
rather than assuming a fibre-degree conclusion. Whether the precise
statement is already implied explicitly by the general classifications is
an open novelty-audit obligation.

## WCQ1:quartic-profile — polarization-dependent cover properties

**Domain:** actual hyperelliptic genus-3 covering Jacobians from the two
branches, with rational maximal isotropic 2-torsion kernels and nonsingular
principally polarized targets. Record whether each target is the Jacobian
of a plane quartic, a hyperelliptic curve or a product.

**Finite conjecture:** among the complete prescribed p7 covering models,
there are two such coverings from different branch origins with the same
Frobenius polynomial and the same rational 2-torsion Frobenius cycle type,
but different target-type multisets over rational maximal isotropic kernels.

**Distinctive property:** identical coarse arithmetic data, different
polarization-dependent geometry. This would justify storing destination
geometry separately from trace and torsion labels.

**First test:** enumerate the mathematical kernel/type data only for the
fixed prescribed covering models; obtain independent target-type
certificates. One verified pair proves the finite statement. A complete
oracle with no such pair refutes it; incomplete target classifications
remain unknown.

[Smith](https://arxiv.org/abs/0806.2995) already studies hyperelliptic-to-
non-hyperelliptic isogenies. The proposed study concerns conditioned target
profiles of these particular covering populations, not a new isogeny method.
