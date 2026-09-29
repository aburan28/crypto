# Native descendant factor bases: an m≥3 pullback control

Status: **conditional protocol and counting preflight only**. No new relation,
solver, ECDLP, or ECC2K-130 timing measurement is reported here. Freeze a
separate input manifest, target stream, source hashes, and resource caps before
any measured successor run. Do not use this note to dispatch a leaf run while
the implicit full-point PDP exporter and its negative-result semantics remain
unadmitted.

## Question and existing boundary

The source ECC2K-130 curve has endomorphism discriminant −7. Its certified
degree-263 descending leaves have conductor 263 and discriminant −484183;
they do not inherit a cheap *geometric* self-Frobenius. The maps and normalized
dual are already certified by [#743](https://github.com/aburan28/crypto/pull/743)
and [#750](https://github.com/aburan28/crypto/pull/750). The ring change is a
real structural difference, but by itself says nothing about easier point
decomposition. The class search's leading-form and D* claims are already
examined in [the class-search note](../notes/ecc2k130/RESEARCH_ISOGENY_CLASS_SEARCH.md);
this protocol does not repeat that search.

The [#753 factor-base gate](../ecc2k130_factor_base_replication_20260925/PROTOCOL.md)
paired source, transported, descendant-native and source-pullback bases at
**m=2**. Native/pullback hits and rank matched exactly; the degree-7 toy native
arm lost its cold field-multiplication component to original in all four
holdout cells. Its two degree-263 challenge smokes had zero hits at a
precomputed expectation of at most 6.40×10⁻³⁶ per arm. The
[#751 three-summand gate](../ecc2k130_pdp_chain_holdout_20260925/RESULT.md)
paired source and transported bases, not descendant-native and pullback bases.
It found exact relation parity but a slower chained affine prescreen. Thus the
unmeasured question is **whether an m≥3 native leaf encoding can repay its
construction and action costs against the exact source pullback on the same
point set and targets**. A native-versus-unrelated-source comparison cannot
isolate this question because its support sets differ.

An m=3 degree-7 toy control would test the mechanism, **not** admit m=3 on
ECC2K-130. Reusing the only existing exact leaf smoke's `s=16` projected
points at n=131 gives at most `C(18,3)=816` distinct unordered three-sums.
With `R=680564733841876926932320129493409985129`, 32 uniform subgroup
targets have at most `32×816/R ≈ 3.84×10⁻³⁵` expected hits. Even a 1%
*necessary* target-support ceiling would need at least
`s=3,443,553,994,135` physical points in one unordered common base; this
number says nothing about whether such a base can be represented or solved.
No larger native leaf base or n=131 m≥3 implicit representation is admitted,
so transfer from a positive toy cell remains unset. A separate higher-arity,
actual-count and memory gate is mandatory.

## Exact paired invariant

Let φ:E→E′ be a certified degree-ℓ isogeny over the challenge field and
φ̂:E′→E its normalized dual, with φ̂φ=[ℓ] and φφ̂=[ℓ]. Work on the
prime-order subgroup of order R, where gcd(ℓ,R)=1. For every nonzero useful
native leaf point b′, define its source pullback

    b = [ℓ⁻¹ mod R] φ̂(b′).

Then φ(b)=b′. For a native set B′ and B=φ⁻¹(B′), and each source target T,
the ordered, repeated-point, signed and exceptional-branch **group-law**
witnesses of `T=Σ b_i` correspond bijectively to those of `φ(T)=Σ b′_i`.
Consequently exact supported-target membership and per-target witness counts
are equal. If both sides label columns by corresponding *physical points* and
emit corresponding rows, row vectors and rank modulo R are equal up to a
column permutation. Any orbit compression must use a separately verified
conjugate action and a matching column partition on both sides; native leaf
coordinate squaring is not a self-map and may not be counted as free.

This is an algebraic positive control, not a performance result. A leaf can
still offer a better coordinate equation, solver ordering or implicit
representation. Those gains have to beat map, dual, base, action and failed
query costs. For raw unprojected full-group points, replace `R` by the exact
full-group order in the inverse formula when ℓ is coprime to that order;
the prime-subgroup comparison below deliberately uses projected points.

## Admission before any measured pair

1. **Freeze one common instance.** Start with a fully solvable small rung
   with a certified descending map; the n=21 degree-7 code from #753 can be
   a semantic control but is ramified and cannot model the split n=131 yield.
   Freeze a second, independently generated base and disjoint train/holdout
   target streams. After the small-rung gate, use the two nonconjugate saved
   degree-263 lines `[1,0]` and `[1,4]` for bounded n=131 representation
   checks. No leaf line, base, field basis, arity, solver or target may be
   selected from holdout outcomes. Publish the exact parent commit, source
   hashes, map certificate, curve/field model, useful point count and seed
   derivation before dispatch.
2. **Admit an implicit full-point producer first.** It must cover infinity,
   `x=0`, inverse and finite-fibre branches; lift every positive model to an
   exact signed group sum; and distinguish checked UNSAT from timeout, memory
   cap, unsupported input and unknown. The same producer semantics and solver
   budget apply to native and pullback. If one arm cannot represent the base
   within the frozen memory cap, record representation failure; do not infer
   a point-decomposition or end-to-end speed difference. Do not enlarge a
   materialized-root table merely to make this comparison run.
3. **Pass the necessary support/rank screen.** Use
   [`admission.py`](admission.py) on the **physical** point counts, arity,
   exact uniform target support size, target cap and required independent
   rank. For one common unordered base of size `s`, the m-sum support has at
   most `C(s+m−1,m)` elements. For distinct ordered slots, use `∏s_i` instead;
   the rotated-slot corpus is of this second type. With `N` target probes
   uniform marginally on a group or coset of size `U`, and at most `b` emitted
   rows per target, `Pr(rank≥K) ≤ min(1,N b min(1,C/U)/K)`. This is Markov's
   inequality applied to the number of supported target probes. Correlated
   probes are allowed; deliberately choosing supported targets is not. Set
   `N` to **all** probed torsion translates and choose `U` for their actual
   marginal distribution. One fixed coset has `U=R`; a genuinely uniform
   four-coset full-group target has `U=4R`. Do not silently multiply a
   one-coset result by four. If several independent rows can be emitted per
   target, declare and enforce their maximum `b`; otherwise the rank bound
   is inapplicable. If `N b < K`, full rank is deterministically impossible.
   A small upper bound rejects an underpowered proposed run, not the
   underlying isogeny hypothesis. A large bound admits only a test, not a
   hit-rate prediction.
4. **Freeze paired arms.** Use (A) the current original source construction,
   (B) its exact transported leaf image, (C) a native leaf base selected by
   the same predeclared construction budget, and (D) the exact source pullback
   of C. Match useful **physical-point** counts and, where a symmetry is
   claimed, independently verify the orbit partition. Use the same natural
   source subgroup Q stream and map each target exactly for B/C. Charge
   rejected target draws and all failed or negative PDP probes. Never compare
   only planted positives. Verify A/B and C/D support, witness and rank parity
   with a separate point-law replay; keep any mismatch as a correctness
   failure, not a timing datum.

## Cost ledger and falsification

For each arm publish one cold, common-unit ledger from the public P,Q through
the first verified full rank and at least one independently verified held-out
Q logarithm, together with the matched same-Q automorphism-aware batched rho
reference. Preserve separate native field multiplication, squaring, inversion,
group-addition, modular-R, CPU and RSS counters until a measured conversion
puts all arithmetic in one unit; do not add overlapping group and field
counters. The cold ledger includes torsion discovery; kernel and dual/map
construction; field/model conversions; every base candidate, cofactor lift
and rejected point; forward/pullback maps; actual or certified tagged-conjugate
leaf action and orbit labeling; target generation and all torsion translates;
every PDP positive, checked negative and unknown; model lifting; dependent
rows; sparse linear solve; target descent; and scalar verification. Report
shared setup once for a batch and arm-specific setup separately. Also publish
the post-setup paired solver cost for C/D, since D must build C before its
pullback and a cold comparison alone can hide a coordinate effect.

A C/D support or rank-parity mismatch is a **correctness failure** and
quarantines the cell until explained; it is not evidence for a mathematical
leaf advantage. Once both arms pass independent replay, the m≥3
**representation hypothesis fails on the frozen cells** if C cannot attain
the same verified full rank/logs within the cap or does not lower the fully
charged common-unit cost against both D and A on both disjoint holdouts for
both base seeds. Passing a solver-stage C/D comparison alone is only an
engineering lead.
Passing the toy degree-7 cells does not imply a degree-263 or n=131 gain;
passing a bounded n=131 PDP smoke without full rank does not imply an ECDLP
cost or a crossover. Preserve failed cells and a null speedup field.

## Theory sources

The dual identities and isogeny homomorphism are standard in
[Kohel's thesis, Chapter 2](https://www.i2m.univ-amu.fr/perso/david.kohel/pub/thesis.pdf).
The descending-conductor rule and orders in an ordinary isogeny class are
given in [Sutherland, *Isogeny volcanoes*, §§2–3](https://msp.org/obs/2013/1-1/obs-v1-n1-p25-s.pdf).
These sources justify the algebraic controls; the repository's saved
certificates establish the particular maps and curves named above.
