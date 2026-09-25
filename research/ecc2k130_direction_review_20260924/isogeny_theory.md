# ECC2K-130 isogeny, endomorphism-order and factor-base direction review

Date: 2026-09-24. Scope: public Certicom challenge and generated controls.
The exact-target probe and independent replay are archived alongside this
note. The existing fixed-presentation degree-search note and Rust module
carry corrected interpretation; the numerical sweep artifacts are unchanged.

## Main result and decision

**OBSERVATION, EXACT TARGET:** The 263-isogeny neighborhood can be constructed
directly using rational torsion on the quadratic twist. The accompanying
preflight constructs every one of the 264 degree-263 kernel lines, finds
262 distinct descending codomain j-invariants and two horizontal loops,
and verifies explicit x-coordinate transport on representatives of all
four Frobenius orbits. Both deterministic runs passed, in 3.27 and 3.78
seconds on the recorded host. These are feasibility timings, not comparative
performance evidence. This removes an implementation obstacle in the
existing research notes; it does not establish easier point decomposition.

The recommended research direction is a **joint curve, factor-base,
target and solver experiment** on actual transported instances, scored by
cost per verified independent relation and then full ECDLP cost. The
existing observations about unsuccessful toy curve selection remain useful.
They do not imply that every curve/basis/target combination has equal
cost, or that a selector must find a curve that is easy for every target.

## Exact parameters, orders and what descent changes

**THEOREM, target-specific derivation.** The challenge is

\[
 E_0:y^2+xy=x^3+1,\qquad q=2^{131},\qquad \#E_0(\mathbf F_q)=4r,
\]

with

\[
 r=680564733841876926932320129493409985129.
\]

The field polynomial and public subgroup order agree with
[Certicom's challenge parameters](https://www.certicom.com/en/curves-list).
The base-field Frobenius satisfies
\(\tau^2+\tau+2=0\). Hence its field is
\(K=\mathbf Q(\sqrt{-7})\), and
\(\mathbf Z[\tau]=\mathcal O_K\), the maximal order, because the
polynomial has fundamental discriminant -7. Since this order is already
contained in End(E_0), End(E_0)=O_K.

Applying the exact recurrence
\((A,B)\mapsto(-2B,A-B)\), initially (1,0), gives

\[
 \pi=\tau^{131}=A+B\tau,
\quad A=8123678690673902378,
\quad B=38531015900842053623
      =263\cdot146505763881528721.
\]

The trace is \(2A-B=-22283658519494248867\). Thus the conductor of
Z[pi] is B, while **the conductor of End(E_0) is 1**. These are different
orders and must never be assigned the same discriminant.

| Object | Conductor in O_K | Discriminant |
|---|---:|---:|
| End(E_0)=O_K | 1 | -7 |
| End of a degree-263 descending codomain | 263 | -484183 |
| Z[pi] | 38531015900842053623 | -7 B² |

Isogenies preserve the rational endomorphism algebra K. Descent replaces
an order by a suborder: a prime-degree descending step multiplies the
conductor by that prime and the order discriminant by its square. It does
not change the CM field or supply more endomorphisms. The order-theoretic
rules are standard in [Sutherland's Isogeny Volcanoes, §§2.6–2.10](https://arxiv.org/pdf/1208.5370)
and [Kohel's thesis bibliography and original thesis link](https://www.i2m.univ-amu.fr/perso/david.kohel/pub/index.html).

For rational descent from this maximal order, a prime must divide B.
The smallest is 263. The ramified 7-isogeny and the horizontal split-prime
isogenies return to j=1 because h(-7)=1. For 263, the quadratic symbol
(-7/263)=+1, giving two horizontal loops and 262 descending edges. The
263-volcano has depth one because v_263(B)=1. This conductor arithmetic
was already implemented in `src/cryptanalysis/isogeny_class_search.rs`
and `isogeny_degree_search/class_number.rs`; it is confirmed, not novel.

## The twist construction that opens the neighborhood

**THEOREM, with computational checks below.** A mod 263 is -1 and B mod
263 is zero, so pi acts as -1 on E_0[263]. Therefore all its torsion
abscissae are F_q-rational, although its nonzero 263-torsion points are
not E_0(F_q)-rational. On the quadratic twist

\[
 E_1:y^2+xy=x^3+x^2+1
\]

the q-Frobenius action is negated; all 263² torsion points are rational.
The exact twist order and projection cofactor are

\[
 \#E_1(\mathbf F_q)=2722258935367507707684713200934651442782,
\]
\[
 v_{263}(\#E_1)=2,\qquad
 H=\#E_1/263^2=39356632817700237211535705315020478.
\]

Since H is coprime to 263, multiplying a uniform random twist point by H
gives a uniform point in E_1[263]. The probe samples nonzero abscissae,
so its actual sampling distribution is slightly different; independence
is checked explicitly. Two independent points generate the
entire two-dimensional torsion module. The 264 cyclic kernels are the
lines \(\langle U+sV\rangle\) for s in F_263, and \(\langle V\rangle\).
Independence needs only a lookup over the 263 multiples of U.

For each line, enumerate the 131 distinct x-coordinates of ± pairs.
Writing this half-kernel set S and \(t=\sum_{u\in S}u\), the normalized
characteristic-two Vélu quotient has

\[
 b'=1+t+t^2,\qquad
 X(x)=x+A(x)+A(x)^2,\qquad
 A(x)=\sum_{u\in S}\frac{u}{x+u}.
\]

These formulas already have a derivation and small-degree modular-
polynomial controls in `src/cryptanalysis/binary_velu.rs`. Pairing ±Q
in Vélu's x sum gives
\(x(P+Q)+x(P-Q)=ux/(x+u)^2\), which equals
\(u/(x+u)+(u/(x+u))^2\). The codomain coefficients before normalization
are a4'=t, a6'=1+t; the change y→y+t removes a4' and adds t² to a6'.
The same abscissa kernel gives the isogeny on the untwisted target, with
a2=0 on the codomain. This does not confuse the quadratic twist's DLP
with the target's DLP: the twist is used only to construct the kernel.

The preflight uses the existing Python `FastGF2m` and `Koblitz` arithmetic,
which already supports a2=1, plus an extended-Euclidean inverse checked
against the existing exponentiation inverse. It writes dependency and
source hashes and exact coordinates for independent replay.

| Verified item | Seed 20260924 | Seed 20260925 |
|---|---:|---:|
| Independent rational 263-torsion basis | yes | yes |
| Kernel lines | 264 | 264 |
| Distinct torsion abscissae across lines | 34584 | 34584 |
| Kernel-line Frobenius orbit lengths | 1,1,131,131 | 1,1,131,131 |
| Distinct descending j-invariants | 262 | 262 |
| Horizontal loops | 2 | 2 |
| Full abscissa and codomain Frobenius covariance | yes | yes |
| Representative x-map scalar checks | 1,2,3,17,263 | 1,2,3,17,263 |
| Public P,Q images lift into order-r subgroup | yes | yes |
| Generic oriented transport interface | not implemented | not implemented |
| PDP or full-ECDLP improvement | not measured | not measured |

The explicit x maps are checked on one representative of each orbit,
including both horizontal loops. Codomain parameters and abscissa
covariance are checked for all 264 kernels. A successful x check means
the relation is preserved up to sign; an oriented relation collector
must implement the y map or certify every sign choice by group addition.
The four stored representative kernels provide reproducible controls for
a generic transport implementation.

**Independent red-team check, reported after this run:** a separate agent
reconstructed the four saved kernels from the first seed using independent
bit-polynomial arithmetic and evaluated the full Vélu X,Y sums over the
formal quadratic extension s²+s=1. For P, Q, P+Q and 2P it checked the
normalized codomain equation, group addition and agreement with the saved
X coordinates. This supports the map derivation independently of the
imported arithmetic; the executable preflight itself remains x-only.
The red-team artifact records its exact scope and reproducibility limits.

**Correction to the old construction estimate:** the target's division
polynomial does not have 264 *irreducible* degree-131 factors over F_q.
It splits into 34584 linear factors, grouped into 264 kernel polynomials
of degree 131. Computing that entire division polynomial is unnecessary
for this target. The twist construction builds those groups directly.

## Cheap Frobenius: lost on one leaf, available between conjugate leaves

**THEOREM, distinction relevant to accounting.** On a genuine descending
codomain E', tau is not in End(E')=Z+263 O_K. Coordinate squaring sends
E' to its conjugate E'^(2), not back to E'. The convenient endomorphism
of the original Koblitz model is therefore absent as a degree-two
endomorphism on a fixed leaf.

Nevertheless, an isogeny phi transports the action on the order-r
subgroup by

\[
 T'= [263^{-1}\bmod r]\,\phi\tau\widehat\phi.
\]

This is a perfectly valid action on the finite subgroup. It is not a
free coordinate squaring: charge the dual map, forward map and scalar
operation, or prove and measure a cheaper implementation. Globally,
phi tau phihat corresponds to 263 tau and has degree 2·263². Statements
that the group action is permanently unavailable on a leaf overreach;
statements that it is equally cheap overreach in the other direction.

The 262 descending j-invariants form two Frobenius orbits of size 131.
There can be no shorter nontrivial orbit because 131 is prime and none
of those j is in F_2. Thus structural screening needs only two
representatives **if the factor bases, targets and field bases are
conjugated along with the curve**. Holding a polynomial subspace fixed
while conjugating b tests basis alignment instead and need not give
identical solver behavior. The preflight separately checks kernel-line
orbits and actual quotient j values; it does not infer the latter from
group action alone.

## What the current code and notes settle, and where claims overreach

1. `binary_isogeny.rs` provides the older j-only walk with degrees 2 and
   3. `binary_velu.rs` and `examples/koblitz_isogeny_transport.rs` add
   real x-coordinate maps on small fields. The class-cost module's
   division-polynomial/kernel discussion predates these explicit maps.
   Its same-secret-on-each-curve benchmarks remain **untransported**
   independent instances, as its own scope correctly states.
2. `isogeny_class_search.rs` already distinguishes the homogeneous
   top-part statistic from affine refutation degree. Preserve its toy
   holdout negatives. The headline of
   `RESEARCH_ISOGENY_CLASS_SEARCH.md:13–47` is stronger than its evidence:
   lack of a curve that is easiest on every target does not rule out
   lower average cost, higher relation yield or charged target routing.
3. The corrected `RESEARCH_ISOGENY_DEGREE_SEARCH.md` gives a qualified
   fixed-input statement: fixed field representation, target and
   subspace, only b's constant vector changes. Under the stated
   definition, a regularity statistic computed from identical top
   homogeneous generators is identical. This does **not** identify
   actual solving degree, affine refutation degree, Gröbner runtime,
   SAT runtime, solution count, liftability, or relation rank. Higher
   summation-polynomial encodings require their own qualification;
   their b dependence is not merely constant.
4. The previous degree-search note correctly noted that S3 ignores a2,
   but incorrectly spoke of twists as halves of one isogeny class. Here the quadratic twists have opposite nonzero
   traces and are different F_q-isogeny classes. Their rational lifting
   constraints differ although their x-only S3 equation is identical.
5. The previous degree-search note contained two concrete
   mathematical errors: ordinary Verschiebung is separable, and the
   rational 7-isogeny is a self-loop at this class-number-one crater.
   The first distinct j over F_q is reached by degree 263.
6. The field invariant-subspace calculation {0,1,130,131} is correct:
   ord_131(2)=130. It excludes useful small Frobenius-stable **linear**
   F_2-subspaces and the low-magic-number standard GHS construction.
   It does not exclude nonlinear Frobenius-stable sets, non-invariant
   subspaces with orbit unions, arbitrary coefficient linearized
   polynomials, or every rational-map factor base. Broad statements
   about *all* quasi-subfields require the precise definition and
   coefficient restrictions used in the local census.

Suggested replacement for the overbroad class-search conclusion:

> For the fixed encodings examined, the highest homogeneous components
> are independent of the curve parameter under the stated target and
> factor-base choices. Toy curve selection did not retain a useful
> holdout advantage in the reported experiments. These results do not
> settle joint curve/factor-base selection or target-conditioned
> routing scored by verified relation yield and full pipeline cost.

Suggested replacement for the small-degree paragraph:

> In characteristic two, Frobenius is inseparable while Verschiebung
> is separable for an ordinary curve. At j=1 both lead back to the same
> j. Degree 7 supplies a ramified rational horizontal self-isogeny;
> since h(-7)=1 it does not reach another j. The first rational
> descending degree is 263, with 262 distinct descending neighbors
> and two horizontal loops.

### A small formal counterexample to confusing top parts with refutation degree

Work in the Boolean quotient F_2[x,y,z]/(x²-x,y²-y,z²-z). Let

\[
 f_1=xy+1,\quad f_2=xz+1,\quad f_3=yz,\quad
 f_4=xy+xz+yz+c.
\]

The positive-degree parts are identical for c=0 and c=1, and both
systems are unsatisfiable. At c=1, the sum of all four equations is 1,
so the affine Macaulay refutation degree is 2. At c=0, no constant
linear combination of the quadratic equations is 1: the only
dependency of the three quadratic monomials is f1+f2+f3+f4=0. A
degree-three refutation is

\[
 (y+1)f_1+yf_2+xf_3=1.
\]

Consequently D*=3 for c=0. This is a restricted Boolean-Macaulay
counterexample to the inference about affine D*, not a counterexample
to the tautological invariance of a statistic defined solely on the
top homogeneous generators. It does not claim an ECC curve advantage.

## Three candidate research programs

### A. Conservative extension: select for useful relation rate, with structure-aware routing

**HYPOTHESIS, UNTESTED:** Joint selection of curve, subspace/basis,
target features and solver reduces the cold cost per independent
verified relation on held-out instances. The mechanism is variation in
solution density and affine structure; equal top forms do not preclude it.

Start with the existing frozen two-/three-point toy corpus and the
current normal-basis support implementation. Add cheap, pre-solve
features: actual factor-base point count, Boolean rank of linear and
quadratic coefficient maps, bilinear ranks after fixing one block,
predicted decomposition multiplicity, trace-lift compatibility, and
previously measured orbit duplication. Use Gray/FES for tiny residual
dimension, SAT with native XOR handling for sparse parity structure,
and the existing Gröbner/crossbred options as a portfolio. Every
selector feature computation is charged.

Compare against an unconditioned best fixed solver and a shuffled-label
selector, with frozen targets and seeds. Separate true SAT, true UNSAT,
extraneous x solutions and rank-increasing group relations. A low
UNSAT-refutation degree alone is not an attack success. Require at
least a 10% paired reduction in total work per independent relation on
holdouts before expanding the campaign; require complete ECDLP
comparison before claiming a method improvement. Falsify the selector
if the gain vanishes after accounting for rejection rate, feature
cost, or selection on held-out seeds.

### B. Representation change: degree-263 transported bases versus native leaf bases

**HYPOTHESIS, UNTESTED:** A leaf and a native simple factor base can
provide a cheaper representation of the same subgroup DLP after map,
lift and orbit costs are included.

Use the preflight's two descending orbit representatives. Build an
oriented phi on at least one representative, or use x-only transport
with group-verified sign resolution. On four tractable toy sizes with
real isogenies, compare (i) original curve/native base, (ii) leaf with
the exactly transported original base, (iii) leaf/native base of
matched cardinality, and (iv) original curve with pullback of (iii).
The transported-base arm is a positive equivalence control. Random
basis changes and Frobenius-conjugated copies are negative controls
for accidental coordinate gains.

Preserve phi(P), phi(Q) for the same planted input. Count all setup,
kernel construction, transport, failed attempts, lifting, relation
verification, rank and final recovery. On the real target, gather only
bounded structural/solver-stage diagnostics until a viable pipeline
exists. A meaningful first gate is a persistent improvement in
charged relation cost with identical verified workloads; no arithmetic
above predicts that such an improvement exists. Falsify the specific
representation if its native base is equivalent to an already cheaper
pullback base or if losing cheap fixed-curve Frobenius dominates.

### C. Speculative direction: semilinear relation collection across conjugate leaves

**CONJECTURE, HIGH RISK, UNTESTED:** Keeping a tagged family
E',E'^(2),...,E'^(2^130) might retain useful cheap Frobenius relations
while using a non-Koblitz coordinate representation.

For phi_i=Frob^i phi tau^(-i) on the subgroup, coordinate squaring
relates logs against generators phi_i(P) by the known scalar lambda
of tau. This offers a way to share labels and preprocessing between
conjugates. It does not permit adding points on different curves for
free, and does not give an automatic factor-131 relation gain: fully
conjugated instances are isomorphic copies of the same problem.

Minimal test: build the tagged action on a toy descending orbit,
verify every group/log label, and compare cold setup and useful rank
against a single-leaf collector. Check whether any proposed cross-curve
relation reduces to transported versions of already counted relations.
Falsify a claimed advantage if canonicalization yields duplicate
systems/rows or if combining points requires map work that exceeds the
savings. The useful outcome even on failure is a precise statement of
which Frobenius accounting can be transported and which cannot.

## Literature overlap and limitations

Isogeny-assisted GHS was already established by Galbraith, Hess and
Smart, *Extending the GHS Weil Descent Attack*. Their construction
uses isogenous curves with different descent behavior over composite
extension fields; it is precedent for representation search, not
evidence of a low-genus path for degree 131. The author's
[Weil descent page](https://nigelsmart.github.io/weil_descent.html)
links the original works. The exact invariant-subspace calculation
above is stronger for the narrowly specified standard GHS magic
window on this target, and weaker than an exclusion of all descent.

Galbraith–Granger–Merz–Petit, [*On Index Calculus Algorithms for Subfield Curves*](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf),
already studies Frobenius-invariant factor bases and orbit reduction.
It also examines nonlinear/rational-function constructions from
Couveignes–Lercier. Its experiments show these rational-function
encodings can cost much more despite symmetry; the authors explicitly
leave structural solver exploitation open. Thus invariant linear
subspace obstructions cannot be promoted to exclusions of every
Frobenius-invariant factor base. Their favorable linear-subspace
experiments use parameters different from degree 131. No novelty is
claimed for these general mechanisms.

The matching rho reference for this public subgroup includes negation
and 131-fold Frobenius, giving the heuristic iteration count
sqrt(pi*r/(2*262)), approximately 2^60.81. This is an iteration-count
reference, not a measured hardware time. The primary target-specific
discussion is [Bailey et al., *Breaking ECC2K-130*](https://ecc-challenge.info/anon.pdf).
No claim of beating that reference is made here.

## Reproduction and immediate next action

```sh
python3 research/ecc2k130_direction_review_20260924/twist_torsion_preflight.py --out /tmp/twist_torsion_replay.json
```

Preserved artifacts: `twist_torsion_contract.md`,
`twist_torsion_preflight.py`, and `twist_torsion_results.json`.
The result is exact-target structural evidence, with 129-bit subgroup
membership checks and explicit degree-263 x maps. It is not toy
key recovery, a PDP speed measurement, or a challenge solution.

**Next concrete action:** use the saved independent representative replay
as a control while implementing one full oriented y map and compare
native versus transported factor bases on matched tiny ECDLP
instances. Keep the proposed selector work separate from the map
feasibility result so a cheap UNSAT result cannot be mistaken for
productive relation collection.
