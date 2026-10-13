# Isogeny and torsion fibers: rational-kernel correction and first control

**Date:** 2026-10-12. **Status:** structural correction and bounded native
checks. **Classification:** accounting correction to the earlier fiber proof;
no algorithmic performance result. The protocols are [PROTOCOL.md](PROTOCOL.md)
and [MAP_PROTOCOL.md](MAP_PROTOCOL.md). The [run report](reports/run-20261012T1947Z/README.md)
records commands, outputs and validation. The [diagram](fiber.svg),
[census graph](census.svg), and [PDF](report.pdf) accompany this result.
The [original local validation receipt](validation/README.md) records the
release-library gap as it stood for this first control; the later
[paired-round validation](paired_round/validation/README.md) records a
complete controlled full-suite pass on the same source tree.

**Follow-on, 2026-10-12:** The [paired same-target round](paired_round/REPORT.md)
now executes explicit transport, source and destination relation
collection, linear recovery, target descent, verification, and matched
negation rho on 16 public targets. It has 48 verified arm records and a
separate native arithmetic checker. Its group-operation comparison is a
diagnostic; both destination models are also base-field isomorphic to
their sources, so an isomorphism-only relation arm is the next mechanism
control. Full native-work pricing, six-class rho, scaling and
independent-host replay remain open. The sections below document this
earlier structural control and retain its original run scope.

## Question and exact setting

Can fibers of isogenies or torsion maps make an ECDLP cheaper? A necessary
first distinction is between the geometric kernel `M = ker φ`, of size
`deg φ` for a separable isogeny, and its rational part
`K = M ∩ E(F_q) = M^π`. The preceding [fiber report](../../RESEARCH_ISOGENY_FIBERS_ECDLP.md)
correctly uses `K` to describe rational fibers, but its executive summary
and proofs of Props. 2.7 and 4.1 said that `deg φ | h`, where
`#E(F_q)=h r`, `r` is the target prime and `r ∤ deg φ`. That sentence is
false for kernels with non-rational points.

The subgroup and field in the explicit control are:

| object | exact identity and status |
|:--|:--|
| source | `E/F₁₃: y²=x³+2`, `#E(F₁₃)=19`, ordinary, prime-order rational group |
| geometric kernel | `M={O,(0,√2),(0,−√2)}` over the algebraic closure; order 3; Frobenius exchanges the nonzero points |
| rational kernel | `K={O}`, order 1, because 2 is nonsquare in `F₁₃` |
| quotient | `E'/F₁₃: y²=x³+11`, `#E'(F₁₃)=19`; the degree-3 formula is in [MAP_PROTOCOL.md](MAP_PROTOCOL.md) |
| target subgroup | all of `E(F₁₃)`, `r=19`, `h=1`; all of `E'(F₁₃)` is its transported image |

This is an exact field/model/subgroup identity for the public control.
No ICV1 or EC1 registry entry was added, and these model parameters must
not be conflated with a registered curve. The first witness in the frozen
census is `E/F₇: y²=x³+3` with order 13; the `F₁₃` pair was chosen for the
separate explicit-map check because its source and target equations differ.

## Correction and proof

**Corrected statement.** Let `φ:E₁→E₂` be a separable `F_q`-isogeny of
degree `d`. Let `#E₁(F_q)=h r`, with `r` prime and `r ∤ d`. Put
`K = ker φ ∩ E₁(F_q)` and `k=|K|`. Then `k | h`. Consequently, if
`ΣP_i−R ∈ K`, then `[h]ΣP_i=[h]R`. If cofactor projection is intended to
identify the nonzero order-`r` group, it additionally needs the usual
`r ∤ h` assumption.

**Proof.** `K` is a subgroup of `E₁(F_q)`, so Lagrange gives `k | hr`.
It is also a subgroup of the geometric kernel of order `d`, so `k | d`.
Since `r` is prime and `r ∤ d`, `gcd(k,r)=1`; Euclid's lemma gives
`k | h`. Thus `[h]` kills every difference in `K`. This repairs the
relation-projection step without asserting `d | h`.

The explicit `F₁₃` pair has `d=3`, `r=19`, `h=1`, and `k=1`; hence
`d ∤ h` while `k | h`. Over `F₁₃`, the isogeny restricts to a bijection
on rational points. Geometric fiber size three therefore supplies no
threefold rational choice. The bounded map check confirms its formula
maps all 19 source points bijectively and preserves all 361 ordered sums.
The quotient and homomorphism laws for arbitrary fields use the standard
isogeny construction described by [Shumow](https://arxiv.org/abs/0910.5370);
the finite check alone does not prove them.

## What ran

The frozen [census](PROTOCOL.md) covers every prime `p≤199` with
`p≡1 mod 3` and every nonsquare `b` in `E_b:y²=x³+b`. It scanned 1,014
curves. Direct point enumeration and independent quadratic-character
counts agreed for every row, and every trace passed the Hasse check.
There were 349 prime-order rows, recorded in [census.tsv](census.tsv);
[summary.tsv](summary.tsv) gives all 21 base primes, including three with
zero prime-order rows. The count is descriptive of this deliberately
selected family and finite range, not a density estimate.

The [map check](map_check.txt) on the named `F₁₃` pair passed 19 image
membership checks, 19 distinct-image checks, and 361 ordered-pair
homomorphism checks. Both sides had exactly 19 rational points. There was
no target-logarithm computation or performance measurement in either run.

| scope | source group order | `deg φ` | `|K|` | cofactor `h` | observation |
|:--|--:|--:|--:|--:|:--|
| first witness, `F₇`, `b=3` | 13 | 3 | 1 | 1 | `3 ∤ 1`, `1 | 1` |
| explicit-map control, `F₁₃`, `b=2` | 19 | 3 | 1 | 1 | rational map bijection and 361/361 pair checks |
| frozen census | 349 prime-order rows | 3 | 1 | 1 | two exact counts agree for all 1,014 curves |

## Proof map

| mathematical statement | Lean declaration or status | gap |
|:--|:--|:--|
| `k | hr`, `k | d`, `gcd(r,d)=1` imply `k | h` | `FiberCofactor.rationalKernelOrderDvdCofactor`, checked with `lake build` | Lagrange and kernel-order facts are explicit hypotheses; elliptic-curve objects are not encoded |
| `3 ∤ 1` | `FiberCofactor.degreeThreeNotDvdCofactorOne`, checked | none for the arithmetic statement |
| `E:y²=x³+2 → E':y²=x³+11` is a degree-3 isogeny over `F₁₃` | written standard formula and exhaustive rational-point checks | no full Lean formalization of Vélu's construction or geometric kernel |

The Lean project pins `leanprover/lean4:v4.34.1`, imports only `Init`, and
contains no `sorry`, custom claim axiom or native evaluation in the proofs.
`#print axioms` reported only `propext` and `Quot.sound` for the divisibility
lemma, and `propext` for `3 ∤ 1`. `lean4checker --fresh` was unavailable.
This is partial formalization of the corrected geometric statement, with
the exact missing bridge identified above. The retained [build output](lean/build.log)
records the theorem checks.

## Mechanisms and next cost test

The distinction narrows what a non-rational kernel can contribute over
`F_q`: its rational fibers are singletons whenever `K={O}`. It does not
rule out a cheaper coordinate-based relation oracle on an isogenous
representative. That would require a measured destination-solver advantage,
including the map cost and all phases of the solve. The older [fiber
rounds](../isogeny_fibers_20261010/round4/README.md) measured rational
kernel saturation, subgroup checks, vertical transport, and related
controls; those records remain the relevant baseline, with their stated
curve families and run limits.

| candidate route | possible source of cost change | first missing obligation |
|:--|:--|:--|
| small rational torsion quotient | factor-base coordinates and symmetrized relation systems | complete cost on the named workload; compare with the existing rounds |
| non-rational kernel and extension-field lifts | a different coordinate oracle after lift/trace, proposed only | explicit extension map, lift success rate, extension arithmetic, and independent relation checks |
| non-rational kernel with base-field isogeny | destination model or solver behavior, proposed only | same-target complete solve and matched rho, including transport |
| order-`r` torsion fiber | descent or pairing into another group | explicit connecting map, embedding field and destination DLP cost |

No row in this table has a measured end-to-end advantage in this study.

The `[r]`-torsion fiber is a different case: its kernel meets the hard
subgroup. A descent through it needs an explicit connecting map and a
destination logarithm calculation. Voloch's [descent paper](https://web.ma.utexas.edu/users/voloch/Preprints/disclog3.pdf)
relates the prime-to-characteristic case to multiplicative-group descent
and the characteristic-primary case to Artin–Schreier–Witt descent.
[Faugère et al.](https://www.iacr.org/archive/eurocrypt2014/84410158/84410158.pdf)
use small rational torsion to change summation-polynomial systems; their
result does not turn this non-rational kernel into rational multiplicity.
[Galbraith, Hess and Smart](https://research-information.bris.ac.uk/en/publications/extending-the-ghs-weil-descent-attack/)
show that isogeny transfer to a different coordinate solver can matter over
specified composite-degree binary fields. These are distinct mechanisms.

The next discriminating performance test is a preregistered native,
same-target pair: transport a public prime-order challenge through an
explicit non-rational-kernel isogeny; compare source and destination
relation generation and complete solves under equal inputs and caps;
include failed attempts, construction, map evaluation, linear algebra,
target recovery and verification; run matched rho and an independent
replay. A stage-only difference remains a stage diagnostic. This requires
a supported native prime-field complete-solve path and an appropriate
size where both arms and rho can be audited. The unresolved obligation is
tracked in [the problem record](../../problems/nonrational-isogeny-fiber-dlp-cost-gap.md).

No new DLP cost or isogeny-class edge was certified for the canonical
scoreboard or curve registry. `docs/index-calculus-scoreboard.html`,
`docs/ic-current-state.html`, and `docs/curves/` were checked and left
unchanged for that reason.

## 2026-10-13 follow-on: isomorphism-only and six-class controls

The earlier mechanism table and next-step paragraph above describe the
state before this follow-on. The [new result report](paired_round/iso_control_20261013/REPORT.md),
[run PDFs](paired_round/iso_control_20261013/reports/run-20261013T1634Z/README.md),
and [raw v3 records](paired_round/iso_control_20261013/raw-v3.tsv) now provide
complete relation solves after the public scaling isomorphism, alongside
fresh degree-three, source, negation-rho, and six-class-rho arms on the same
16 frozen source targets. An independent native checker verified all 80
answers and the orbit/multiplier action; v1/v2 six-class cap failures are
retained separately.

At order 211, the degree-three and isomorphism destination arms used 359
and 358 relation trials and 4,871 and 5,046 counted group operations over
the same eight targets. The isomorphism/degree-three ratio is 1.035927 in
that diagnostic unit. At order 19, the respective operation totals were
563 and 691. The fixed six-class rho control verified every target in v3,
with three charged restarts at order 19; it is not a calibrated best-rho
claim. Both destination curves are base-field isomorphic to their sources,
so these data still cannot isolate a degree-three-specific effect. The
next mechanism test needs a base-field non-isomorphic pair; an admitted
performance comparison additionally needs growing sizes, complete native
pricing, and independent-host replay. The [problem record](../../problems/nonrational-isogeny-fiber-dlp-cost-gap.md)
tracks the remaining work. The canonical scoreboard graphs were checked
and left unchanged because this is a bounded counted-operation diagnostic,
not a new complete-cost or scaling finding.
