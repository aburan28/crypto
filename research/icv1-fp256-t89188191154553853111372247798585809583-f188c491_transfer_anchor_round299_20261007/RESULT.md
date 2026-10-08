# P-256 transfer-anchor and structured multiplicity, round 299: result

## Verdict

No route reaches rho parity, but two ambiguities are now resolved.

First, a nonzero group-homomorphic lift of the prime-order P-256 subgroup is
injective and preserves the discrete-log scalar exactly.  A cover, isogeny or
Jacobian representation can still change decomposition cost, but transport
alone cannot quotient the 164 factor-base logarithm classes.  A zero
restriction kills the subgroup and cannot support recovery; every
endomorphism of the nonzero cyclic image is again scalar.

Second, the known-target multi-row threshold is 164 distinct normalized rows,
not 165.  If a fixed nonzero target logarithm is known, the first decomposition
is already an inhomogeneous equation and 164 coefficient rows can have rank
164.  Their homogeneous differences have rank at most 163 because the true
factor-base log vector is in the kernel.  If the target logarithm is unknown,
the augmented system retains an overall scalar ambiguity and does not recover
an absolute DLP.  This is an additive correction to Round 298, not a speedup:
the equivalent raw folded threshold is still 328 preimages and the P-256 mean
distinct occupancy is only 1.005109899.

The complete toy screen cleanly separates generic scalar structure from
geometric selection.  Affine-progression and pair-sum bases, both selected
using known scalar labels, produce sparse rank-6 buckets.  The two geometry-
defined families evaluate 1,848 complete six-column bases without using those
labels and produce zero rank-6 buckets.  No P-256 relation was attempted and
promotion remains negative.

## Requirement-to-evidence status

| requested direction / gate | status | Round-299 evidence or gap |
|:--|:--:|:--|
| pursue fundamentally different transport | verified within stated theorem scope | zero/nonzero homomorphic restrictions and image endomorphisms classified exactly |
| pursue many independent rows per event | verified on bounded prime-subgroup controls | six registered `B=6,m=3` families; every domain complete |
| screen different factor-base structures | verified within registered families | 71,992 globally unique bases; scalar-defined and geometry-defined controls separated |
| exact relation replay | verified | 39,133,744 raw toy rows; zero FP/FN |
| field/group operation, RAM and timing accounting | verified for the bounded experiment | exact scalar, modular-rank and transition costs plus isolated process telemetry |
| cover/Jacobian implementation | not implemented | no explicit construction or formulas were supplied; remains an open obligation |
| P-256 geometry-defined rank-164 event | not established | toy geometry trigger is 0/1,848; no extrapolation is made |
| structured residual degree `<=5` | not established | no P-256 candidate reached the solver gate; Round 296 remains degree 6 |
| cost per relation, full collection and storage gates | not established | no promoted P-256 oracle or yield model |
| actual unplanted P-256 relation | correctly not attempted | prerequisite projection gates fail |

Thus the bounded Round-299 protocol is complete, while the broader end-to-end
parity goal remains open and unpromoted.

## Prime-order anchor certificate

Let `H=<G>` have prime order `n` and let `phi:H->A` be the group homomorphism
induced on the relevant subgroup by a proposed correspondence.

1. `ker(phi) intersect H` is a subgroup of `H`.  Since `n` is prime, it is
   either `{0}` or all of `H`.
2. If `phi(G)=0`, then `phi([x]G)=0` for every `x`; recovery is impossible.
3. If `phi(G)!=0`, the kernel is trivial, `|phi(H)|=n`, and
   `phi([x]G)=[x]phi(G)`.  Therefore
   `log_{phi(G)}(phi([x]G))=x mod n`.
4. Since `phi(H)` is cyclic of prime order, every endomorphism of the image is
   multiplication by one scalar modulo `n`.

For same-target relation rows `r_i`, every difference satisfies

```text
(r_i-r_0) . z = 0 mod n,
```

where `z` is the actual vector of factor-base logarithms.  Consequently the
homogeneous rank ceiling is `B-1`.  For a known nonzero target `t`, the rows
instead satisfy `r_i.z=t` and their coefficient matrix may have rank `B`; the
first row supplies the affine offset and the remaining `B-1` dimensions come
from differences.

This certificate is complete only for group-homomorphic restriction to the
P-256 order-`n` subgroup.  It does not classify nonlinear set maps, prove that
no useful cover exists, or reject a destination decomposition algorithm.  The
open obligation is an executable correspondence with explicit fields,
dimensions, formulas, exceptional points, kernel, inverse recovery, relation
pullback and paired end-to-end cost.

## Corrected P-256 event boundary

| quantity | exact / modeled value |
|:--|--:|
| factor-base columns | 164 |
| distinct normalized rows for a known nonzero target | **164** |
| raw folded preimages before global-negation deduplication | 328 |
| homogeneous same-target rank ceiling | 163 |
| exact normalized domain | `58191887563875222213524008097554601544558227237900511180641613175619196026880` |
| canonical target space | `57896044605178124381348723474703786764998477612067880171211129530534256022185` |
| mean distinct occupancy | 1.0051098993154817 |
| random-map any-bucket Chernoff union diagnostic at 164 | `2^-713.831` |
| random-map Poisson expected-bucket diagnostic at 164 | `2^-720.277` |
| one free-oracle event / rho | 0.964088 |
| two free-oracle events / rho | 1.446132 |

The tail values are explicitly a random-map null, not a P-256 measurement or
proof of randomness.  They omit the cost of finding and materializing every
preimage of a selected bucket.

## Complete order-149 subgroup screen

The runner counts

```text
#E(F_1151) = 1192 = 8*149
H = [8](0,376) = (360,129)
[149]H = O
```

and enumerates all 149 subgroup points.  Negation leaves 74 nonidentity
columns and 75 canonical targets.  Each `(B,m)=(6,3)` base has 160 raw signed
rows, 80 distinct normalized rows and fixed distinct mean occupancy `80/75 =
1.0666667`.

| family | proposals | admissible | unique bases | bases with rank 6 | bases with homogeneous rank 5 | max distinct occupancy | best rank pair `(known,hom)` | selection uses logs? |
|:--|--:|--:|--:|--:|--:|--:|:--:|:--:|
| affine progression | 21,904 | 20,424 | 5,032 | **740** | 740 | **10** | 6,5 | yes |
| geometric progression | 21,904 | 21,312 | 2,664 | 0 | 0 | 6 | 5,4 | yes |
| pair-sum closure | 64,824 | 59,139 | 58,497 | **1,401** | 1,401 | **10** | 6,5 | yes |
| hash scalar | 4,096 | 4,096 | 4,096 | 0 | 0 | 8 | 5,4 | yes |
| low-coordinate geometry | 924 | 924 | 924 | **0** | 0 | 6 | 5,4 | no |
| hash-coordinate geometry | 924 | 924 | 924 | **0** | 0 | 6 | 5,4 | no |

After cross-family deduplication, 71,992 bases produce 11,518,720 raw rows and
5,759,360 distinct normalized rows.  All raw rows replay through an exact
elliptic-addition transition table with zero false positives, false negatives,
kernel failures or rank-ceiling violations.  The primary screen charges
34,556,160 scalar-label additions and the same number of logical group
additions.  Its six modular rank controls charge 57,100,246 additions or
subtractions, 2,095,739,750 multiplications and 28,370,872 inversions across
`F_149` and the formal `2^61-1` field.  Constructing the complete transition
table costs 22,052 exact toy group additions: 130,388 field additions or
subtractions, 457,320 multiplications and 21,756 inversions on `F_1151`.

The strongest affine base has folded labels `{1,3,5,7,9,11}`, only 14 nonempty
targets, occupancy 10 at target 1, known-target rank 6 and homogeneous rank 5.
That is real multiplicity, but it is manufactured by defining the base through
known discrete-log labels.  Reproducing that construction on P-256 is generic
log transport, not a geometric factor-base improvement.  The pair-sum result
has the same issue.  Neither coordinate-selected family triggers the
preregistered next-stage gate.

Ranks were computed both over the actual subgroup field `F_149` and modulo
`2^61-1` as a formal coefficient control.  No bucket in any registered family
showed a rank divergence.  Only the `F_149` ranks carry the toy anchor/kernel
interpretation.

## Exhaustive reference

The implementation reference exhausts every `C(74,4)=1,150,626` four-column
base at arity two.  It replays 27,615,024 raw rows and checks 13,807,512
distinct normalized rows, charging 55,230,048 scalar-label additions and the
same number of logical group additions.  Its rank controls charge 1,686,368
additions or subtractions, 3,964,764,098 multiplications and 56,354,404
inversions.

| reference property | result |
|:--|--:|
| maximum distinct occupancy | 3 |
| maximum known-target rank | 3 |
| maximum homogeneous rank | 2 |
| bases reaching rank 4 | 0 |
| formal/subgroup rank divergences | 0 |
| replay failures | 0 |
| false positives / false negatives | 0 / 0 |

Across the primary and exhaustive reference screens, 39,133,744 raw rows were
replayed with 89,786,208 charged scalar-label additions and 89,786,208 logical
group additions.  Modular rank accounting totals 58,786,614 additions or
subtractions, 6,060,503,848 multiplications and 84,725,276 inversions.  These
are toy/control-field operations, not P-256 field-equivalent projections.  The
enumeration digests are preserved in the JSON artifacts.

## Typed transfer assessment

| path | subgroup result | log result | measured advantage | status |
|:--|:--|:--|:--|:--|
| `phi|H=0` | entire subgroup in kernel | all logs erased | recovery impossible | refuted as useful transport |
| nonzero homomorphic restriction | injective, image order `n` | scalar exactly preserved | no destination solver supplied | transport theorem supported; advantage unknown |
| image endomorphism | scalar `[a] mod n` | known scalar orbit or zero map | generic / none | non-generic quotient refuted |
| nonlinear or extra destination relations | construction-specific | unknown | unmeasured | open obligation |

The full obligation table, typed graph, controls, exploration boundaries and
weakest open obligation are in `transfer-assessment.json`.  The transfer
skill's main workflow was available; its referenced methodology and JSON
template were unavailable, so that limitation is recorded in the artifact.

## Gates

| gate | status | evidence |
|:--|:--:|:--|
| dependency hash | pass | corrected Round 298 SHA-256 pinned |
| toy group and order-149 subgroup exact | pass | 1,192 ambient and 149 subgroup points enumerated |
| registered families complete within scope | pass | proposal/admission/unique counts frozen |
| all-base `B=4` reference complete | pass | 1,150,626 / 1,150,626 bases |
| exact replay and zero FP/FN | pass | 39,133,744 raw rows; zero errors |
| subgroup rank ceilings | pass | homogeneous `<=B-1`, augmented `<=B` everywhere |
| anchor eliminated by homomorphic transfer | **fail** | nonzero restriction preserves scalar; zero restriction kills recovery |
| geometry-defined sparse rank-6 candidate | **fail** | 0 / 1,848 bases |
| P-256 construction and scaling | **fail / unset** | no geometry-defined base or destination correspondence |
| structured residual degree at most 5 | **fail / unset** | no promoted P-256 system; Round 296 remains degree 6 |
| parity at or below rho | **fail** | no qualifying rank-164 P-256 event or quotient |
| usable relation below `2^103` | **fail / unset** | no P-256 oracle |
| collection below `2^120` | **fail / unset** | no P-256 yield/cost |
| storage below `2^50` | **fail / unset** | no promoted implementation |
| promoted | **no** | parity-critical gates fail |

## Resources and reproduction

The final canonical and independent runs were isolated, uncontended and pinned
to CPU 4 of the AMD EPYC 9V74 host.  The canonical run used 31.629179 s wall,
31.608254 s user, 0.015936 s system and 39,044 KiB peak RSS.  The independent
run used 31.872951 s wall and 39,108 KiB peak RSS.  Both result and assessment
artifacts are byte-identical.  Two earlier successful pairs are preserved in
the six-row isolation ledger; they predate the final operation counters.

- result: 39,954 bytes, SHA-256
  `e6b78d729e785b1831a061cd39e6e57f10c46da8b8e98418321dc72f56752f35`;
- result semantic evidence SHA-256:
  `8b82e8d7eda1e65737ffd0617103861a0ea9b102ef830b80a31057d4a56a17ed`;
- transfer assessment: 6,725 bytes, SHA-256
  `1925eb6bdc6dcd100d205ecc1451283e32cb85dfcbab828b27dac80b808537d7`;
- assessment semantic evidence SHA-256:
  `06ce7180c57727528755e5f6598990bd0f56739f0ece89956fc1dda5d76146df`;
- six-run isolation ledger: 12,729 bytes, SHA-256
  `0c0da8ef44648ac7bae35bd8910321b64e4ff4f054980fc3790f036d67b58356`.

```bash
cargo test --release --bin p256_transfer_anchor_screen
cargo clippy --release --bin p256_transfer_anchor_screen -- -D warnings
cargo build --release --bin p256_transfer_anchor_screen --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-transfer-anchor-round299-canonical-v3-ops \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_transfer_anchor_round299_20261007/isolation.jsonl -- \
  target/release/p256_transfer_anchor_screen \
  --round298 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_transfer_anchor_round299_20261007/transfer-anchor-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_transfer_anchor_round299_20261007/transfer-assessment.json
```
