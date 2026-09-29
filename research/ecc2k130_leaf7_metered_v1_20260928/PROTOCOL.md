# Held protocol: meter the degree-7 route to a paired ECC2K-130 leaf

**Status: draft protocol only.** No explicit degree-7 map, metering
producer, cost receipt, PDP relation, or logarithm is produced by this
PR. The frozen `run_metered.py` path is absent and its source hash and
commit are null. Implementing and independently reviewing an instrumented
producer is an additional release prerequisite; this checker cannot turn
self-declared counters into evidence. The structural producer in
[PR #810](https://github.com/aburan28/crypto/pull/810) is itself frozen with
`release_main_head: null`; it has no accepted map outcome. This successor
stays held until that producer is refrozen after its parent merge, its exact
head passes CI and peer review, and its producer, Fq replay, and extension
replay yield the required `PRODUCER_PASS`, `FQ_REPLAY_PASS`, and structural
`PASS`. Because that refreeze changes a source hash here, this protocol
must then be refrozen and reviewed **before** any metered outcome. The
`release_main_head: null` in [FROZEN.json](FROZEN.json) is intentional.

The question is narrow and falsifiable. Once the first descending
degree-263 quotient `L_C` is available, can its unique horizontal degree-7
bridge reach the paired leaf more cheaply than the already implemented
second degree-263 quotient? The frozen hypothesis is an **incremental**
bridge/direct cost ratio at most 0.90 and a **cold route** ratio below 1.00
on the same twelve full points, with all setup and extension arithmetic
charged. Success is an engineering improvement in reaching a second
representation. It is not a PDP-support, relation-rank, rho, complete
index-calculus, or ECC2K-130 DLP improvement.

## Facts, prior evidence, and limits

For `E: y²+xy=x³+1` over `F_(2^131)`, `tau²+tau+2=0` and the geometric
source endomorphism ring is the maximal order of discriminant `-7`.
The [degree-263 certificate](../ecc2k130_endo_ring_263_20260925/README.md)
computes `tau^131 = 8123678690673902378 +
38531015900842053623*tau`, classifies all 264 source kernel lines
(two horizontal and 262 descending), and checks two descending
coordinate-squaring orbits of length 131. By ordinary-isogeny-volcano
theory ([Sutherland](https://arxiv.org/abs/1208.5370)), a non-tau-stable
degree-263 line from the maximal source descends to conductor 263 and
endomorphism discriminant `-7*263² = -484183`. The Frobenius-generated
order has conductor `38531015900842053623`; its discriminant is not the
leaf endomorphism discriminant. Equal trace and curve order alone cannot
certify the exact leaf ring.

The [structural bridge note](../notes/ecc2k130/LEAF_HORIZONTAL_SEVEN_BRIDGE_20260925.md)
derives the unique horizontal degree-7 edge from a leaf, the norm-7
source action `alpha=1+2*tau`, and the commuting square
`psi * phi_C = iota * phi_Cprime * alpha`. These are mathematical
deductions from the order structure and quotient universal property.
They are **not yet an observed 7-map**. Its full-point pullback theorem
also says an isogeny preserves every labelled m-sum decomposition and
relation rank when the base is pulled back exactly. A conductor change
alone cannot yield more support. The [factor-base replication](../ecc2k130_factor_base_replication_20260925/README.md)
measured zero hits in the tiny exact degree-263, m=2 panels; those zeros
do not decide an implicit m>=3 PDP. This protocol does not run one.

## Frozen instance and input panel

Use the exact binary field polynomial and public `P,Q,r` from
[FROZEN.json](FROZEN.json) and the source hashes recorded there. Its
`r=680564733841876926932320129493409985129` and
`#E(F_q)=4r`; every Fq-rational degree-263 or degree-7 isogenous
leaf has that order. Require `[r]P=[r]Q=O` and the same identity after
both routes, rather than inferring point order from coordinates.

Use saved torsion seed **20260924** and exactly these projective lines:
`C=[1,0]`, the other-cycle witness `C2=[1,4]`, and
`Cprime=alpha(C)=[1,74]`. The saved twist-basis matrix has columns
`(186,51),(73,78)`; source `tau` acts as its negative on this twist
basis. Thus the required norm-7 line action is `A=I-2*M_twist`, not
`I+2*M_twist`. Require `A²=-7I mod 263`, two fixed horizontal
lines, two disjoint closed 131-cycles, `A` exchanging those cycles,
and `M_twist^k(C2)=Cprime` for exactly `k=23` in `0..130`.
The wrong-sign action is a failure control.

Transport exactly twelve full points, in this order:
`P,Q,P+Q,2P,target-0,...,target-7`. For `i=0..7`, define
`u_i=SHA256("leaf-seven-bridge-v1|i|u") mod r`,
`v_i=1+SHA256("leaf-seven-bridge-v1|i|v") mod (r-1)`, and
`target-i=[u_i]P+[v_i]Q`. The canonical eight-pair projection hash
is `b39da29224f3d503d2409195ae1586b2198dcda62a9e70e4f35a47e516aa22bd`.
These are the existing [#808 panel](../notes/ecc2k130/LEAF_HORIZONTAL_SEVEN_BRIDGE_20260925.md),
not selected after timing. The preflight recomputes every coefficient,
all 264 projective-line orbits, matrix sign, trace, order, conductor
fields, and input hashes without Sage or map construction.

## Structural admission before costing

The prior [#810 producer](../ecc2k130_leaf7_bridge_20260925/PROTOCOL.md)
must independently certify the 7-division factorization
`(degree 3)^1 * (degree 21)^1`, the irreducible cubic kernel and
its six order-7 exceptional points, the dual's six exceptions, infinity
and off-curve controls, and `dual*psi=[7]` in both directions. It must
produce a unique Fq codomain isomorphism and the full-point commuting
square on all twelve points, not just matching `j` or an x-coordinate.
The Fq replay and separate extension replay must both pass. The
opposite-cycle `[1,4]` map is a structural control for unique
23-conjugacy; it is not charged to either route arm.

The cost producer must bind its receipt to the exact hashes of the
structural producer, independent Fq replay, and extension replay.
The checker requires the Fq replay to report `FQ_REPLAY_PASS` on all
twelve labels and the extension replay to bind that exact Fq file.
A `MEASURED` receipt must archive both arms' ordered, labelled
full affine `(x,y)` arrays; the checker recomputes each digest and
compares every point with the independently checked parent outputs.
A `STOP` receipt may lack outputs and still remains STOP. These
receipt checks do not replace independent review of the metering
implementation.

## Two matched route arms and one unit

The frozen direct comparator is the saved-basis construction used in
#810: derive `Cprime` from the paid torsion basis and build the second
`BinaryVeluMap.from_generator` quotient, with its normalized dual
and source `alpha` on every panel point. A later faster direct
construction must be a separately frozen comparison; do not choose a
favorable baseline after seeing costs. Both arms supply the same
forward-and-dual capability and the same twelve oriented leaf outputs.

| Cost block | Exclusively charged work |
| --- | --- |
| Shared cold `S` | Recover the frozen 263-torsion basis from seed 20260924; build the first degree-263 map and normalized dual; materialize the twelve source points. **No** `phi_C(T)` panel evaluation is shared: the direct arm does not need one. Archived torsion is a verification reference, not a free cold input. |
| Direct incremental `D` | Derive `Cprime`; construct the second degree-263 map and dual; apply source `alpha` to the panel; normalize the model; evaluate `phi_Cprime(alpha(T))` on all twelve full points. |
| Bridge incremental `B` | Discover/factor the 7-kernel including extension arithmetic; construct the 7-map and dual; charge twelve `phi_C(T)` first-leaf evaluations; apply the 23 conjugacy squarings; normalize the model; evaluate `psi(phi_C(T))` on the same twelve full points. |

The only headline unit is an `F_q` multiplication equivalent:
`mul + w_sq*sqr`. Count every underlying Fq multiply and square,
including Sage polynomial/kernel discovery, extension arithmetic,
inversion internals, model conversion, validation needed by the route,
and all failed work. Record inversion calls separately as diagnostics;
do not add their internals twice. The frozen input hashes include the
field primitives, both 263-map/dual implementations, and the structural
7-map producer. A metering producer must be added in a new reviewed
commit; its own SHA-256 and commit must be frozen before release.

Measure `w_sq` on the **same host and code checkout** with eleven raw
batches of 10,000 SHA-fixed nonzero operand pairs labelled
`leaf-seven-cost-v1`. For batch `b`, index `i`, suffix `a` or
`b`, the operand is `1+int(SHA256("leaf-seven-cost-v1|b|i|suffix"))
mod (2^131-1)`, encoded as 17 big-endian bytes. Each multiplication uses `a*b` and each square uses `a²` from the same indexed pair. The ordered pair-stream
SHA-256 for each batch is frozen in the manifest. Archive each batch's
pair hash, operation count, multiply seconds, square seconds, and host
digest. The checker recomputes median per-operation times, square
weight, and both relative median absolute deviations from those raw
batches. It requires a host manifest with OS, architecture, CPU,
Python/Sage versions, checkout commit, and the exact field-primitive
hashes; calibration and all three children must bind its digest.
If either MAD exceeds 5%, or any shared/arm phase is unmetered, leave
both ratios null. Report CPU, wall, peak RSS and failed work as
secondary data. This first protocol makes no wall-time speed claim.

Compute `R_incremental=B/D` and
`R_cold=(S+B)/(S+D)` in that same calibrated unit, with the
common block paid **once per route**, not omitted or double-counted.
The cost receipt lists exactly the frozen phases. The
[cost-scope checker](cost_scope.py) rejects missing/extra categories,
nonexclusive or unmetered declarations, missing or mismatched archived
full-point arrays, an unbound Fq replay or parent certificate, raw
calibration inconsistencies, cross-host records, and caps. It recomputes both
ratios rather than trusting a reported headline. A green schema check
is not proof that the producer counted its internals; that code needs
separate pre-outcome review.

## Decision and stop rule

Do not start a cost child until structural admission passes. Stop each
common, direct, or bridge child at **300 wall seconds or 2 GiB peak
RSS**; preserve a STOP, timeout, crash, or OOM receipt, with no
favorable ratio inferred from a censored run. A STOP receipt does
not need target outputs or calibration arrays; a missing STOP receipt
remains unresolved and cannot become a pass. If structural checks fail, retire the proposed
bridge map and do not measure. If all checks pass but
`R_incremental>0.90`, retire this bridge-cost lead. If the
incremental ratio passes but `R_cold>=1`, report only a stage
diagnostic, not a full cold route win. Promote this *route engineering*
lead only if complete coverage, exact outputs, all caps, and both
frozen thresholds pass. Preserve negative and inconclusive receipts.

This is a **full cold cost of reaching and using the second leaf for
the fixed panel**, not a full ECDLP cost. No relation producer,
factor-base construction for a useful m>=3 attack, failed PDP queries,
linear solve, or same-Q rho reference is in this comparison. Those
remain unset. A later complete IC claim requires a separately frozen,
equal-target leaf-native versus exact-pullback comparison with all
those phases and a recovered scalar.

Preflight and checker self-test are safe before release:

```sh
python3 research/ecc2k130_leaf7_metered_v1_20260928/ci_preflight.py
python3 research/ecc2k130_leaf7_metered_v1_20260928/cost_scope.py --self-test
```

The preflight runs arithmetic and source-hash checks only. The checker
self-test uses synthetic counts only. Neither command imports Sage,
runs an isogeny, or dispatches a measured experiment.
