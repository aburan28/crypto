# The exact degree-83 signed-orbit base and four-point witness

Descriptive source review, 2026-10-12. This note explains the existing
implementation and retained verification receipts; it makes no new algebraic
discovery claim and runs no arithmetic checker during the timed panel.

## The instance and the three spaces

The exact model is `icv1-f2m83-t6151469093347-9eba6cd5`:

\[
E:\quad y^2+xy=x^3+x^2+1,\qquad
F=\mathbb F_2[z]/(z^{83}+z^7+z^4+z^2+1).
\]

The field has 83 coordinates over F2. The DLP subgroup has prime order
`r=8569786107849059`, approximately 53 bits; its cofactor is
`1128547018`. Field degree 83 is not subgroup bit length.

The retained base has 600 representatives and 99,600 distinct usable signed
orbit points. It is a **point-defined union of orbits**, not a small
83-coordinate linear subspace V. Its vectors of relation coefficients live
in `(Z/rZ)^600`; those are a third space, separate from the field coordinates
and the actual elliptic-curve point set.

The 600 formal coefficient coordinates do **not** mean 600 independent
group generators. The rational prime-order subgroup used here is
`H=<G>`, cyclic and one-dimensional over `F_r`. Its representatives obey
`R_j=[L_j]G`; the coefficient-to-point map is

\[
\Phi:(\mathbb F_r)^{600}\longrightarrow H,\qquad
v\longmapsto\sum_j[v_j]R_j=[v\cdot L]G.
\]

The formal coordinates name 600 unknown scalar values L_j, each relative
to the same G with known logarithm one. Guided preparation solves a
600-column **inhomogeneous** system to recover these values; it is not a
declaration that the points R_j form an independent group basis. This
description concerns H, not the full geometric r-torsion over an algebraic
closure. Source: coefficient accumulation and known right-hand side at
[row construction:1919](../../../examples/koblitz_orbit_dlp_fast.rs#L1919),
augmented echelon system at [line 808](../../../examples/koblitz_orbit_dlp_fast.rs#L808),
and scalar recovery relative to G at [line 2280](../../../examples/koblitz_orbit_dlp_fast.rs#L2280).

Identity correction: the historical alias
`icv1-f2m83-t6151469093347-cdcc5432` belongs to modulus low terms
`[0,1,2,45]`, not this retained fixture's `[0,2,4,7]`. The explicit-field
identity computed by the fresh strong producer is
`ICV1:f2m-83-a5fd0d39:6151469093347:9671406556910881928556062:0x1:-7:unk:r:9eba6cd5d359`.
The original native archive receipt and its PDF contained the incorrect
historical label; their exact field, subgroup, G, Q and arithmetic checks
remain valid. Frozen evidence is preserved and corrected by the separate
[metadata erratum](CURVE_ID_ERRATUM.md), not rewritten.

Sources: [fixture](../../../experiments/koblitz-single-target-n83-20261006/frozen/fixture.json)
and [native archived receipt](receipt.json). The field, subgroup and base
checks in that receipt have already passed independent native arithmetic
replay. The frozen upstream snapshot is
`a8194d616ed8c04489ec46177d369f1b1670360a`.

## Signed Frobenius invariants and vector labels

Here `tau(x,y)=(x²,y²)` and `-P=(x,x+y)`. On the prime subgroup,
`tau(P)=[lambda]P`, where `lambda=2998004823950457`,
`lambda²-lambda+2=0 (mod r)`, and lambda has exact order 83.
The signed orbit of each nonidentity subgroup representative R_j has 166
members:

\[
P_{j,k,\epsilon}=\epsilon\tau^k(R_j),\qquad
\log_G(P_{j,k,\epsilon})=\epsilon\lambda^k L_j\pmod r,
\quad L_j=\log_G(R_j).
\]

The coefficient vector attached to this point is
`v(P)=epsilon lambda^k e_j`, with e_j the j-th standard unit vector in
`(Z/rZ)^600`. In the stored base the ordering is `166j+2k+s`, with s=0
positive and s=1 negative. Negation changes y but leaves x unchanged;
therefore an x-coordinate root alone never chooses the correct sign.

The independent archived checker verified every representative's subgroup
membership, all 99,600 coordinates and labels, distinctness and orbit
closure. It also checked `tau(G)=[lambda]G`. Source: [check.cpp:69](check.cpp#L69),
[check.cpp:76](check.cpp#L76).
The existing builder explains cofactor projection, representative selection,
signed images and labels at [koblitz_orbit_dlp_fast.rs:1852](../../../examples/koblitz_orbit_dlp_fast.rs#L1852).
The actual retained run loads its full point-defined header rather than
reconstructing a base: [same file:2016](../../../examples/koblitz_orbit_dlp_fast.rs#L2016).

## What the root index does

The characteristic-two third summation polynomial is

\[
S_3(x_1,x_2,u)=(x_1x_2+x_1u+x_2u)^2+x_1x_2u+1.
\]

For fixed x_1,x_2 this is a quadratic `A u²+B u+C`, with
`A=(x_1+x_2)²`, `B=x_1x_2`, `C=(x_1x_2)²+1`. Existing half-trace
arithmetic solves it; degenerate linear and single-root cases have separate
branches. Source: [koblitz_fast_arith.rs:301](../../../src/cryptanalysis/koblitz_fast_arith.rs#L301)
and its wide implementation at [line 598](../../../src/cryptanalysis/koblitz_fast_arith.rs#L598).

The index considers ordered representative pairs and relative Frobenius
shifts `(j,l,k)`. It solves `S3(x(R_j),x(tau^k R_l),u)=0`, retains the two
candidate roots, and stores their Frobenius-canonical normal-coordinate
keys with the source pair and canonical shift. The archived summary records
29,880,000 states and 29,878,182 distinct root-table entries; those are
stored-state counts, not independently verified log results or useful-row
counts. Source: [index construction:1637](../../../examples/koblitz_orbit_dlp_fast.rs#L1637)
and the archived `ic_summary.jsonl` R1 row.

The normal basis is an arithmetic/index representation. This particular
wide implementation deterministically searches increasing polynomial words
for the first normal element, builds inverse coordinate maps, and makes
squaring a cyclic rotation. It returns the least numeric rotation and its
shift. It does not change the point set or introduce new linear factor-base
membership. Sources: [normal-basis construction:936](../../../examples/koblitz_orbit_dlp_fast.rs#L936),
[canonical key:969](../../../examples/koblitz_orbit_dlp_fast.rs#L969).
Retain the construction/source identity or exact maps when exporting keys.

For a target Q, extraction rotates a first indexed pair root u, solves
`S3(u,x(Q),v)=0`, and looks up v's canonical key as a second indexed pair
root. The second pair's absolute phase is reconstructed as
`stored_shift-partner_shift (mod83)`; its relative phase is then added.
This yields four candidate x coordinates. The root key match is only a
candidate. `lift128` enumerates actual base points above those abscissae
and accepts only a full point identity `P1+P2+P3+P4=Q`. Sources:
[scan and phase transport:1746](../../../examples/koblitz_orbit_dlp_fast.rs#L1746),
[phase reconstruction:1789](../../../examples/koblitz_orbit_dlp_fast.rs#L1789),
[point lift and re-add:1693](../../../examples/koblitz_orbit_dlp_fast.rs#L1693).

## The retained four-vector witness

The three archived target repeats give the same point indices:
`[87298,27061,36601,66281]`. Their already verified orbit labels are:

| Point index | Column j | Shift k | Sign | Stored scalar coefficient |
|---:|---:|---:|:---:|---:|
| 87298 | 525 | 74 | + | 2610222652704134 |
| 27061 | 163 | 1 | - | 5571781283898602 |
| 36601 | 220 | 40 | - | 5545029831725595 |
| 66281 | 399 | 23 | - | 6086825800742164 |

Thus the target coefficient vector is

\[
v_Q=\lambda^{74}e_{525}-\lambda e_{163}
-\lambda^{40}e_{220}-\lambda^{23}e_{399},
\qquad d=v_Q\cdot L\pmod r.
\]

This table is an extraction from the retained header and target row, not a
new timed computation. The native archived check verified the four-point
sum, each x code and both pinned intermediate pair-sum abscissae. Source:
[archived IC target](../../../experiments/koblitz-single-target-n83-20261006/runs/N83A1K600W12We202610061R21R1/ic.jsonl#L1),
[native witness check:83](check.cpp#L83).
The existing recovery computes the dot product of stored coefficients with
the solved factor logs and then checks `[d]G=Q`: [recovery:2280](../../../examples/koblitz_orbit_dlp_fast.rs#L2280).

## Guided rank preparation and what is still missing

Guided collection picks an unresolved column j and a known setup scalar a,
decomposes `T=[a]G-R_j` into four base points, and therefore obtains

\[
L_j+\sum_{i=1}^{4}\epsilon_i\lambda^{k_i}L_{j_i}=a\pmod r.
\]

It accumulates coefficients of repeated columns and performs incremental
echelon elimination modulo r. A dependent row does not increase rank.
The solver requires a pivot in every one of the 600 columns before
back-substitution. Sequential collection targets the next missing pivot;
parallel collection uses per-column deterministic setup scalars and admits
completed rows as they arrive. Sources: [row construction:1919](../../../examples/koblitz_orbit_dlp_fast.rs#L1919),
[guided parallel relation:2131](../../../examples/koblitz_orbit_dlp_fast.rs#L2131),
[guided sequential relation:2190](../../../examples/koblitz_orbit_dlp_fast.rs#L2190),
[echelon insertion and solve:817](../../../examples/koblitz_orbit_dlp_fast.rs#L817).

The archive's rank=600 and zero-failure counters are producer assertions;
the summary rows do not retain every rank query, decomposition witness,
coefficient row, echelon matrix or solved L vector. The independent checker
can certify the final point witness and scalar answer without proving that
entire preparation transcript. Fresh target/strong producer JSON can support
the same final-answer replay and outcome checks, after execution completes;
it does not fill this full-rank/factor-log transcript gap by itself.

For a complete algebraic preparation certificate, retain all admitted
relations with setup scalars, point indices and labels; the matrix digest and
rank evidence; the solved 600 factor logs; and independent `[L_j]G=R_j`
replays. Keep its reusable setup cost separate from target extraction.
Independent correctness replay does not certify wall-time isolation,
resource equivalence, timing accounting, performance across independent
targets or a speedup. No checker was launched while the fresh timed panel
was active.

## How the strong rho comparator uses the same invariant

Strong rho works on the same exact field, subgroup, G, Q and scalar
Frobenius action. Its state has **two scalar coefficients**, rather than
600 factor-log coordinates:

\[
P=[a]G+[b]Q,\qquad a,b\in\mathbb Z/r\mathbb Z.
\]

When canonicalization chooses `C=epsilon tau^k(P)`, it applies the same
multiplier to both scalars:

\[
(P,a,b)\longmapsto
(C,\epsilon\lambda^k a,\epsilon\lambda^k b).
\]

This preserves the point/coefficient identity. This is a forward shift
convention: do not substitute the inverse exponent used when reconstructing
an original point from a stored representative. Source: [walk state:414](../../../src/cryptanalysis/koblitz_strong_rho.rs#L414),
[canonicalization and scalar transport:1056](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1056).

The routine converts x and y into normal coordinates, chooses the least x
rotation, chooses between y and y+x, and returns the polynomial-coordinate
point with transported coefficients. Tied x rotations use an exhaustive
signed-coordinate choice; infinity is returned unchanged. The underlying
normal field basis remains an 83-dimensional F2 coordinate representation,
as in IC, rather than a new factor base or an additional logarithm system.
The strong implementation uses a least-rotation selector based on runs of
zero bits in the cyclic coordinate word, with an exact canonical-key result;
the retained wide IC target path uses its full least-rotation scan. Both
implementations search deterministic normal elements for their exact
field. Do not infer timing equivalence merely from the shared mathematical
representation. Sources: [rho normal basis:608](../../../src/cryptanalysis/koblitz_strong_rho.rs#L608),
[rho least rotation:652](../../../src/cryptanalysis/koblitz_strong_rho.rs#L652),
[IC normal basis:936](../../../examples/koblitz_orbit_dlp_fast.rs#L936),
[IC target canonicalization:1784](../../../examples/koblitz_orbit_dlp_fast.rs#L1784).

Here a cyclic **zero run** is a property of the coordinate word. It says
nothing about the partition of elapsed time into phase timers. Both fresh
producers retain positive timer residuals. In formal matched pair 1, rho
records `online_interval_ns=4516159970`, phase sum `4516159569` and
residual `401 ns` (`0.000401 ms`). IC records `target_ms=6437.368839`,
five-phase sum `6437.338883` and residual `0.029956 ms`. The complete
online intervals include these residuals; no exact phase partition is
claimed. Retain them explicitly in the report rather than silently dropping
them or describing either producer's timing partition as exact. Sources:
[rho formal pair 1](../records/panel/matched_1_rho/producer.json),
[IC formal pair 1 semantics](../records/panel/matched_1_ic/execution.json).

The r-adding jump is a prepared known multiple `[s_i]G`, so an ordinary
step adds it and updates a to a+s_i before signed-orbit canonicalization.
The current fresh producer prepares this G-only jump table before Q's
online interval. Target-dependent starts have the form `[c]G+Q`, generated
using a random stride of G; b initially equals one. The default advances
32 lanes in lockstep, sharing an inversion across the independent additions;
32 lanes do not mean 32 CPU worker threads. Sources: [G-only jumps:1090](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1090),
[target starts:1110](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1110),
[fresh preparation boundary](../../../examples/n83_strong_point_eval.rs),
[batched additions:1219](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1219).

At a distinguished-point collision, let the current stored pair be
`(a1,b1)` and the earlier hash-table hit be `(a2,b2)`. The identical
canonical point gives

\[
(a_1-a_2)+(b_1-b_2)d=0\pmod r,
\qquad d=\frac{a_2-a_1}{b_1-b_2}\pmod r.
\]

The code first attempts to invert `b1-b2`. On the prime subgroup, a
noninvertible denominator is zero; that meeting is counted as a wasted
merge and supplies no recovered scalar. The lane is already retired at
the distinguished point, so it subsequently receives a fresh start. A
nonzero-denominator candidate is accepted only after `[d]G=Q`; an invalid
candidate increments failed-collision accounting. The fresh wrapper then
performs a separate generic-arithmetic scalar replay inside its online
interval. Sources: [collision recovery:1183](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1183),
[independent producer replay](../../../examples/n83_strong_point_eval.rs).

Fruitless-cycle handling is distinct from the zero-denominator collision.
The walker detects a next point equal to its current point or one of its
four retained previous points. With the default `escape_fruitless=false`,
it counts a fruitless cycle and abandons that lane. An optional escape
instead doubles the next point and both coefficients, canonicalizes again,
and clears recent history. This option is not enabled by the fresh
producer's defaults. Each lane is also retired after `8·2^dp_bits` steps;
a separate overall step cap returns no outcome, which remains a failure
record rather than a verified solve. Default `dp_bits=4`, lane count32 and
step-cap factor2000 are existing configuration values, not new measured
advantages. Sources: [defaults:443](../../../src/cryptanalysis/koblitz_strong_rho.rs#L443),
[caps:1128](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1128),
[cycle policy:1258](../../../src/cryptanalysis/koblitz_strong_rho.rs#L1258),
[fresh producer default parameters](../../../examples/n83_strong_point_eval.rs).

The comparison therefore shares the signed-Frobenius coefficient invariant
and field representation, while doing different target work. IC first pays
for its reusable root index and 600-column factor-log preparation and then
extracts a four-point witness. Rho prepares a small known-G jump table,
builds a target-dependent distinguished-point table and solves a collision
equation in two coefficients. Rho does not construct the IC root index or
solve its 600 unknown factor logs. Both final-answer paths require scalar
replay. This description explains the existing methods; it certifies no
fresh timings, asymptotic improvement or new algebraic identity.
