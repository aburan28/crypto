# Rotated normal-basis support: exact toy gate and n131 density admission

**Decision: the preregistered support gate passes on n=13.** On the exact
Koblitz curve E/F_(2^13), five distinct Frobenius-rotated summand spaces
cover 1,591 of 2,003 projected prime-subgroup targets, versus 61 for the
same-size repeated-F0 control. Six rotated spaces cover all 2,003 targets
for both independently selected normal generators, versus 85 or 369 for
the corresponding repeated control. The full-group histogram, all 18,501
positive full-point witnesses and all 45,595 negative full-point decisions
across the eight arms were independently replayed with different finite-field
multiplication and inversion code. This passes a *small-curve support*
hypothesis. It does not admit a practical S6/S7 solver, relation rank, an
n=131 PDP success rate, a Certicom logarithm, or a speedup over rho.

The experiment follows [GGMP §3.1](https://eprint.iacr.org/2020/1315.pdf):
`V_i = span(beta^(2^(m*j+i)))`, with `m*d <= n` and `F_i = tau^i(F_0)`.
The control repeats `F_0` for every summand; both policies use the same
number of point choices in each slot and the same projected log-column
source. A distinct translated space need not have the same PDP-support
behavior as a repeated space, even if its individual factor has the same
size. The literature's possible `m!` relation/symmetry advantage is
conditional on the cost of a complete decomposition solver; this
experiment measures only exact group-law support and lift density.

## Exact n=13 support

The field polynomial is `x^13+x^4+x^3+x+1` (`0x201b`), the curve is
`y²+xy=x³+1`, and the verified group order is `8012=4*2003`. The frozen
subgroup generator is `H=(4793,2429)`. Both beta values (3 and 7) were
selected by field-only normality rules before observing support. All
2,003 Q=[k]H, including O, were evaluated against each of the four
rational torsion shifts `Q+T`, `T∈{O,(0,1),(1,0),(1,1)}`. Each projected
count was also built separately by mapping every distinct full sum S to
`[4]S`, then equated to the sum of its four coset counts. The table counts
distinct targets, not labelled tuple multiplicity.

| beta | m | F_i size | Compressed `[4]F_0` points | Policy | Labelled tuples | Distinct full sums | Raw H hits / 2,003 | Projected H hits / 2,003 | Rotated / repeated projected support |
| ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 3 | 5 | 5 | 5 | Rotated | 3,125 | 1,799 | 404 | 1,591 | 26.082 |
| 3 | 5 | 5 | 5 | Repeated | 3,125 | 61 | 10 | 61 | 1.000 |
| 3 | 6 | 5 | 5 | Rotated | 15,625 | 6,307 | 1,659 | 2,003 | 23.565 |
| 3 | 6 | 5 | 5 | Repeated | 15,625 | 85 | 21 | 85 | 1.000 |
| 7 | 5 | 5 | 5 | Rotated | 3,125 | 1,799 | 404 | 1,591 | 26.082 |
| 7 | 5 | 5 | 5 | Repeated | 3,125 | 61 | 10 | 61 | 1.000 |
| 7 | 6 | 7 | 7 | Rotated | 117,649 | 8,012 | 2,003 | 2,003 | 5.428 |
| 7 | 6 | 7 | 7 | Repeated | 117,649 | 377 | 97 | 369 | 1.000 |

The n=13 nonzero target counts are one less than the projected counts in
every row because O is supported. The two beta-3 and beta-7 m=5 factors
are genuinely different: their F0 point sets share only `(0,1)` and their
1,591-member projected target sets intersect in 1,245 targets (union 1,937).
Equal aggregate counts therefore do not mean that the same targets were
chosen twice. The exact m=5 rotated miss sets contain 412 targets each,
with only 66 common misses. Both m=6 rotated rows have *zero* projected
negative H targets, so this n=13 panel cannot supply a negative-instance
S6 refutation benchmark. For beta=7/m=6, all 8,012 full-curve points occur
as sums; for beta=3/m=6, 6,307 full points occur, yet projection still
covers all H targets. The four torsion shifts are essential: for
beta=3/m=5, raw H support is 404 but `[4]`-projected H support is 1,591.

Target multiplicity is distinct from support. For beta=3/m=5 rotated, a
supported H target has 1–6 labelled witnesses after projection (median 2);
the repeated control has 1–221 (median 20) while covering only 61 H targets.
For every nonzero H target in this complete census, `Tr(x_Q)=0`,
`tau`-orbit size is 13, and the exact cofactor class `[2003]Q` is O. Those
cheap public features have **no within-H variation to stratify** easy
versus hard PDP targets here. Witness multiplicity varies, but computing
it is the decomposition problem, not a free target feature. A new
held-out hardness classifier would need other preregistered features and a
larger rung.

## n=131 structure, lift sample and count gate

The n=131 polynomial model is `x^131+x^13+x²+x+1`; it is field-isomorphic
to the public challenge field but does not import target coordinates. The
Weil trace recurrence from `#E(F_2)=4` independently gives
`#E(F_(2^131))=2722258935367507707729280517973639940516=4q` with
`q=680564733841876926932320129493409985129`. The frozen normal beta=3
has rank 131 and trace 1. A deterministic rational point with x=3 has
nonzero H=`[4]P`; independent group-law replay verifies `[q]H=O`,
`tau(H)=[lambda]H`, `[4]tau(P)=tau([4]P)` and `lambda^131=1 mod q`, where
`lambda=196511074115861092422032515080945363956`. Thus the structural
Frobenius transport can use one projected F0 log-column family. Its actual
size is not measured at n=131.

For each cell, 256 identical coefficient masks transported across all `m`
spaces passed exact Frobenius/liftability covariance. Another 16,384
unique, disjoint nonzero masks sampled the F0 lift condition; all 81,920
sample decisions were independently replayed. The Wilson intervals below
are descriptive under the frozen SHA-mask pseudorandom sampling model,
not rigorous full-space bounds or target-support estimates. The ideal
ceilings assume `|F_i|≈2^d`; the sample-calibrated ratios extrapolate
from the measured lift rate. Every value is only a *necessary tuple-count
ceiling* for uniform target coverage, not a predicted hit rate.

| m,d | Solvable x / 16,384 | Wilson 95% lift fraction | Ideal raw / full-curve ceiling | Ideal projected / H ceiling | Sample-calibrated projected ratio | Descriptive upper projected ratio | Ideal physical choices `m*2^d` | F0 log-column proxy `2^d` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5,24 | 8,070 | [0.48490, 0.50021] | 0.0004883 | 0.0019531 | 0.0018120 | 0.0019572 | 83,886,080 | 16,777,216 |
| 5,25 | 8,143 | [0.48935, 0.50467] | 0.015625 | 0.0625 | 0.0606530 | 0.0654706 | 167,772,160 | 33,554,432 |
| 5,26 | 8,228 | [0.49454, 0.50985] | 0.5 | 1 (ratio exceeds 1) | 1 | 1 | 335,544,320 | 67,108,864 |
| 6,20 | 8,117 | [0.48777, 0.50308] | 0.0004883 | 0.0019531 | 0.0018483 | 0.0020264 | 6,291,456 | 1,048,576 |
| 6,21 | 8,198 | [0.49271, 0.50802] | 0.03125 | 0.125 | 0.1255503 | 0.1375250 | 12,582,912 | 2,097,152 |

Under the preregistered descriptive 1% count-admission rule, m=5,d=24
and m=6,d=20 fail even at the Wilson-model upper rate; m=5,d=25 and
m=6,d=21 clear this *necessary* screen. m=6,d=21 is the first smaller-base,
higher-ceiling candidate, conditional on an implicit solver and memory
feasibility; m=5,d=25 is the fallback. The displayed F0 columns are nominal
proxies, not upper bounds: with two possible lifts for every nonzero x,
`|[4]F0| <= |F0| <= 2^(d+1)-1` (4,194,303 at d=21 and 67,108,863 at
d=25). Cofactor projection may create further duplicates. Merely
materializing a three-summand half at 24 bytes per tuple would project
about **221 EB** for m=6,d=21 and **907 ZB** for m=5,d=25 before any
deduplication; these are conditional storage projections, not lower bounds
for an implicit algorithm. No n=131 target decomposition or relation rank
was attempted.

## Costs, receipt and next falsifiable step

The successful runner ran 2026-09-25 11:45:34–11:45:56 UTC on the isolated
local Mac after the source/input freeze at commit `a6be357` and green
freeze-only CI. The toy producer's one-time field/curve/group/target setup
took 0.402 s and 40,134 counted point additions. Individual arm totals,
including factor construction, the exact histogram, projection, witnesses
and all 8,012 Q+T right-hand sides, took 0.080–0.803 s and 17,298–161,867
counted point additions. The n=131 shared field/normal/trace and lambda
point setup took 0.238 s; each lift-density cell took 0.754–0.838 s.
These timings are implementation and host diagnostics, not calibrated
ECC2K-130 attack costs. The separate independent replay took 8.435 s for
the exact toy arms and 6.750 s for density; it reports its own CPU, peak
RSS and arithmetic counters. Primary peak RSS stayed below 36 MiB and
verifier peak RSS below 50 MiB. All eight toy arms and five density cells
completed below their preregistered deadlines and RSS admission cap.

The [frozen protocol](PROTOCOL.md), [source/input hash map](FROZEN.json),
[evidence receipt](evidence/receipt.json) and
[compressed raw archive](evidence/raw.tar.gz) preserve every target row,
full-point witness index, factor list, covariance mask, density x decision,
phase cost and independent replay report. The archive SHA-256 is
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`
(2,680,225 bytes); the receipt SHA-256 is
`8b35cf4ea223b112d987abeed04a0ba8bf04e2f509ec0943eead87bdcdce9eac`.
The focused CI extracts it, verifies all raw hashes and repeats the
independent census and density replay; host-specific wall/CPU/RSS fields
are excluded from deterministic comparison.

The next bounded experiment should preregister an **implicit m=6
architecture and memory gate** at the m=6,d=21 n=131 count setting,
with m=5,d=25 as a fallback only if complete solver costs warrant it.
The merged [PR #763 solver-admission audit](ROTATED_M5_M6_SOLVER_ADMISSION_20260925.md)
shows that the current in-process Boolean polynomial engines are limited
to 64 variable bits while these direct systems start at about 126/125
summand bits, respectively; a canonical multiword or streaming exporter
and semantic corpus are prerequisites to any solver timing.
A separate exact PDP corpus/negative-refutation correctness gate can
start at n=13,m=5 and n=19,m=6 before solver timing: at n=19 the prime
subgroup has order 130,873, while each d=2 factor has at most seven
rational points, so `7^6=117,649<130,873` forces at least 13,224
projected-negative H targets before tuple collisions. No n=13,m=6
projected-negative corpus exists under these measured bases. Promotion
to an n=131 relation experiment requires a bounded implicit solver,
measured memory/work per SAT and UNSAT target, cofactor-correct group
witnesses, independent rank and scalar recovery, and matched same-Q rho.
