# Torsion-coordinate relation solves over a quadratic extension: cold DLP cost

Status: four frozen public elliptic-curve DLP instances measured and independently replayed on one shared host, L0. The direct Semaev and nonfactorable χ polynomial oracles in the original protocol remain unimplemented. This is a measured comparison of the explicit pair-table, χ orbit-index, and degree-three isogeny controls against same-target rho, with the original oracle question still open.

## Question, identity, and boundaries

The group being solved is the prime-order subgroup of `E(F_(p²))`, where `E:y²=x³+1`, `F_(p²)=F_p[s]/(s²−ν)`, and the public rational three-torsion point is `T=(0,1)`. The tested coordinate is `χ=(y+1)/(y−1)`, with `χ(P+T)=1/(1−χ(P))`. The public degree-three map is `φ:E→E':y²=x³−27`. [The frozen protocol](PROTOCOL.md), [selector amendment](PROTOCOL.md), [fixture records](FIXTURES.jsonl), [EC1 and ICV1 identity records](IDENTITIES.jsonl), and [paired run spec](RUN_SPEC.json) pin the experiment. The original χ polynomial decomposition oracle was not substituted silently: [the pre-run amendment](AMENDMENT_V3.md) labels the implemented χ **index** and exact raw-pair control.

The four tuples `(p,ν,r,h)` are `(61,2,19,192)`, `(271,3,61,1200)`, `(787,2,211,2928)`, `(3319,3,811,13584)`. Each source generator and target has a non-base-field coordinate. The destination is a separate explicit degree-three image, checked at 64 points per instance. The additive finite-field embedding degrees `ord_r(p)` are respectively `9,10,105,810`; relative to `q=p²`, they are `9,5,105,405`. These degrees screen a possible pairing route into `F_(p^ord_r(p))`; no pairing construction or multiplicative-field DLP was timed here. The measured reference is the ECDLP itself over `F_(p²)`, with the exact same `G,Q,r` in negation rho, duplicate negation rho, and charged six-class rho.

The predeclared generic two-class count floor is `S=√(π/4)=0.886226925` group operations divided by `√r`. The table's unit is instead `S_cpu_gae = cold current-thread CPU nanoseconds / measured source-addition nanoseconds / √r`. It charges native roots, field arithmetic, allocations, table construction, relation solving, and verification, but the conversion is local to this host and gives **no universal group-operation floor**. Its floor-ratio column is a local equivalent-cost comparison with the fixed generic reference, not a proof of a generic-group bound. The admissible speed ratios are whole-cold-solve paired CPU ratios; no phase-only speed ratio appears.

## Complete cold solves on the frozen targets

The [corrected V3 records](RUN_V3.jsonl) have three paired rounds per size and arm. The [native summary](SUMMARY_V3.jsonl) calculates medians and the minimum–maximum run range. The ratio to negation rho is the median of same-target, same-round total CPU ratios, so it need not equal the ratio of displayed medians. `cap` rows spent their full query budget but have no answer or floor/reference ratio. The table uses one CPU-converted unit for every variant. `χ index` means exact projected-pair membership indexed by the χ orbit; every proposed hit is checked by full subgroup point equality.

| `p` | `r` | arm | cold `S_cpu_gae`, median [min, max] | local ratio / fixed floor | paired cost / negation rho | scalar replay |
| ---: | ---: | :--- | ---: | ---: | ---: | :--- |
| 61 | 19 | chi-orbit-index | 195.05 [194.2, 199.78] | 220.09 | 9.74 | 3/3 |
| 61 | 19 | isogeny-projected-pair | 172.81 [169.15, 523.88] | 194.99 | 8.49 | 3/3 |
| 61 | 19 | projected-pair | 164.19 [161.49, 188.12] | 185.27 | 9.44 | 3/3 |
| 61 | 19 | raw-pair | 136.94 [129.22, 211.71] | 154.52 | 9.18 | 3/3 |
| 61 | 19 | rho-negation | 19.93 [14.92, 72.41] | 22.49 | 1 | 3/3 |
| 61 | 19 | rho-negation-aa | 16.27 [12.56, 36.04] | 18.36 | 0.63 | 3/3 |
| 61 | 19 | rho-six-class | 23.93 [20.83, 44.65] | 27 | 1.2 | 3/3 |
| 271 | 61 | chi-orbit-index | 156.62 [148.01, 188.97] | 176.73 | 12.59 | 3/3 |
| 271 | 61 | isogeny-projected-pair | 140.38 [134.08, 143.92] | 158.4 | 11.02 | 3/3 |
| 271 | 61 | projected-pair | 138.43 [129.92, 146.06] | 156.2 | 11.47 | 3/3 |
| 271 | 61 | raw-pair | 572.3 [564.88, 592.09] | — | — | 0/3; cap |
| 271 | 61 | rho-negation | 12.74 [11.22, 15.02] | 14.37 | 1 | 3/3 |
| 271 | 61 | rho-negation-aa | 9.48 [9.33, 13.06] | 10.7 | 0.83 | 3/3 |
| 271 | 61 | rho-six-class | 16.93 [13.21, 54.09] | 19.11 | 1.51 | 3/3 |
| 787 | 211 | chi-orbit-index | 190.65 [189.66, 208.76] | 215.12 | 19.47 | 3/3 |
| 787 | 211 | isogeny-projected-pair | 160.28 [157.37, 169.73] | 180.86 | 15.83 | 3/3 |
| 787 | 211 | projected-pair | 159.01 [152.31, 190.68] | 179.42 | 17.78 | 3/3 |
| 787 | 211 | raw-pair | 866.26 [863.4, 910.86] | — | — | 0/3; cap |
| 787 | 211 | rho-negation | 10.72 [6.97, 11.62] | 12.1 | 1 | 3/3 |
| 787 | 211 | rho-negation-aa | 8.55 [5.96, 12.01] | 9.65 | 0.86 | 3/3 |
| 787 | 211 | rho-six-class | 13.26 [9.86, 15.64] | 14.96 | 1.46 | 3/3 |
| 3319 | 811 | chi-orbit-index | 284.96 [270.42, 531.87] | 321.54 | 51.89 | 3/3 |
| 3319 | 811 | isogeny-projected-pair | 217.45 [214.89, 355.69] | 245.37 | 48.71 | 3/3 |
| 3319 | 811 | projected-pair | 202.56 [199.07, 204.41] | 228.56 | 27.74 | 3/3 |
| 3319 | 811 | raw-pair | 1961.89 [1707.48, 2355.93] | — | — | 0/3; cap |
| 3319 | 811 | rho-negation | 7.3 [4.4, 10.25] | 8.24 | 1 | 3/3 |
| 3319 | 811 | rho-negation-aa | 6.26 [5.86, 8.77] | 7.06 | 0.86 | 3/3 |
| 3319 | 811 | rho-six-class | 9.53 [8.66, 10.16] | 10.76 | 1.39 | 3/3 |

The corrected run contains 84 arm rows, 75 independently replayed scalar recoveries and nine capped raw-pair rows. [The independent arithmetic replay](REPLAY_INDEPENDENT_V4.json) checked `[k]G=Q` for every recovered scalar and checked the cap statuses against the byte-hashed run. The record SHA-256 is `ca9a0ec7e31de59c58a78ae5c269350b37963b62ff9de815c805c36c2ec26693`; the fixture SHA-256 is `1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`. Relation rows pass exact source group equality before elimination; the decomposition oracle and collector receive `G,Q,r,seed` through the context. The framework runner receives the known scalar solely for its final equality assertion after recovery; the independent arithmetic verifier checks `[k]G=Q` again.

At `r=811`, the projected-pair median cold cost is `202.56 S_cpu_gae`, the χ-index median `284.96`, the isogeny-control median `217.45`, negation rho `7.30`, and six-class rho `9.53`. The same-round projected/rho and χ/rho ratios are `27.7369` and `51.8898`, respectively. The χ index, projected pair table, and isogeny table used the **same** search outcome, rank, and recovered scalar in every paired `p×round` case; the extra indexing and map work did not buy relation yield in these fixed rows. At each size the χ-index median exceeded the projected-pair median. The three-round ranges show noise, notably the `r=811` χ maximum `531.87` and isogeny maximum `355.69`; these are descriptive L0 comparisons, not isolated wall speed claims. The six-class implementation is charged, but its measured cost is not a certificate that it is the best automorphism rho implementation.

The [precommitted replay protocol](INDEPENDENT_REPLAY_PROTOCOL.md) led to a second native implementation of `F_(p²)` arithmetic and the chord/tangent group law. It shares neither the run's field operations nor its scalar routine. It checked generator and target order, unique coverage of all 84 `(p,round,arm)` slots, exact seeds and status, and the same 75 scalars. The independent receipt SHA-256 is `785c5e615a18742302f05f6c2a219b260d48d5e2b63e984cec945222bd51ac5c`. The first independent [V3 receipt](REPLAY_INDEPENDENT_V3.json) checked arithmetic and slot integrity but did not assert the frozen input hashes; V4 adds that assertion and a one-byte [negative control](REPLAY_HASH_NEGATIVE.stderr) failed as expected. The earlier same-arithmetic [replay](REPLAY_V3.json) is retained separately. This strengthens arithmetic and row-integrity verification; it does not independently reimplement relation collection.

## Finite support and cost audit

The [precommitted census protocol](SUPPORT_PROTOCOL.md) and [exact native support records](SUPPORT.jsonl) enumerate every unordered signed-base pair. For `p=61,271,787,3319`, the raw exact subgroup support under shifts by `0,T,2T` has sizes `[2,0,0]`, `[0,0,0]`, `[0,1,1]`, `[0,0,0]`; the union sizes are `2,0,2,0`. Thus the `p=787` shifts do add **two** raw subgroup support points on this exact base, but its frozen target is in none of those supports. The projected supports are `18/18`, `58/60`, `196/210`, and `742/810` nonidentity subgroup points. The degree-three image supports have those same cardinalities. All `128,128,450,1682` enumerated nonidentity pair sums passed the χ torsion-action identity. The χ orbit key counts are `6,20,70,270`. These finite counts explain why the raw exact-pair control caps at three sizes and why shifts do not create a new projected pair sum; they do not test a χ polynomial-system oracle.

Each cold relation run includes the native line lifts and square-root tests, cofactor projection `[h]`, signed factor-base construction, pair table, all failed and successful queries, exact relation checks, rank and modular elimination, scalar recovery, and internal scalar verification. The isogeny arm additionally evaluates its public map on the projected base, generator, target, and query controls; the six-class rho interval includes eigenvalue discovery. The per-phase `factor_base`, `decomposition.setup`, `decomposition.cost`, `linear_algebra`, and `verify` objects and field-multiplication/inversion counters are in each [raw V3 row](RUN_V3.jsonl). For example, the `p=3319`, round-zero projected solve recorded `1043.695`, `1719.959`, `895.149`, `3.111`, and `13.000` counted GAE in those phases, respectively; its cold CPU was `913208 ns`, and it recorded `90110` base-field multiplications and `3689` inversions. That counted decomposition omits native costs represented by the CPU interval, so the primary table uses CPU conversion. The measured source-addition calibrations were `169.76625`, `147.45333`, `162.66167`, and `156.87292` current-thread ns for the four fields. Calibration is outside each method's cold interval.

The pair-table entry counts are `128,128,450,1682`. Their serialized 24-byte payload lower bounds are `3072,3072,10800,40368` bytes, excluding hash-table capacity, point objects, allocator overhead, and linear-algebra storage. A per-arm peak RSS trace was **not** captured, so an equal-memory or memory-adjusted comparison remains open. Setup reuse count is one target; no amortized column is claimed. No direct Semaev resultant, χ descended polynomial system, quotient-invariant polynomial oracle, or finite-field pairing transfer is part of these timings. This is a complete CPU price for the implemented arms only, not a complete price for the original oracle proposal.

## Corrected runs, proof scope, and decision

The first 84-row output [V1](RUN_V1.jsonl) completed, but its `/usr/bin/time -l` wrapper exited 1 when sandboxed `sysctl kern.clockrate` failed; [the stderr](RUN_V1.stderr) is retained. [V2](RUN_V2.jsonl) exited cleanly and its scalar replay passed, but the driver timed JSON serialization of the much larger relation reports inside the method interval. [The correction note](CORRECTION_V4.md) fixes that interval for V3 without rewriting V1/V2. The corrected version stops the clock before report serialization, preserving the same fixtures, seeds, caps, and algorithms. The first selector screen [V1](SCREEN_V1.jsonl) also remains preserved; [the pre-solve v2 amendment](PROTOCOL.md) imposed rising subgroup sizes before the paired run.

The algebraic identity and conditional support theorem are documented in the [prior proof map](../torsion_orbit_20261013/PROOF_MAP.md): `torsion_orbit_pair_support` in `OrbitSupport.lean` proves the quotient-support equality assuming a homomorphism and kernel membership; `cubicAction_order_three` in `EquivariantOrder.lean` proves the stated integer-linear action has order three. Both are reported there as compiled without custom axioms. The finite-field curve-law bridge is supplied for these four instances by the native point checks, not by a new Lean theorem. The χ group-law and Möbius identities have the existing certificates `IDC1h41adc1f4e1c5522c` and `IDC1ha0c3d83775e32754` in the [certificate index](../torsion_orbit_20261013/CERTIFICATES.md). The exact support counts and measured CPU ratios are computational observations, not universal mathematical claims. The curve-level Lean formalization remains partial; this run does not add a new general theorem.

**Decision:** this concrete χ orbit index does not improve complete-solve cost over the unindexed projected pair table or matched rho on the four frozen groups. Classification: an **engineering regression** for the index; an **accounting correction** from V2 to V3 timing. The two extra raw support points at `p=787` are a precise torsion-fiber lead, but they do not change projected support or this target's outcome. No exponent, sub-rho, deployed-curve, or key-recovery claim follows. The more interesting possible mechanism is the still-unimplemented nonfactorable χ polynomial oracle; the [next gate](NEXT_GATE.md) specifies how to distinguish it from a mere index or quotient relabelling.

The [visual](GRAPH.svg) shows every arm's median and range on a logarithmic `S_cpu_gae` axis, with capped raw rows marked separately. The [run report](reports/run-20261014T1541Z/RUN_REPORT.md), [PDF](REPORT_V3.pdf), and [benchmark copy](../../../benchmarks/20261014-torsion-extension/) carry the same data. No canonical algorithm frontier or security graph is promoted because these rows do not establish a new measured advantage; this run adds a dedicated result graph.

### Requirement status

- Frozen ordinary `F_(p²)` instances, explicit three-torsion and degree-three map, source-target identities: **verified complete** for the four recorded instances by native freeze checks and `IDENTITIES.jsonl`.
- Cold pair-table, χ-index, isogeny-control, matched rho and six-class rho, exact scalar replay: **verified complete** at L0 for the measured arms; `RUN_V3.jsonl` and `REPLAY_INDEPENDENT_V4.json`.
- Finite exact torsion-shift and projected support: **verified complete** for the four frozen bases; `SUPPORT.jsonl`.
- Direct Semaev, χ polynomial-system and quotient polynomial oracles from the original protocol: **not attempted**; no implementation or oracle timing exists.
- Per-arm peak memory, L2 isolation, independent-host replay, and portable full primitive accounting: **partial/open**; the CPU conversion and counted records are local.
- Pairing route to a multiplicative extension-field DLP: **screened only** by embedding degree; pairing construction, transfer, multiplicative solver and full costs are open.
