# Prebuilt factor bases at degree 83: review for m = 2, 3, and 4

Date: 2026-10-16. This review uses the exact primary curve
`icv1-f2m83-tm6151469093347-debefd74`, whose prime subgroup order is
`2417851639230796216685689`. The `a=1` degree-83 curve has a different,
53-bit subgroup and is discussed separately. Degree alone does not identify
the factor-base workload. The new local measurement follows [PROTOCOL.md](PROTOCOL.md)
and is reported in [RESULT.md](RESULT.md).

## What has already been built and retained

The earlier [native degree-83 sweep](../koblitz_n83_factor_base_sweep_20261008/RESULT.md)
constructed 57 signed-Frobenius bases, including the 54-object pilot panel,
with 45 distinct point sets. It checked 3,617,306 point records and 21,791
representatives with a separate generic arithmetic path. The object bytes
were uploaded to S3, fetched and hash-matched; the study's later standard
subspace archives brought its documented archive count to 74. These are
actual prebuilt bases, with curve, subgroup, source seed, orbit labels and
content hashes. The stored objects are reusable inputs to a relation
producer once that producer's representation and correctness gates pass.

The [lifecycle diagram](symmetry-flow.svg) marks the archived construction
and replay steps separately from the still-open natural relation oracle.

The signed-Frobenius closure is the useful symmetry here. The subgroup
Frobenius action is multiplication by a fixed scalar `lambda`, so one
representative encodes its 83 Frobenius phases and both signs with known
column coefficients. For the retained full orbits, `K` representatives
encode `166K` distinct points but only `K` relation columns. This is the
same mechanism analyzed by Galbraith, Granger, Merz and Petit in
[On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf),
whose rank and solver assumptions must be checked for the actual base.

There are two importantly different structured choices in the repository:

| Family | Primary-curve stored example | Algebraic advantage | Open cost |
| --- | ---: | --- | --- |
| Explicit signed-Frobenius orbit union | K=2,048; 339,968 points | 2,048 columns and exact compact orbit labels | Finite-domain encoding, natural relation yield, useful rank |
| Projected dimension-12 source, then Frobenius/sign closure | 2,001 columns; 332,166 points | Source coordinates and phase structure as well as orbit folding | Phase-aware relation solving and full cold pipeline |
| Standard dimension-18 source subspace | 130,722 sign-folded columns; 261,444 points | Direct small source-coordinate domain for Weil descent | Many more columns; cofactor-aware lifting and natural rank |

The [dimension-12 object](../koblitz_n83_factor_base_sweep_20261008/verification/structured-d12-orbit-20261010/RESULT.md)
is 10,058,074 compressed bytes and passed source-point, subgroup, orbit,
label and S3 round-trip replay. Its [earlier exact census](../f6_n83_mixed_base_20261004/ORBIT_RESULT.md)
records 4,054 source points becoming 332,166 closure points while the
column count changes from 2,027 sign-only to 2,001 signed-Frobenius. A
phase-aware five-summand SAT model for it was constructed under a 4 GiB
cgroup, but its first bounded search ended at `UNKNOWN_solver_cap`. This
does not measure its natural yield. The [standard-subspace frontier](../koblitz_n83_factor_base_sweep_20261008/verification/standard-d1-d17-s3-20261015/RESULT.md)
also has stored and checked point sets, but cofactor projection generally
moves a point outside the original x-subspace; its solver needs the source
to projected-column map.

## What arity changes on the primary curve

For `B` distinct subgroup points and repetition allowed, the exact number
of unordered m-summand multisets is `M = binomial(B+m-1,m)`. Every multiset
has one group sum. For a uniform nonidentity target, the full-smooth hit
probability is at most `min(1,M/(r-1))`. This bound counts no actual
relations; repeated sums and cancellation can make true coverage lower.
It also assumes uniform query targets. The table copies the exact-count
screen in [SUPPORT_MOMENTS.md](../koblitz_n83_factor_base_sweep_20261008/SUPPORT_MOMENTS.md)
and its [machine receipt](../koblitz_n83_factor_base_sweep_20261008/pilot-01/support-moments.json).

| Orbit columns K | Distinct points B | m=2 ceiling | m=3 ceiling | m=4 ceiling |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 10,624 | 2.334304e-17 | 8.268104e-14 | 2.196628e-10 |
| 256 | 42,496 | 3.734622e-16 | 5.290466e-12 | 5.620988e-8 |
| 600 | 99,600 | 2.051462e-15 | 6.810989e-11 | 1.695987e-6 |

Under the same count screen, the first `K` where the ceiling reaches one
is 13,247,128,046 for m=2, 1,469,216 for m=3 and 16,627 for m=4.
These are thresholds for a bound to stop being informative, not adequate
base sizes or predicted relation rates. The m=4 K=16,627 case would need
22,945,941,707 retained root states in the current materialized index
design, according to [S5_CAPACITY.md](../koblitz_n83_factor_base_sweep_20261008/S5_CAPACITY.md).
An implicit or segmented design needs its own correctness and memory
measurement. The dimension-12 closure's earlier exact four-summand
count/order ratio is only 0.000209791390 even though it has 332,166
points. Its count advantage over the 4,054-point source does not measure
solver speed.

The primary sweep executed 1,728 fixed-public-point two-summand probes
over its original panel with zero relations; these are dependent,
fixed-fixture diagnostics. Its K=64 m=2 cold preflight and one completed
trial, and its capped m=3 attempt, are recorded in the
[sweep result](../koblitz_n83_factor_base_sweep_20261008/RESULT.md).
There is no comparable natural m=4 relation-and-rank panel for these
stored primary bases. The 81-bit subgroup also requires full-width row
coefficients and rank arithmetic; the source audit lists the remaining
[solver transfer gates](../koblitz_n83_factor_base_sweep_20261008/SOLVER_GATES.md).

## The distinct a=1 four-summand evidence

The separate `EC1N83Ckb1a1haacb69ee9f37` curve has subgroup order
`8569786107849059`. Its [cold K1024 screen](../n83_cold_factor_base_20261016/RESULT.md)
reports 150.173–163.404 seconds for completed K1024/root22 arms across
three cloud allocations, with 169,984 usable points. The K2543 A/A control
exceeded its 10% spread gate in two allocations, and the selected base
remains open pending reserve-point confirmation. A prior
[five-point online panel](../n83_factor_base_panel_results_20261014/RESULT.md)
reported 14.4220794 ms mean target-online time for K3240 after 616.3 s
setup on another cloud run. Those clocks answer different cost questions
and cannot be transferred to the primary 81-bit curve. They do show the
practical distinction between paying setup cold for one target and reusing
an already prepared base across targets.

## How to choose a stronger prebuilt base

The next candidate must be a *base plus an oracle and index policy*.
Prebuilding has value only if the saved point set and its checked labels
can be imported and the remaining root index, relation generation, rank
and target query costs are charged. For equal point count, the multiset
expectation is identical across base shapes, but actual unique sum support,
solver search effort, row dependence and index memory can differ. A
selection by construction time or orbit count alone misses those terms.

A discriminating next comparison is the primary K=2,048 explicit-orbit
object against the near-size-matched dimension-12 closure on the *same*
public query corpus. Freeze both import hashes, the m=4 oracle and index
representation, query coefficients, search caps, host and memory limit,
then measure natural verified relations, new rank per query, solver cost,
base and index preparation, and query time. Preserve every capped result.
The checkable implementation and experiment gate is tracked in
[the N83 prebuilt-base problem record](../../problems/n83-prebuilt-base-quality-m2-m4.md).
The current four-summand backend and full-width adapter must first pass
the source and independent group-lift gates in `SOLVER_GATES.md`. A later
K=16,627 experiment should use a proved equivalent implicit/segmented
index and report its memory and construction costs; a source count alone
cannot admit it. For an end-to-end improvement, compare the complete cold
pipeline with the frozen current index-calculus baseline and same-point
rho under the repository's isolation and replay rules. Any speed claim
also needs the fastest applicable external implementation on the same
inputs and host with equal checked outputs.

No new admitted end-to-end number changes the canonical scoreboard,
progress timeline, curve graph or performance-gains graph in this review.
