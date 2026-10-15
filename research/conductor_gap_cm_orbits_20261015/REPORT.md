# Native conductor-gap controls and larger explicit floors

15 October 2026. This report is cumulative; sealed run directories and the
historical evidence archive are immutable. The 100-pair solve campaign and
all 30 replication sessions are complete. The larger-memory 96-bit retry is
separate work in progress; no unconstructed floor is counted as a result.

## Completed observations

The primary campaign verified all **500 attempts**: 300 orbit-folded IC and
200 independent same-point rho runs on the original 200 endpoint targets.
Floor/crater native-backend CPU is **0.979702**, with paired percentile 95%
interval **[0.866023, 1.111363]**. The query-count ratio is **0.878396
[0.764263, 1.009140]**. Ten thousand bootstrap resamples retain whole gap
pairs. The timing interval crosses one; this solver does not establish a
systematic crater/floor mean timing difference in the tested population.
Native-backend CPU includes orbit/table construction and native recovery;
CM-parameter derivation and independent replay are included separately in
the complete measured-process receipts. All primary/replication timings are
co-scheduled L0 diagnostics, with incomplete common operation prices.
Crater and floor shards ran mostly in separate queue waves. Host/wave
confounding is not removed by the paired bootstrap; these timing ratios
are descriptive rather than a controlled causal comparison.

All **265 replication attempts** also verified: 155 IC and 110 rho. The
same-crater CM/scalar scientific outputs agree exactly in **45** further
comparisons, as well as all **100** primary comparisons.

| Replication cohort | Gap clusters / paired targets | Floor/crater query ratio, 95% cluster interval | Floor/crater native CPU, 95% cluster interval |
|---|---:|---|---|
| Fresh targets | 10 / 20 | 0.988 [0.845, 1.149] | 1.059 [0.912, 1.221] |
| Isomorphic coordinates | 10 / 10 | 1.397 [0.882, 2.056] | 1.465 [0.974, 2.049] |
| D = −43 | 10 / 10 | 1.375 [0.998, 2.011] | 1.339 [1.006, 1.887] |
| Pair 46, five query-seed rounds | 1 / 5 | 0.977; descriptive | 0.951; descriptive |

The D = −43 CPU interval narrowly excludes one in this exploratory subset;
its extra CPU accompanies extra queries. These intervals are unadjusted for
multiple comparisons, and timing is L0. This is a follow-up lead requiring
more fresh targets and interleaved endpoint blocks, not an intrinsic-hardness conclusion. The 20 additional
floor vertices completed another 40 same-point IC/rho attempts. Per-attempt
records and whole-gap cluster resampling are in
[replication-metrics.json](results/replication-metrics.json).

### A generic recurrence changes the baseline

The canonical action satisfies λ² − λ + 5 = 0 modulo the prime group order.
Consequently P[k+2] = P[k+1] − [5]P[k], after one initial multiplication,
generates the same orbit on either endpoint using generic arithmetic. This
control needs no degree-5 floor endomorphism. The original comparison against
repeated full scalar multiplication overstates the advantage available for
this particular linked orbit workload.

The additional 12,000-trial campaign uses six arms over all 100 gaps and 20
shuffled repetitions. Mean thread CPU per 8,192-point construction is 5.265 ms
for the crater CM map, 6.643 ms for crater recurrence, 6.644 ms for floor
recurrence, and 41.494 ms for direct floor scalar multiplication. All 4,000
recurrence/direct comparisons agree, as do 2,000 crater map/scalar comparisons.
The backend compares complete vectors in memory, with supplementary FNV64
fingerprints. These are co-scheduled L0 observations, not isolated speed claims.

All 100 same-crater full-solve pairs have now completed. Their factor-base
SHA-256 hashes, query counts, retained relations, ranks and recovered targets
agree exactly. Native aggregation asserts that equality rather than inferring
it from similar times. The matrix may pin the target before all five unknowns
are determined; final verification checks [d]G = Q independently.

Another 12,000 trials ran pinned to CPU 0. The local monitor recorded no
contention, 0.34 other CPU seconds during 300.43 wall seconds, and 51 moved
threads. This is a local reservation receipt, without a physical-host or
frequency-isolation claim. Means are 5.464 ms for the CM map, 6.859 ms for
crater recurrence, 6.857 ms for floor recurrence, and 42.141 ms for direct
floor scalar multiplication. First observed model/arm passes and later
shuffled passes are reported separately; each arm was not restarted in its
own process. See [receipt](results/isolated-receipt.jsonl) and
[summary](results/isolated-summary.json).

Twenty further fresh-process batches produced 12,000 cold-process orbit
records, using the unchanged frozen backend arithmetic. Means are 5.104 ms
for CM, 6.658 ms for crater recurrence and 6.666 ms for floor recurrence.
These batches overlapped validation/construction at L0 and are not a causal
cold-versus-warm timing comparison. Each fresh process contains all 600
model/arm jobs; process-level CPU includes parameter loading and output,
while individual orbit CPU excludes prior CM-map derivation. The complete
single-target campaign charges that derivation in its cold process window.
The construction campaigns total 44,000 trials, including the initial 8,000.
[Cold summary](results/cold-summary.json), [fresh-process times](results/cold-passes/).

Raw data: [recurrence trials](results/recurrence-orbits.jsonl),
[native summary](results/recurrence-summary.json), and
[current campaign snapshot](results/campaign-progress.json). The snapshot
contains complete sessions and verifies each record-file hash against its seal.
The original checkpoint report and figures are preserved in
[reports/initial-checkpoint](reports/initial-checkpoint/).

### Independent structural witnesses

PARI point-order and Weil-pairing checks passed for 136 cases: the original
100 pairs, 15 newly constructed pairs, 20 extra floor vertices, and one small
D = −43 control. On quadratic twists, the crater has invariant factors
ℓ and qℓ, while the floor has a point of order ℓ² and is cyclic of order qℓ².
Here q is 5 or 11. The original ECDLP endpoints both retain the same prime
order N; this structural observation does not substitute a different subgroup
for the ECDLP comparison. The module structures imply ℓ+1 versus one rational
degree-ℓ kernel lines, using Frobenius on ℓ-torsion. Those counts are derived,
not an enumeration of every isogeny. Explicit witnesses and pairings are in
[torsion-structure.jsonl](results/torsion-structure.jsonl).

The 20 extra floors have distinct j invariants within their respective gaps
and independent certificates in
[extra-floor-certificates.jsonl](results/extra-floor-certificates.jsonl).
The 30-spec replication queue completed 265 verified attempts: fresh targets,
isomorphic coordinates, repeated query seeds, additional floor vertices and
ten D = −43 pairs. All requested attempts are retained.
The coordinate control transports G and Q by a degree-one isomorphism but
reruns the public prefix-seed selection rule in the new coordinates. It can
detect representation sensitivity; it does not isolate arithmetic cost while
holding the geometric seed points fixed. The second CM family also changes
the sampled gaps and is a transfer control rather than a one-variable test.

Independent PARI replay also checks recovered public targets and rejects a
scalar incremented by one. Its initial exporter failed on controls whose
target was generated rather than embedded in the curve spec. The corrected
exporter reads the actual workload target; the failed receipt is retained as
`target-replay-first.stderr`. The final replay independently verifies all
797 completed recoveries, including controls, prime-order/Hasse witnesses
and rejected scalar mutations. A separate native post-hoc replay rebuilt
and hashed every base, checked relations and ranks, and rejected corrupted
relation signs for all 473 IC records. The six initial fixture-transport
preflight failures remain preserved; they are not counted as solves.

### Original 96-bit candidate: an explicit prerequisite barrier

For the separate original candidate with ℓ = 129149576694539 and
p = 79228162511798403263235114731, the least extension degree for a nonzero
kernel point on the selected eigenspace is 2,226,716,839,561. The twist point
degree is 4,453,433,679,122, and the x-coordinate degree is again
2,226,716,839,561. Native PARI checks certify the scalar's multiplicative
orders. A dense 96-bit coefficient representation alone would need about
26.7 TB per field element. This explains why the present rational-kernel-point
constructor does not apply directly; it does not exclude other algorithms
using rational kernel polynomials. Its floor remains unconstructed. See
[prerequisite certificate](results/nonrational-kernel-prerequisite.json).

### Initial four-arm observations

All 8,000 construction-only trials finished: 100 gaps, 20 shuffled repetitions,
three scientific arms and a duplicate scalar A/A arm. The backend compared
entire crater point/column/coefficient vectors, in addition to recording
supplementary fingerprints. All 2,000 paired map/scalar comparisons agreed.
Each trial constructed 8,192 points from four seed columns and 2,048 levels.

| Arm | Runs | Mean thread CPU per construction |
|---|---:|---:|
| Crater, degree-5 CM map | 2,000 | 0.005104398 s |
| Same crater, scalar action | 2,000 | 0.040143001 s |
| Floor, scalar action | 2,000 | 0.040216660 s |
| Crater scalar A/A | 2,000 | 0.040501120 s |

These are descriptive L0 observations under concurrent work, not an isolated
wall speedup. Native operation counters retain unpriced modular exponentiation,
hash, allocation and elimination work. The construction-stage difference does
not establish a difference in total ECDLP hardness.

Raw data: [orbits.jsonl](results/orbits.jsonl), native-derived
[summary](results/orbit-summary.json), [conditions](results/CONDITIONS.md),
[pre-run hashes](results/orbit-pre-run.sha256), and [native source](native/folded.cpp).

The new square-root Vélu codomain implementation matched both coefficients of
all 100 previously constructed direct codomains. The reference check retained
every row in [sqrt-velu-reference.jsonl](prepared/sqrt-velu-reference.jsonl).
A separate independent PARI run proved primality, checked both endpoint
orders using Hasse-unique prime-order witnesses, and checked the prime-kernel
witness for all 100 pairs: [certificates](results/pair-certificates.jsonl).

The independent checker also certified 15 newly constructed floor models:
five D = −19 fields through 91 bits and ten D = −43 controls. These include
prime conductor gap 17,592,186,219,863 over
p = 1547425079971670554692473983. Exact coefficients, subgroup generators,
orders and discriminants are in [floor-certificates-v2.jsonl](results/floor-certificates-v2.jsonl).
Classification uses the ordinary CM order lattice and the class-number-one
maximal-order j invariant; this mathematical bridge is not claimed as a Lean
formalization. No explicit full isogeny point map or dual is certified here.

The first upper-ladder process was stopped
after the 91-bit result when combined memory reached 15,982,964,736 bytes of
a 17,179,869,184-byte cgroup limit. [Stop receipt](results/larger-construction-stop.json).
Its first output schema marks an absent direct reference as true; those values
are not validation receipts. The corrected version emits null and an explicit
`reference_supplied` field. All new endpoint claims use independent PARI checks.

The batched version matched all 100 reference codomains, including forced
small evaluation batches. It bounds the remainder-tree memory by evaluating
giants in blocks, trading memory for repeated polynomial work. Under a 6 GiB
address-space cap, the 95-bit case reached two of four blocks before its
1,800-second wall cap; both 96-bit cases aborted on FLINT allocation failures.
The launcher failure caused by missing GNU time is retained separately and
does not count as an algorithm attempt. The same near-limit rational-family
96-bit input reached 5,242,880 of 5,590,170 giant points under an 11 GiB
address-space cap before its 3,600-second wall cap (exit 124). A separate
precomputed-tree retry uses a 14 GiB address-space cap and a 3,600-second
wall cap, after all solver tables have exited. Its exact binary first matched
all 100 direct-reference codomains. This is distinct from the older non-rational-kernel
96-bit candidate. See the versioned [launchers](native/run_large_batched.sh),
[expanded retry](native/run_96_expanded.sh),
[precomputed retry](native/run_96_precomputed.sh) and their raw time/exit receipts.

## Native end-to-end harness

`ic.cm_folded_native_v1` runs the existing compiled native relation backend
through ecbench, with an executable-content hash required before launch.
Rust independently checks every factor-base scalar, every retained relation,
matrix rank, and recovered public point. No validation scalar is an input.
All reusable setup is included in supplementary cold cost. The five exclusive
online phases include independent replay and sum to the complete online clock.
Unpriced work stays explicit; zero charged GAE is not a complete zero cost.

The first small control failed all six jobs before solving because measured
children lost the supplied public fixture. That session is retained in
`sessions/control`. `CurveSpec::for_child` fixes measured execution, profiling
and replay. The next control verified all six jobs, including exact map/scalar
query-count agreement in both rounds: `sessions/control-fixture-fix`.

The primary campaign completed 300 IC attempts and 200 independent same-point
rho solves across the original 200 public regression fixtures. Its eight
shards used four co-scheduled workers at L0, with no shared target tables.
Source, executable and spec hashes were saved before launch in
[campaign-pre-run.sha256](results/campaign-pre-run.sha256). Do not interpret
a running shard, a censored target or an incomplete cost vector as a win.

## Requirement ledger

| Requested work | Current evidence/status |
|---|---|
| A: 100-pair repeated orbit comparison | 44,000 constructions, including recurrence, A/A, locally pinned repeats and fresh-process batches. Exact vectors agree. First/later pass and process boundaries are explicit; no physical-host isolation claim. |
| B: 200-target folded campaign plus controls | All 500 attempts verified: 300 IC and 200 rho. No complete-price IC/rho speed ratio is admitted. |
| C: fresh targets, repeated query seeds, coordinate/floor/family replication | All 265 attempts in 30 sessions verified. Ten second-family pairs and 20 additional floor vertices certified; 45 same-crater output comparisons agree exactly. |
| D: separate structural and generic effects | All 145 same-crater comparisons agree exactly. Generic recurrence removes most direct-scalar orbit cost; 136 torsion-module cases certified. No prime-field F4/SAT complexity claim. |
| E: native fast codomain construction and larger gaps | All 100 direct references match; 15 new floor models independently certified through 91 bits. 95/96-bit rows and the original non-rational-kernel 96-bit floor remain open. |
| PR, reports, required checks and canonical indexes | In progress. Passing controls and opening a PR do not complete the research ledger. |

## Next steps supported by these results

1. Repeat the D = −43 comparison on more fresh targets and gaps, with the
   generic degree-11 recurrence control. The ten-pair exploratory CPU signal
   accompanies query-count variation and needs independent confirmation.
2. Transport the four seed points as well as G and Q in the isomorphic-model
   control. The completed prefix-seed control changes geometric seeds, so it
   measures representation sensitivity rather than arithmetic alone.
3. Price the missing native work before publishing an IC/rho work ratio.
   Retain complete cold process costs and separate the five online phases;
   use earned reservation/isolation for any primary wall-time claim.
4. Finish and independently certify a 96-bit rational-family floor, retaining
   every cap and allocation failure. The original 96-bit candidate needs a
   different kernel representation; its extension-degree prerequisite is now
   explicit. Full isogeny point maps and dual-point-map certification remain
   separate work.
5. Implement an independently checked prime-field algebraic backend only
   when testing algebraic IC complexity. The completed MITM experiments do
   not measure Semaev/F4/SAT regularity or establish intrinsic hardness.

The full Rust library run passed 4,116 tests with 117 ignored; updated fixture
validation and CM-map tests also pass, and pinned Rust 1.98 strict Clippy
passes. Native ecbench audits check the sealed sessions. The browser includes
all 47 study sessions and keeps unpriced zero charged bounds out of its S
rankings. Per-session reports and 150-dpi plots are in [reports](reports/).

## Proof and provenance

The three evaluated polynomial identity certificates are
`IDC1hd41ecc586e59b223` (CM norm numerator),
`IDC1hc6bfc34b0ebbfa79` (paired second moment), and
`IDC1hd387bbbb1058c248` (paired third moment). Each has 64 issuer points and
64 additional domain-separated checker points. See [certificates](certificates/).

Lean 4.34.1 checks conditional orbit-equality, two-state recurrence and norm-bound lemmas with no
axioms. [Proof map](lean/README.md) lists the unformalized CM bridges explicitly.

The original six-page PDF and source evidence remain in
[historical-evidence-summary.zip](inputs/historical-evidence-summary.zip).
The initial source revision was e15e29d. Conductor coordination was attempted
but returned `not logged in`; no shared task lease is claimed. Work is on an
isolated feature branch. No remote scientific registry or paper is modified.
