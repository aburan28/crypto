# F6-IC T7 profile: Macaulay construction dominates

The preregistered stage diagnostic in [PROTOCOL.md](PROTOCOL.md) completed
on October 6, 2026. All 12 fresh-process runs exited zero, recovered and
replayed the same scalar for each public point, and had five exclusive
online phases summing exactly to online wall. T1 is `(73407,129763)` with
scalar 4785 and workload `146a1e9ee3c8`; T7 is `(98625,98119)` with
scalar 2391 and workload `ced1677f0976`. Both use the frozen n17 curve,
62 usable base points, 29 folded columns, one prepared log state, and
one Rayon thread. The same binary
`faa1d65d394e90b153ea4c1858a122c14fab59cc8744ab602b4e79b1d570171f`
ran inherited F4, original F6-IC, and attempt-local pair F6-IC.
[FREEZE.tsv](FREEZE.tsv), `candidates/`, [measurements.jsonl](measurements.jsonl),
and `runs/` retain exact IC1 candidate/workload/run IDs, source and input
hashes, complete raw outputs, stderr, exit status, and timestamps.

The timers below are **nested diagnostics inside target PDP**, not a new
set of exclusive online phases. Macaulay is the existing F4 profile's
build, elimination, and readback sum. Oracle is time inside F6's exact
node callback. Each run had 11 target PDP attempts; F4 made 1,635 matrix
calls, while each F6 arm made 691.

| T7 arm | Rep | Online ms | Target PDP ms | Macaulay build ms | Macaulay total ms (% PDP) | F6 oracle ms (% PDP) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Inherited F4 | 1 | 205.967 | 205.903 | 171.884 | 173.696 (84.36%) | 0 |
| Inherited F4 | 2 | 207.657 | 207.593 | 172.534 | 174.308 (83.97%) | 0 |
| Original F6-IC | 1 | 157.573 | 157.508 | 128.996 | 130.900 (83.11%) | 2.474 (1.57%) |
| Original F6-IC | 2 | 153.409 | 153.340 | 126.719 | 128.653 (83.90%) | 2.393 (1.56%) |
| Attempt-local pair | 1 | 171.230 | 171.162 | 137.302 | 139.372 (81.43%) | 3.999 (2.34%) |
| Attempt-local pair | 2 | 152.263 | 152.200 | 123.901 | 125.687 (82.58%) | 3.514 (2.31%) |

F6's matrix construction took 98.55% and 98.50% of its measured
Macaulay stage on T7; elimination was only 1.904 and 1.934 ms, and
readback rounded to zero in this path. The F6 oracle was called 1,013
times per T7 run. T1 independently showed the same direction: original
F6's Macaulay total was 2.478 and 2.956 ms of 3.149 and 3.658 ms target
PDP, while oracle time was only 0.043 ms in both runs. Full T1 values
and all counters remain in the 12 [measurement rows](measurements.jsonl).

For this fixed T7 search path, even eliminating every nanosecond timed
inside the original F6 oracle would save only about 1.6% of its online
wall interval. This is an inference from the nested timers, conditional
on leaving the matrix work and branch decisions unchanged. A better
oracle that changes which branches are explored could save more; the
pair-index variants kept the same 691 reductions and did not do that.

**Decision:** The registered Macaulay-majority condition passes in both
T7 F6 repetitions. The oracle-majority condition fails. The next
optimization target is the Macaulay **build**, especially row generation
and packing; the existing `KIC_F4_BUILD_SUBPROFILE` diagnostic can split
those further in a separately frozen profile. Another point-addition
microoptimization is unlikely to meet the 2× complete-call goal on this
frozen n17 path. The Mac host is unisolated, so these wall fractions are
exploratory, and instrumentation overhead is included in the interval.
No controlled CPU speedup, n83 transfer, ordinary n83 relation,
IC-versus-rho result, or general 2× F6 bound is claimed.

The focused native pair-closure test and release worker build passed.
Run `sh derive.sh` to reproduce the row and phase checks; the machine
receipt is [DERIVATION_CHECK.json](DERIVATION_CHECK.json).
