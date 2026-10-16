# N83 cross-panel final factor-base tournament

This protocol freezes a fresh, same-host comparison after the direct,
root22/root21, root4 Bloom, and root11 high-K screens. It selects the best
measured complete pipeline within the current deterministic Frobenius-orbit
construction and the declared memory envelope. The compared factor-base
point sets, index layouts, and Bloom policies remain distinct candidate
identities.

## Candidate admission before the final measurement

The source runs are fixed as direct six-arm `38768544528`, root22/root21
challenger `38770004923`, root4 Bloom `38772262240`, and root11 high-K
`38773883285` in `aburan28/crypto`. Do not replace a run after seeing its
outcome; a failed run remains a failed source panel.

The two standing controls are K3300/root11/wide3 and K3596/root4/Bloom-off.
For each of the four source screens, independently replay its uploaded raw
artifact. Admit only arms with two complete source-pinned process receipts,
full factor-base and target scalar replay, correct phase accounting, zero
new CPU-throttled periods and OOM kills, and observed peak RSS at most that
screen's frozen guard. From each screen, add its eligible lowest arithmetic
mean target-online pipeline, regardless of whether that screen's 99% interval
selected it. Its within-screen status remains recorded. Resolve equal means
by the pipeline label's bytewise order. Deduplicate exact candidate
configurations in standing-control, direct, legacy, root4 Bloom, root11
high-K order. A screen with no eligible pipeline contributes no candidate,
and its failure remains in the selection record. Preserve the raw ZIP and
analysis SHA-256 for every admitted source. A source screen still in progress
does not satisfy this admission step.

The resulting list has two to six candidates. Freeze it with each complete
IC1 manifest and base digest before launching the final allocation. The
producer source is SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`.
For every arm use the exact source flags of its source screen, symmetric S3,
sparse canonical selector, rank seed 29, four rank threads, and one target
thread. Each process begins with an empty reusable base/index/log state.

## Fresh paired comparison and decision

`generate.py` fixes 4,000 public N83 subgroup points from domain
`n83-final-tournament-fresh-4000-20261015-v1`. It pins and excludes all
previous N83 panel and reserve points listed in that file. The independent
scalar sidecar remains outside every producer input and timing interval.
The public-point SHA-256 is
`cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424`;
the sidecar SHA-256 is
`1dbd4250796de4434badc2458089d7ad9fe63cc59d7b9d3f6825ad431476a13b`.
Use one cloud allocation with at least four visible CPUs, 14.5 GiB observed
peak RSS guard, at least 256 MiB of host memory above that guard, and at
least 2 GiB free disk. Every child has a 1,500-second cap; the full job has
a 330-minute cap. Record exposed CPU topology, affinity, cgroup CPU and
memory limits, throttling, steal time, calibration noise, host memory,
source/binary hashes, exact base and root-key counts, phase costs, all
attempts, stderr, and failure status. This is a cloud-instance exploratory
many-point selection panel.

Run candidates in frozen list order, then reverse order, each in a fresh
process on the same public points. For each target and candidate, average
its two online intervals only after both pass independent scalar and
four-point relation replay. Form ordered pairwise differences on those
per-target averages. Bootstrap 80 contiguous blocks of 50 targets with
10,000 resamples at seed `20261015960+i*16+j` for candidate positions
`i,j`. Select the lowest-mean candidate only if the upper endpoint of its
99% interval against every other candidate is below zero. Otherwise record
`BOUNDED_UNKNOWN` and retain all pairwise intervals and arm receipts. An
ineligible or incomplete arm also makes the overall choice
`BOUNDED_UNKNOWN`; a complete subset may be analyzed descriptively.

The selected pipeline advances to at least five new, previously unseen,
same-point single-target IC/rho pairs under the Tier 1 cloud-instance
contract, alternating arm order and retaining all calibration and noise
receipts. Report each target's verified online interval and same-point
rho/IC ratio; preparation stays separate. A Tier 1 result is scoped to its
frozen instance and resource envelope.
