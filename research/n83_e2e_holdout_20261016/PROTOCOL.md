# N83 one-target cold factor-base holdout fixture

Five new public points are frozen here before the exposed K-grid screen
`38940108434` reports any timing. The first 14-value job `38938203526` was
cancelled while queued, before a producer started; the replacement grid
brackets the exact K1102 root-table allocation boundary. The executed grid
commit is `551b9dfc632dce7a5ab8218b1b1dff6ae4b4d785`. The generation law is
`1 + (little_endian(SHA256(domain || u32_le(index))) mod (r-1))`, with domain
`n83-e2e-cold-holdout-five-20261016-v1`, subgroup order
`r = 8569786107849059`, and indices 0 through 4. The generator checks the
prior public-point manifests, subgroup membership and scalar replay. The IC
producer receives only one public point per fresh process; the answer sidecar
is reserved for the independent verifier.

The point file SHA-256 is
`722bc60cd010e533cd715a9bfb321fba52d950352c99684cd422b23af7f93d26`.
The answer sidecar SHA-256 is
`b4c2ae858770b7c8842703558984326bc0002b0f1deb3c898d5ac30562de8b96`.
The fixture and domain are immutable once this branch is pushed. All five
points are used; no point is replaced after its timing is known.

If the exposed grid passes every replay, resource and K1152 A/A gate, the
holdout candidates are K1152/root11 dynamic plus the three lowest cold-time
Ks **other than K1152**. If the grid gate fails, retain its rows and do not
pick finalists from its descriptive timings. On each holdout point, run a
fresh K1152 process first and last, with the three other finalists in rotated
order; alternate rho before and after IC across points. Require the same
source, host, base/index identity, scalar replay and resource gates as the
screen. Each fresh IC child has a 480-second cap and each rho child has a
120-second cap; the allocation must expose at least four CPUs, 15.5 GiB RAM
and 2 GiB free disk. Retain each target-dependent online interval separately
from its full cold interval.
The host receipt includes visible CPU core, package and sibling topology,
online NUMA nodes, effective cpuset and `lscpu` output. Every IC and rho arm
also retains before/after `/proc/stat` CPU-steal and interrupt counters,
available CPU-frequency readings and CPU pressure. These fields describe
the observed cloud envelope; the existing K1152 A/A and resource gates
remain the timing acceptance rule. Missing required topology or counter
receipts makes the independent replay incomplete.

Before selecting challengers, `materialize_selection.py` checks the exact
GitHub artifact SHA-256, ZIP integrity, frozen run ID and executed commit,
then replays the complete grid locally and requires byte-identical analysis.
It seals the raw analysis and the three challenger Ks in versioned files.
The runner checks those files and their digests before opening any holdout
point.

Rank the four pipelines by their mean full cold one-target time over the five
points, using the mean of K1152 A/B on each point. Report all per-point
differences and a paired two-sided 95% Student-t interval
(`t_(0.975,4) = 2.7764451051977987`). Select a winner within this finite
grid only if its lower paired bound is positive against each other pipeline,
and its advantage exceeds half that point's K1152 A/A difference on all five
points. Otherwise retain a null selection and all complete and failed rows.
The finite grid does not exclude other base geometries, index layouts or
construction algorithms; those need separate candidates and comparisons.
