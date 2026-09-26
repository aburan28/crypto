# Development reference qualification, frozen before measurement

This is reference selection, not one of the three improvement rounds. The
hypothesis is that the restored `both`, `scaled` and `pairinv` pipelines can be
compared correctly on identical supplied public targets and that their measured
rankings may depend on the online versus cold accounting boundary. No 20% gain,
family-wide result, confirmation result or global-optimality claim is permitted.

## Frozen production panel and limits

Run the existing `tournament.py` with `--qualification --profile pilot`, seed
`2026092541`, development cells `17a1,19a0,23a0,23a1,31a0` and the reserved
unmeasured holdout cell `29a1`. Qualification prepares and executes only A/A,
smoke and development. It creates no selection or confirmation workloads.
Field degree is not subgroup bit length; exact subgroup orders, curve records,
generators and actual public points are in the frozen fixtures/manifests.

Use three process repetitions, one target per job, one pinned Linux amd64 CPU,
8 GiB address-space cap, 180 seconds per child, Rust 1.94.1 with the prepared
musl/ISA flags, and Valgrind 3.22.0. The job budget is 1400 profiled/native pairs;
the declared schedule contains 1290: 30 A/A, 315 smoke and 945 development.
The workflow has a 90-minute wall limit; interruption, timeout and OOM retain
artifacts and cannot produce a completed qualification. No automatic retry of a
failed measured job is allowed. OS caches are uncontrolled; cold means the whole
worker process, not a machine reboot. Hardware and build provenance are frozen.

All IC sources use summands=3, pair_table, tiny_gauss, batch_trials=1 and
max_trials=65536. Factor-base policy is the declared comparison variable: actual
support, usable point count, folded columns and source-specific construction
remain in each canonical candidate. Only the source policies already present
in the archived pipelines are measured. Sources are derived from hash-checked
archives with the reviewed instrumentation and allocator correctness patch;
the original archived bytes are unchanged. Independent row-arithmetic release
controls run before qualification. Full scalar-field elimination remains the
archived algorithm; its zero-prefix arithmetic is an optimization opportunity,
not evidence that a different kernel has already been measured.

Rho requested widths 1, 2, 4, 8, 16 and 32 run from each distinct source snapshot: eighteen
reference configurations. Preserve the actual effective widths on every cell;
clipped widths are not independent algorithms or independent target samples.
The intermediate widths prevent the old 1/8/32 screen from skipping a better
small-width reference when larger requests collapse to the same effective width.
Every reference is interleaved with IC in the same case/repetition blocks and
randomized order. Reusable rho arithmetic/Frobenius preparation is excluded from
its online interval. The public point, resources and target-dependent work are
matched. Current `icx` remains outside this comparison because its entry point
constructs target scalars internally and lacks this admitted public-point stage
adapter; the scope is the strongest *qualified archived* implementation here.

## Accounting, correctness and decision

The primary reported metric is verified single-target online native wall time,
after reusable preparation through scalar replay. Complete cold native time and
complete user-space guest instructions are separate columns/boundaries. Fixture
construction is outside measured jobs. All unsuccessful target-dependent attempts
remain charged. Exclusive phases, ordinary-query yield, novel rank, final matrix
work, target descent, certificates and raw failure records remain in the existing
per-run evidence. Unknown costs stay null. Instruction normalization is
`S = Ir / sqrt(r)`; the K-instruction full-rank-collector floor is the existing
weak implementation-specific bound, not a generic lower bound for all IC.

A/A must pass the existing instruction gate. Every declared smoke and development
job is retained; any failed smoke or incomplete development workload disqualifies
that reference. A complete reference needs independently verified group relations,
rank, factor logs, nonempty IC descent, scalar replay, phase closure and native/
profile agreement. Rho has its own identity and never fictitious IC stages.

Selection uses development data only. First take the median of the three process
repetitions per target; use an equal-cell geometric mean of paired target ratios.
Keep separate IC and rho leaders for complete instructions and online time.
Native cold ratio, then lexical alias, breaks exact ties. All other references,
cell tradeoffs and failures remain in the table. Do not manufacture a pipeline
from the cheapest phases or substitute one metric's leader for another's.
The 95% paired bootstrap intervals are descriptive development uncertainty;
selection is not a held-out significance claim. Before improvement rounds, lock
the appropriate complete reference sources/settings and the familywise 95%
confirmation rule, retain both leaders when different, and exclude this run's
points along with prior campaign points from fresh confirmation targets.

Success here means an auditable complete comparison and explicit selected
references or an explicit incomplete result. It does not enable promotion. The
three-round goal still requires at least 60 fresh confirmation targets, the
20% cold gates, per-cell limits, confirmation and replay, and merged evidence.

CI also runs a smaller wiring control with two configurations of one source on
n13a0. That control tests the qualification path and is not the five-cell reference
selection above. Its scope and actual inputs are separately sealed.
