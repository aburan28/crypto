# Frozen order-229 different-j distinct-target comparison

**Status:** protocol fixed before implementation and before any batch
measurement. The hypothesis is that using the full rational image of an
explicit degree-three isogeny can change pair-table relation support or
complete solve cost after setup is reused. A support change alone is a
stage result. The end-to-end question is whether the full-image arm's
cold, same-target total falls below the cheaper direct source control.

## Frozen instance and boundary

Use only the public synthetic instance in the
[screen](../different_j_order_screen_20261015/RESULTS.md) and
[independent full-fiber audit](../different_j_order_screen_20261015/FIBER_REPORT.md):
`F_211[s]/(s²−2)`, source `E:y²=x³−21x+37`, destination
`E':y²=x³−81x−243`, kernel `{O,(3,1),(3,−1)}`, `r=229`, cofactor
`h=195`, source generator `G=((170,0),(0,109))`, and image generator
`G'=((76,0),(0,69))`. The separate audit checked every 3-point fiber
and all 52,441 subgroup-pair identities. The models have different
`j=77,160`; neither has the exceptional `j=0,1728` residue. The
generic negation-orbit collision floor in the stipulated group-addition
equivalent unit is `S_floor=sqrt(pi/(2·2))=sqrt(pi/4)`. The reference
values are the **measured** complete same-target direct source BSGS
and negation-rho controls in this panel, reported with intervals.

## Inputs, controls and execution

Enumerate all nonzero scalars `1..228` only in the target-freeze
process. Hash the public source point tuple with domain
`different-j-order-batch-target-v1:13`, sort by digest then tuple,
and take the first 64 distinct source/image point pairs. Commit a
solver-public JSON file with coordinates and curve identities but no
scalars or answer-bearing fields. Keep the 64 scalar correspondences
in a separate checker-only JSON file; the solver must never open it.
The checker validates target order, map, distinctness and every
recovered scalar from independent field and group formulas.

Select eight source/projected columns by scanning lines `x=t(a+s)`
in ascending `a=0..210`, `t=1..210`, with `y` ordered by public field
packing, retaining the first eight distinct signed subgroup orbits.
Use the same rule on the destination for the direct base. For the full
image, enumerate all 44,655 source points, form and verify the 14,885
distinct rational images with three preimages each, sort images by
public point key, cofactor-project, and retain the first eight distinct
signed subgroup orbits. Commit the selected bases and source
preimage witnesses before measurement. Keep candidate scans and map
work charged during the measured run, including duplicates and rejects.

Run seven arms on the exact same 64 targets: source projected
relation; projected isogeny-image relation; destination direct
relation; destination full-image relation; reusable direct source
BSGS; direct source negation rho; and rho A/A. Use three arm-rotated
rounds, one process per arm/round, with the same committed round/target
random stream in all relation arms. Include prefixes `N=1,4,16,64`.
Fix eight factor-base columns, 1,220 relation queries per target,
rank nine for complete recovery, and at most 32 rho restarts with
4,580 steps each. The BSGS table has 16 baby entries and up to 15
giant lookups. Per process: 120-second wall ceiling and 512 MiB RSS
ceiling. Preserve every failed attempt, capped rank, timeout, exit,
stdout, stderr, seed, source/binary/input hash, and memory receipt.
Do not replace capped rows or lower the target count after inspection.

## Verification and accounting

Require 21/21 arm-round processes and 1,344/1,344 scalar outputs
(`7·3·64`) to verify on both curves. Reconstruct every relation
witness independently, including the signed base pair, map and
linear row; require each relation arm to reach rank nine for every
target. Check full-image base witnesses against the frozen fiber map.
An arm that caps or fails is retained and marked incomplete, not
silently filled from its known-answer fixture.

Charge field construction, point enumeration, map evaluation,
factor-base and pair-table setup, all queries and failed attempts,
linear elimination, DLP recovery and scalar verification, all process
startup/input/output, memory, and the exact target count reusing
setup. Record per-phase native field multiplication/inversion and
group-operation counters, whole-child CPU/wall, and peak RSS. Report
the cold prefix totals in `S=operations/sqrt(229)` against both direct
controls only where every primitive is priced in a measured conversion
with a stated uncertainty; unpriced work makes that `S` bounded. Also
report direct whole-child CPU ratios as shared-host L0 descriptions,
separate from portable counted cost. Calibrate any CPU-to-group-addition
conversion on the same field and host before interpreting results; it
does not turn L0 wall time into an admitted speedup. Equal memory
ceilings are required for all seven arms. No scaling or exponent claim
follows from this one subgroup size. Isolated interleaving and an
independent host remain gates for qualified performance conclusions.

Commit runner, independent checker, exact scripts, source and binary
SHA-256, public/checker input SHA-256 and ceilings before measurement.
The source code and hashes are a separate freeze record following this
protocol. Run a complete release-library suite and touched-example
checks before pushing. Preserve this protocol and any failed run
version verbatim.
