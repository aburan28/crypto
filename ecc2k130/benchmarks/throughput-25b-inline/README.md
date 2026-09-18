# Selective polynomial inlining experiment

Engineering target declared before timing: median >25 billion complete scalar
updates/s on one RTX PRO 6000 Blackwell Server Edition. Reference is the
unmodified slot-unroll-2 preset, previously 17.2976 B/s; rerun it paired with
candidates. Count finished scalar updates, never field operations or progress
lines. Require 300 replayed reports and zero dropped reports per candidate.

The arithmetic is unchanged: 5.3125 products/update at batch 16. The prior
THROUGHPUT-29B.md pipe model puts the CLMAD ceiling near 26.6 B/s; this is a
hardware model, not an achieved rate or a generic-group advance. The generic
rho operation count and its ratio to the generic-group bound remain unchanged.
Inlining can remove calls and expose common subexpressions, but cannot lower
that arithmetic floor. Reject variants with incorrect replays, and retain
regressions. Report ratios to the paired baseline and the 25 B/s target.

PACKED_INLINE_POLY is a bitmask: 1 inlines single products, 2 inlines paired
products, 3 inlines both, and 0 preserves the baseline. Original checkout and
its uncommitted work remain untouched; this workspace contains a snapshot.

## Measured result

Target not reached. Three repetitions per variant, alternating order; third
repetition uses independent run-id 1937. Two verification runs per variant,
run-ids 0 and 1937, each replaying 300 reports against the scalar reference.
Verification uses 1,024 threads; timing uses the automatic 96,256 threads,
16 slots, 1,024 steps and 32 launches (50,465,865,728 scalar updates/sample).
Every accepted timing run passed the GPU-contention monitor. No GPU clock or
power settings were changed for this experiment.

| Variant (engineering) | Median B updates/s | / paired baseline | / 25 B target | Registers | Report replay |
|---|---:|---:|---:|---:|---|
| Baseline | 17.295877 | 1.000000 | 0.691835 | 110 | 600/600, zero dropped |
| Inline single products | 17.296764 | 1.000051 | 0.691871 | 110 | 600/600, zero dropped |
| Inline paired products | 17.940329 | 1.037260 | 0.717613 | 102 | 600/600, zero dropped |
| Inline both | 17.946263 | 1.037604 | 0.717851 | 98 | 600/600, zero dropped |

Inlining both improves the observed median by 3.760%; paired-only by 3.726%.
The tiny difference between those candidates does not establish that inlining
both is better. Single-only is effectively unchanged. The best observed rate
is 71.785% of the requested target, needing another 39.305% increase to reach
25 B/s. This is a walk-kernel engineering gain, not a full-DLP speedup or a
change to rho's generic-group operation bound. Products/update remain 5.3125;
the ratio to the algorithmic floor is unchanged.

All four kernels have zero stack and zero spills. Entry instruction counts
are 4,088 / 4,408 / 5,024 / 5,352, with 21 / 18 / 17 / 14 static call sites.
Those entry-only counts exclude callees and are code-size diagnostics, not
dynamic instruction or whole-walk operation counts. The no-spill reduction
from 110 to 98 registers does not by itself prove an occupancy improvement.

The default remains 0: these results apply to native CLMAD, this preset, and
this GPU. Globally enabling inlining for software arithmetic could enlarge
code substantially. Use PACKED_INLINE_POLY=2 or 3 with the frozen native preset.

## Reproduce and review

Run from ecc2k130/ on an idle GPU:

```sh
python3 benchmarks/throughput-25b-inline/build.py --output /tmp/ecc2k-inline-rerun
python3 benchmarks/throughput-25b-inline/run.py --output /tmp/ecc2k-inline-rerun
```

These scripts require a new output directory and preserve the frozen evidence.
Build arguments are frozen in build-config.json.
optimization.patch contains only this experiment's changes to Makefile and
include/packed131.h, excluding pre-existing edits. Apply from the repository
root with git apply after reviewing against any newer changes.

results.json contains raw outputs, completion validation, configuration,
source/binary hashes and paired samples. host-test.log records independent
arithmetic checks. build-*.log and registers.json record compiler evidence.
The earlier aborted attempt is preserved in attempt-1/: another GPU job
appeared during verification; it contributed no timing samples. The original
/home/ubuntu/crypto checkout has not been edited by this experiment.

The accepted rerun also recorded contention during baseline seed-0 verification;
that replay passed. No accepted timing sample had detected contention.

## PR scope and source provenance

The 17.946263 B/s measurement applies to the frozen tuning snapshot in
measured-source.tar.gz, not to current main with only this knob enabled.
That snapshot includes the existing pair-ILP, L2-persistence and slot-unroll
preset. This PR ports only the opt-in inlining control to main, leaving those
other tuning changes out of scope. The archive hash is in measured-source.json;
its Makefile, packed131.h and packedkernels.cuh match the hashes in results.json.
The build script rebuilds that exact snapshot. The current-main port receives
separate host arithmetic validation; it has no new throughput claim.
The optimization.patch file records the original snapshot change, not a patch
against current main. The Git diff is the authoritative current-main port.

Exploratory paired log-ratio t interval (three pairs, df=2) for inline-both:
95% CI [1.035305, 1.040346], with the
small-sample normality assumption. See paired-comparison.json; this is not
a new measurement or evidence of full-DLP improvement.
