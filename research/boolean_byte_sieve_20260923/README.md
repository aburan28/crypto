# Projected syndromes and compiled Gray schedules

This standalone study continues PR #677 on generated public Boolean systems.
It implements the proposed necessary-condition filter, retains negative discovery
steps, and separately tests compiled full-word schedules. Production solver,
full index-calculus, calibrated-operation and rho costs remain unmeasured.

The primary is complete: every universal dramatic and incremental gate is
**REJECTED**. The compiled 64-point kernel and fixed dispatcher each pass 3/18
dramatic subgroup comparisons; these do not qualify for confirmation. All 192
systems complete with verified results. See [CONCLUSION.md](CONCLUSION.md) for
the complete costs, failures, alias diagnostics and validation history.

## Exact representations

The original system has at most 32 equations. Each equation occupies one bit
of a u32 syndrome. The first filter projects onto bits 0–7. Nonzero projection
rejects a point exactly; a zero projection never authorizes acceptance.

Six low Boolean variables form a 64-point block. High variables follow reflected
Gray order. Fixed low quadratic offsets and affine derivative images are compiled
inside each fresh solve. Every parity cancellation and changing coefficient is
preserved. The universal low-coordinate masks contain no fixture-dependent data.

The single-stage control maintains all six full-word linear coefficients and
checks the full syndrome for each survivor. The two-stage method instead packs
the six coefficients for equations 8–15 into one u64. Its changing coefficients
require one packed XOR per high transition. Survivors of that second filter are
evaluated against the complete original syndrome. Checks are exact even when
either projected subsystem is identically zero or highly dependent. A presumed
1/256 survival fraction is not used as evidence.

Scalar byte, NEON/SSE2 byte and bit-plane representations share the same policy.
The bit-plane method uses native NEON on AArch64 and portable u64 operations
elsewhere. Traced test modes compare complete diagnostics; final quiet timing
arms report null traces and still check every returned model independently.

For high step t=4T+r, the two additional low variables are visited in order
`(r XOR (r>>1)) XOR (2*(T mod 2))`. This preserves the retained four-low-variable
enumerator's first complete model exactly. Every lane computed in a block is
charged, including lanes after an earlier successful point. Partial points,
second-stage checks, full checks and small-domain fallback work are separate.

## Compiled schedules and controls

The lowest four high-coordinate transitions repeat the same Gray flip pattern.
Literal schedules compile those transitions with constant coordinate indices,
retaining general updates at the group boundaries. Each individual 16- or 64-point
block checks its cap before advancing. The compiled and looped implementations
must match models, work, and instrumented traces at every tested cap.

Full-u32 controls cover 16-point and 64-point domains with enumeration checksums
removed. The same compiled schedule is applied to full-word controls, so batching,
instrumentation and schedule compilation are not attributed to byte projection.
Both standard and compiled partial scans also run in the retained 16-variable
leaf search policy, preserving its prefix and original-label recovery. Quiet leaf
arms remove enumeration checksums but retain the prefix's diagnostics.

The dispatcher is fixed before holdout timing: compiled 16-point scanning through
n20, compiled 64-point scanning above n20. It is compared with the prior reference,
not claimed to beat the pointwise minimum of its own components. Exact per-input
model and work equivalence with the selected component is an evidence gate.

## Discovery and primary protocol

Seven immutable discovery probes use only n16/n24, seeds 17/937, three families
and three repetitions. They preserve the initial projection, hot-loop changes,
deferred coefficient updates, packed secondary filtering, bit planes, compiled
64-point schedules and compiled 16-point schedules. Discovery timing alone
does not qualify a method for promotion.

Equivalent wrappers showed timing differences under cyclic discovery orders.
The primary therefore uses four deterministic random permutations per fixture,
each followed by its exact reverse. This reduces fixed-neighbor bias; it is not
a claim of instruction-cache-cold execution or complete position balance.
Code-layout effects remain part of the recorded compiled program's timing and
are not described as new mathematical work.

The primary retains all 168 prior fixtures and adds 24 unused holdouts at
n12/16/20/24. All 32 predecessor methods remain. Forty-nine total methods and
eight repetitions give 75,264 observations. The fixed 39-method reference includes
quiet, full-word, scalar-projection and single-stage controls. A candidate in that
roster removes only itself. New compiled and projection treatments are evaluated
separately; comparisons among their components do not redefine the fixed boundary.

Every candidate must have a paired 95% lower bound above 2.0 in all eighteen
n16/20/24, family and regression/holdout groups, with all 192 cells verified.
Separate >1.0 incremental and >1.05 matched-control thresholds cannot replace
the dramatic gate. Any eligible treatment requires unchanged timed-source
confirmation on unused holdouts. No tuning follows holdout timing.

## Accounting and validation

Cold context totals include encoding, all coefficient images and schedules,
unsuccessful filters, full checks, recovery, destruction and result validation.
`partial_ns` is an exclusive subset. `full_update_words` counts u32 low-linear
coefficient updates; `secondary_update_words` counts packed-u64 updates. These
are distinct work counters, not a calibrated operation unit. Fixture generation,
reference-status preparation and formatting are outside arm timers but inside
process receipts. Candidate-specific peak allocation remains unmeasured.

The Rust suite checks byte images and every zero lane, full 32-bit equations,
pointwise coefficient transport, the exact legacy ordering, false positives
requiring full checks, deferred updates, quiet/traced policies, every small cap,
compiled schedules and recovered leaf models. Prior algebra tests are retained.
An independent Python verifier checks raw receipts, retained controls, semantic
work, model equality, quiet trace nulls, ordering and censored costs. Tests and
CI are producer validation, not independent external review.

Use a fresh output directory:

```sh
python3 research/boolean_byte_sieve_20260923/run.py --out /tmp/boolean-byte-replay
```

Executed source snapshots, protocols, raw samples and manifests are immutable.
The general generated-system driver does not invoke the WDSat/full-curve suite
and cannot establish an index-calculus or rho crossover.
