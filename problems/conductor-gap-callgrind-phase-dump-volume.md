# Conductor-gap Callgrind phase-dump volume

- **Status:** resolved for the v2 instruction census
- **Kind:** profiling implementation obstruction
- **Affected work:** [instruction census](../research/ecbench_conductor_gap_ir_20261014/PROTOCOL.md), [v2 mechanism](../research/ecbench_conductor_gap_ir_20261014/PROTOCOL-v2.md)
- **Next action owner:** native ecbench measurement workstream
- **Last updated:** 2026-10-14

## Problem

Whole-solve instruction profiling of the existing prime-field IC pipeline
produces a profile dump/reset on every phase transition, including repeated
relation-loop transitions, so the original census cannot finish within its
intended bounded output envelope.

## Evidence and uncertainty

- **Observed:** the committed v1 attempt produced 188,865 raw files in its
  first job directory before interruption, including at least 80,938 numbered
  profile parts observed during execution. It emitted no completed-job log
  message. Raw output is retained in `attempt-phase-dumps/` in the study.
  `ic_measurement::Trace::boundary` unconditionally emitted a dump/reset.
- **Derived or proved:** no mathematical obstruction is established.
- **Suspected explanation:** repeated query/PDP/check transitions drive the
  volume; source inspection identifies that emission point.
- **Missing capability:** whole-solve profiling with retained collection and
  outer markers, without per-transition output files.
- **Limits:** interruption does not establish a solver failure, conductor
  effect, absence of an attack, or intrinsic hardness.

## What would unblock it

- [x] A regression demonstrates that internal dump suppression preserves
  exclusive clocks and unconditional external markers, and default behavior.
- [x] A newly pinned native executable emits a bounded number of parts for a
  full IC solve; its scalar and every deterministic phase counter match the
  archived record.
- [x] All 100 frozen jobs execute under the amended protocol, and raw-part
  re-derivation plus native verification admits the complete panel.

All seven module checks pass, including the isolated strict-session check.
The pinned release executable completed all 100 profiles. Every native Ir
receipt has exactly one solve part; every archived scalar, point and
deterministic phase counter matches. [Native analysis](../research/ecbench_conductor_gap_ir_20261014/analysis.json)
and the lossless raw archive record this bounded-output census. These conditions
establish only the measured instruction census; structural hardness remains unresolved.

## Next discriminating step

The complete frozen panel passed with `ECBENCH_CALLGRIND_PHASE_DUMPS=0`,
the v2 pinned binary, unchanged 1800-second child watchdogs and scalar/counter
controls. Recheck the clock/marker regressions and append a fresh bounded-output
census if this mechanism changes. Missing markers, counter mismatches and
watchdog exits must remain failures that block inference.

## Resolution and history

2026-10-14: resolved for the pinned v2 census after all 100 profiles passed
native validation. Both interrupted attempts are retained losslessly and
excluded from primary estimates. The frozen repeated instruction-lead criterion
was not met; both floor/crater 99% intervals include one. This resolves the
profiling obstruction, not the broader hardness question or missing algebraic oracle.
