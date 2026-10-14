# Conductor-gap Callgrind phase-dump volume

- **Status:** investigating
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

- [ ] A regression demonstrates that internal dump suppression preserves
  exclusive clocks and unconditional external markers, and default behavior.
- [ ] A newly pinned native executable emits a bounded number of parts for a
  full IC solve; its scalar and every deterministic phase counter match the
  archived record.
- [ ] All 100 frozen jobs execute under the amended protocol, and raw-part
  re-derivation plus native verification admits the complete panel.

The first condition is implemented and passes standalone module tests. The
release rebuild and full census remain pending. These conditions establish
only the measured instruction census; structural hardness remains unresolved.

## Next discriminating step

Run the same complete frozen panel with `ECBENCH_CALLGRIND_PHASE_DUMPS=0`,
the v2 pinned binary, unchanged 1800-second child watchdogs, and the same
scalar/counter controls. A correct bounded-output solve permits completion of
the panel; a missing marker, counter mismatch or watchdog exit remains a
failure and blocks inference. Recheck when the v2 release build completes.

## Resolution and history

Unresolved. Both interrupted attempts are retained losslessly and excluded
from primary estimates. No instruction result has yet been admitted.
