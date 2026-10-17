# Full debug traits replay asserts in binary-subfield trace

- **Status:** open
- **Kind:** unknown cause / validation gap
- **Affected work:** [new volcano catalogue validation](../research/full_volcano_classes_20261017/README.md)
- **Next action owner:** unassigned, curve-traits maintainer
- **Last updated:** 2026-10-17

## Problem

The full debug `curve_traits build` does not establish a fresh replay of the
catalogue: existing binary-subfield code asserts that an absolute trace lies
in GF(2), at `src/cryptanalysis/curve_traits/subfield.rs:75`.

## Evidence and uncertainty

The native debug binary was built from checkout `fd830abfe` plus the metadata
import argument change. The failed output and deliberate stop are retained
under `research/full_volcano_classes_20261017/evidence/curve-traits-*`.
The run was stopped after the assertion; the queued full check was not run.
The offending binary model is not isolated, and a clean-base reproduction
has not been run. No root arithmetic or assertion was changed for this task.
The newly added six prime-field models do not use that binary-field path.
This evidence does not prove the cause or establish a release-profile failure.

The existing 2,296 trait records were preserved exactly. Only six new prime
model records were computed with the current native module. Other input hashes
were checked unchanged and the registry digest was refreshed. See the native
`refresh_traits.rs` and `evidence/traits-incremental.json`. This incremental
refresh is not a passing full-catalogue replay. The three uploaded graph
censuses and all-map replays passed independently of this catalogue operation.

## What would unblock it

- [ ] Isolate the exact binary model and reproduce its trace failure, including
  its modulus and validity checks, on this head and the unchanged base.
- [ ] Explain or fix the assertion without disabling it or weakening validation,
  then complete full traits generation and byte-identical freshness replay.

These checks establish catalogue replay only; graph-wide ECDLP/IC measurements
and independent isogeny implementation replay remain separate requirements.

## Next discriminating step

Run a native per-record debug diagnostic on the existing binary records,
preserving the first failing model and its field certificate. Stop on the first
assertion, compare the same model on the base, and route a correction to the
responsible input or arithmetic implementation. A successful release-only run
does not resolve an unexplained debug assertion. Recheck on that correction.

## Resolution and history

Unresolved. No assertion was disabled, and the full debug gate is not reported
as passing.
