# Hosted one-shot screen dispatch

The [input freeze](INPUT_FREEZE.md) and source lock are merged. This workflow
exists to execute the [preregistered protocol](PROTOCOL.md) once on the frozen
five-block public-Q set. It does not alter any source or input pin. The
on-pull-request `validate` job independently replays the 5,120 Q and
source-only controls. Only a manual dispatch of the merged workflow on `main`
starts a charged arm.

The `measure` job materializes the pinned 20-file v2 source tree and builds
fresh release generator, compact, and rho executables **before** reserving a
core. It records build logs, compiler and binary hashes, and the synthetic
source materialization. It then reserves one Linux x86-64 CPU and its SMT
sibling, uses a one-second contention sample and fixed 0.10 other-CPU-second
gate, and runs the exact 45-arm schedule once. It may retry only a busy-machine
preflight before the run directory exists, preserving each refusal log. Once
any measured arm begins, a failure or timeout is retained as censored evidence
and does not trigger another arm run.

The raw artifact contains all child commands, output, process CPU/wall/RSS
receipts, source/input hashes, isolation record, and preflight logs. A separate
hosted job downloads that artifact and independently replays every full-rank
trace, four-point witness, and recovered scalar. Relocation affects only
absolute file paths; original command and build-binary paths must still agree.
The replay job applies the preregistered classification if and only if all
children and timing gates pass. It uploads its own artifact even on a replay
error. Artifacts are retained for 90 days; the outcome PR must seal the raw
member bytes or a durable, audited manifest before that deadline.

After this workflow PR merges and its exact-head checks pass, dispatch with:

```sh
gh workflow run ecc2k130-base-window-hosted-screen.yml --ref main
```

Record the run ID, main SHA, raw and replay artifact digests, all failures,
verification result, charged per-arm CPU, fixed-window ratios, A/A intervals,
and decision in a separate outcome PR. A favorable hindsight minimum remains
an optimistic diagnostic, never a deployed selector or ECC2K-130 speedup.
