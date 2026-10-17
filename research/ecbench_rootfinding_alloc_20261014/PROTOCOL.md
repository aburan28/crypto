# Reused buffers in the prime lane relation collector

Frozen before candidate timing. This is a local engineering experiment on the registered 20-bit prime instance. The baseline is the release `ecbench` binary saved at `/private/tmp/ecbench-rootfind-buffer-baseline-20261014` with SHA-256 `ab5764921ad616ae8177ca57c125ec8cd7ba09896f1e0644a6999d1fb33b7159`. The candidate changes only the synchronized lane collector's request, native left/right, and output vectors to reuse capacity across rounds.

## Question and controls

- Prediction: the same target, seed, lane width, and batch setting produce exactly the same ordered relation probes, coefficients, restarts, counts, ranks, and verified answer. Buffer growth should occur only when the required width exceeds capacity, instead of on every round.
- Primary correctness gate: the existing serial-versus-batched stream test at widths 4, 8, 16 and an exact baseline/candidate `ecbench` replay comparison on eight fixed public targets. Any mismatch invalidates the optimization.
- Timing gate: run the release baseline, candidate, then baseline again on one shared macOS arm64 host. Arms alternate within each session and include a rho A/A control. This session order exposes some drift but cannot remove host confounding. Report wall data as L0 descriptive only; an L2 host is required for a qualified speed claim.
- Counted cost gate: compare cold and phase operation counts and peak RSS. The sealed ecbench format omits native field-operation counts and does not price native allocation or lookup work, so its GAE cannot establish a whole-solve advantage.

## Limits and artifacts

Use eight public targets from the existing `prime_search` curve seed 59297, target seed 20261020, B48, `mitm-prime-batch:m=2`, width 8 serial and batch arms, matched `rho.negation`, three measured rounds and one warmup, at most 120 seconds per execution. Keep spec, plan, raw sessions, run logs, replay-all audit, source and binary hashes, comparison script/output, report, graph, and PDF under this directory. Stop on any correctness failure and preserve its raw output. Do not extrapolate to large curves or claim qualified wall improvement from this shared host.
