# Frozen rotated m5/m6 semantic-gate evidence

Successful run source: preregistered #770 head `98c83a0` after the preserved
preflight-only failure in `failure_0/`. The corrected `FROZEN.json` SHA-256 is
`be1cd1db1556a4a36cb79040e576adcbb1b6829e3d34934876d4a427d90c3d8a`.
The #767 corpus archive SHA-256 remains
`39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`.
The successful command was:

```sh
python3.12 research/notes/ecc2k130/rotated_m56_export_gate_20260925/run.py \
  --out /private/tmp/kic-rotated-m56-run-20260925-a1
```

The two children ran sequentially on the local Mac from 2026-09-25
12:43:58.264469 through 12:44:19.554498 UTC. Both exited zero. The exact
stdout, empty stderr, deterministic JSON outputs and `receipt.json` are
committed here byte-for-byte. The receipt preserves each child's command,
start/end UTC, exit code and SHA-256 of each raw file. Its SHA-256 is
`dfe4a769962ba3c026dec33d5f8208cf6bfa459d3c3d0221b33c41d56412ba87`;
`n13-m5.json` is
`d0565f4f6cd01ba275aa6949c3ed76308c64f4223d2978ebcc2ce8acb3cfae40`,
and `n19-m6.json` is
`fd47e2a3de0c79b406ebb5090e1d3bc784870166ba83add0cd29170a262d67b3`.
Each child's JSON records the frozen source and input hash, point/field
operation counts, exact per-mask rational lift counts and witnesses, S3/branch
and parity counts, two structural ordering negative controls, wall, CPU and
peak RSS. No field-operation-to-rho conversion was measured.

A fresh deterministic replay against the committed raw files passed locally
immediately afterward:

```sh
python3.12 research/notes/ecc2k130/rotated_m56_export_gate_20260925/ci_replay.py \
  --evidence research/notes/ecc2k130/rotated_m56_export_gate_20260925/evidence
```

It re-enumerates both exact arms and compares every deterministic result field
including group-law counts, target-model witnesses, parity and chain branch
counts, structural ordering controls and native operation counts. It excludes
only wall, CPU and RSS from cross-host equality; the committed values are
still checked against frozen caps. The focused GitHub workflow repeats this
replay on the exact PR head. This is a replay of one source implementation
using #767's pinned independent group law, not a new independently designed
arithmetic implementation. Its counts also cross-check #767's separately
produced complete full/projected histograms.

The first 12:40:53 UTC attempt is not hidden: `failure_0/` contains its
receipt, child JSON/stdout/stderr and explanation. It ran no point tuple and
is not combined with or substituted for the successful receipt.
