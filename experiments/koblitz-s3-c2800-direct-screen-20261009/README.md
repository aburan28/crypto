# Fresh N53 C1200/C2800 direct-fill screen

The C2800 direct-fill candidate passed the frozen T025 scalar, rank, semantic,
and 32 GiB resource gates. This protocol compares it with a C1200 baseline on
eight previously unused, identical one-target public points T157–T164. Both
arms run from the same frozen binary. The baseline candidate is
`IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0h78c2fc9eed5a`;
the expanded candidate is
`IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0hf4206807ca48`.
Their manifests retain the exact bases, state-generation paths, root-table
policies, source identity, and target schedules.

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) fixes AB/BA order, one run per arm
per target, 14 preparation workers, one target worker, a 600-second per-run
cap, and 32 GiB maximum whole-process RSS. The runner verifies the frozen
source and binary, candidate and workload hashes, public points, full-rank
preparation, phase sums, scalar replay, and exact table/base policies. It
writes an append-only row and raw report for every started run. It stops after
an invalid result or a resource-gate failure, preserving the failing row and
`stop.json`. A complete 16-run screen advances to independent confirmation
only if every scalar and resource gate passes and the geometric mean paired
C1200/C2800 target S3-call ratio reaches 3. The user's final large-gain
threshold is at least 2× verified one-target online wall time under controlled
CPU conditions; this local screen is an algorithmic and feasibility gate.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c2800-direct-screen-20261009/run_paired.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800direct_targetspan \
  --output /Volumes/SSD990/llm/tmp/s3-c2800-direct-screen-20261009
```
