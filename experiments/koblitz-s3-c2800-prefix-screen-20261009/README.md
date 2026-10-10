# Fresh N53 C1200/C2800 prefix screen

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) fixes eight previously unused
one-target public points T149–T156. Each point is paired between the native
C1200 state order and the C2800 root table with an exact embedded C1200
state prefix and complete C2800 fallback. The baseline candidate is
`IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0haeee7ec59767`;
the expanded candidate is
`IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1f41cb1aeaa9`.
The exact manifests distinguish the factor bases, ranks, table capacities,
and frozen source and binary hashes.

The one-repetition screen advances only if all 16 scalar recoveries verify,
every run remains under the 32 GiB peak-RSS envelope, and the target-level
geometric mean C1200/C2800 target S3-call ratio is at least 3. AB/BA order
alternates by target. Each run has 14 preparation workers, one target worker,
and a 600-second cap. `/usr/bin/time -l` records each whole process's peak
RSS, instructions, cycles, and context switches; the runner kills the entire
process group on timeout. These counters cover reusable preparation as well
as target work, whereas the reported online interval contains only the
exclusive target-dependent phases and scalar recovery check. Local online
ratios are exploratory without a host-wide isolation receipt and passing
noise gates.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c2800-prefix-screen-20261009/run_paired.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800prefix_targetspan \
  --output /Volumes/SSD990/llm/tmp/s3-c2800-prefix-screen-20261009
```

The runner verifies candidate/workload hashes, source/binary identities,
public points, full-rank preparation, table and factor-base identities,
search policy, online phase sum, and independently replayed scalars. Its
append-only ledger retains failed rows and complete raw reports. On
completion, [`preserve_screen.py`](preserve_screen.py) archives all rows and
available raw reports, and [`analyze.py`](analyze.py) checks pairing, report
digests, scalar results, resources, and the frozen advancement gate.
