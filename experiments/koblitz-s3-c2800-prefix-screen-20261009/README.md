# Fresh N53 C1200/C2800 prefix screen

The C2800 arm completed T149 with a verified scalar, but its measured
whole-process peak RSS was **36,839,751,680 bytes**, exceeding the frozen
32 GiB envelope by **2,480,013,312 bytes**. The paired C1200 run peaked at
16,147,365,888 bytes. Both arms used **374** target S3 calls on T149; its
single-pair local online ratio was **0.638×**. The resource gate was already
false, so the panel was stopped during the next C2800 run rather than
continuing fourteen more expensive launches. T150 expanded is preserved as
an explicit interrupted row. The aggregate screen call and online ratios
remain **unknown**, and the C2800 configuration does not advance under this
memory envelope. The earlier T025 feasibility run peaked at 30,260,428,800
bytes; the separate per-run receipts show why a single memory control was
insufficient to qualify this panel.

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
the resource stop, [`preserve_screen.py`](preserve_screen.py) archived the
three started rows and available raw reports.
[`analyze_resource_stop.py`](analyze_resource_stop.py) checks their identities,
pairing, report hashes, scalar results, peak RSS, runner exit, and the stopped
gate to regenerate [`analysis.json`](analysis.json). The decompressed ledger
SHA-256 is
`c45d48803188adac2e276d79e0cfe81c1c71e32880e98489f4843d6b6e0ff1ca`;
the raw-report archive SHA-256 is
`421694c939a2cb036a50b59f3524bf9b3637a1bae13167960589396559d59708`.

The next candidate should reduce the C2800 root table's slot payload before
another fresh-target panel. The current table uses 13 bytes per slot across
1,073,741,824 slots; an 8-byte packed slot would save **5,368,709,120
payload bytes** if it preserves the exact key/value and collision semantics.
This is a storage calculation, not a measured peak-RSS or online gain.
