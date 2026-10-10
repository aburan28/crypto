# N53 C2800 prefix feasibility

The C2800 prefix candidate
`IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1f41cb1aeaa9`
passed its frozen feasibility control. It constructed **296,800** distinct
usable factor-base points and **2,800** full-rank relation columns, built
**415,520,000** regular states and a **1,073,741,824-slot** root table,
then recovered and independently verified the T025 scalar with **498**
target S3 calls and zero fallback attempts. Observed child peak RSS was
**30,260,428,800 bytes** under the 32 GiB envelope, and process wall time
including reusable preparation was **97.464 seconds**. These values are a
single correctness and resource control; the next measurement is a frozen
fresh-target paired screen against C1200.

The C2800 candidate searches the exact C1200 target-state prefix against a
larger root table, then traverses the complete C2800 index if needed. Its
source is committed at `a1df596ee` and the release example passes all eight
unit tests, including the prefix/fallback permutation and ordered S3 root
checks.

[`FEASIBILITY_FREEZE.json`](FEASIBILITY_FREEZE.json) fixes the N53 curve,
T025 correctness-control point, 2,800 prepared relation columns, 415,520,000
regular states if all pairs are regular, 831,040,000 root inputs, a maximum
1,073,741,824 table slots, a 32 GiB observed peak-RSS envelope, and a
900-second cap. This run measures feasibility and scalar correctness; its
target timing is a diagnostic. The target point and source/binary hashes are
checked before the process starts. Successful feasibility requires full rank,
verified scalar replay, no fallback, table slots within the declared step,
and observed peak RSS within the envelope. Failed runs and raw reports are
preserved before any candidate identity is generated.

The first two launches ended without result receipts after their controlling
turns were interrupted. [`attempt1_observation.json`](attempt1_observation.json)
and [`attempt2_observation.json`](attempt2_observation.json) retain what was
observed; both results and peak memory remain unknown. The third launch uses
the unchanged frozen source, binary, and runner with receipt paths on the
persistent SSD volume.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c2800-prefix-20261009/run_feasibility.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2800prefix_targetspan \
  --receipt /Volumes/SSD990/llm/tmp/s3-c2800prefix-feasibility-receipt.json \
  --raw /Volumes/SSD990/llm/tmp/s3-c2800prefix-feasibility-report.json.gz
```

Then run [`record_feasibility.py`](record_feasibility.py) to preserve the
receipt/report and, only on success, derive the complete candidate manifest
and fresh T149–T156 workload records for a separately frozen paired screen.
The checked [`feasibility_receipt.json`](local/feasibility_receipt.json),
compressed [`feasibility_report.json.gz`](local/feasibility_report.json.gz),
and canonical [candidate manifest](identities/candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h1f41cb1aeaa9.json)
retain the exact source, binary, base digest, rank, phase costs, and scalar
certificate. The raw report SHA-256 is
`4d65f2fbfcaf3802954db16ebba173edbdaeb2201a57a8e481600fd8547a7ba8`.
