# Duplicate upload recovery, 2026-09-21 Pacific

This is an operational repair and an accounting clarification. It establishes
no change to the walk algorithm, collision probability, or cryptanalytic cost.
The aggregate receipt is in `receipt.json`; no point records or credentials
are included.

The public feed at 2026-09-22 00:40:23 UTC matched the reported 5,885,462
same-seed duplicate records from 17 slots. The ingest journal's last increase
was published at 22:22:01 UTC on September 21; the count stayed constant while
new distinct points arrived. A rolling duplicate count is not a current rate
and does not measure repeated GPU computation: an uploader can resend bytes
without any walk being recomputed.

S3 evidence confirmed overlapping uploads under the same stream identity.
Slot 98000 had 810,268 uploaded record positions covering 475,112 positions,
including two differently sized objects starting at offset zero. Comparing
512 records of those two objects produced identical bytes. Slot 94243 had
four offset-zero objects; a 13,568-record prefix comparison was identical.
These observations establish repeated uploads, without attributing every
duplicate in the dashboard to that cause.

The uploader used only a local JSON offset (in `/tmp` by default), and saved
it after uploading the checkpoint. Losing that file, changing hosts, or a
checkpoint upload failure could replay a corpus prefix under a different
object hash. The repair reconstructs the contiguous remote prefix on every
pass, verifies its hashes against the downloaded corpus, and saves point
progress before checkpoint publication. Missing coverage, a changed corpus,
an older volume snapshot, and incomplete listings stop the upload. A one-shot
sync now exits unsuccessfully when any run reports an error.

Live recovery for GPU runs 8000-8003 and CPU runs 9000-9003 uploaded exactly
1,191,589 pending records in eight objects. The next live pass uploaded zero
records, objects, or checkpoints. All 74 focused uploader, launch-rule, and
driver tests passed, including lost state, stale/foreign offsets, overlapping
objects, pagination, prefix mismatch, gaps, and failure after point upload.

Separately, process inspection found legacy Modal runs 4242-4245 still walking
at weight 34 against the campaign's weight 32. Their uploads had stopped, so
the public checkpoint-age test no longer counted them as walking. Each client
received SIGTERM, all four containers finished successfully, and their full
checkpoints and headers were verified on the volume with a 17:46 PDT timestamp.
The legacy app was then stopped at 17:47:47 PDT. The other campaign app
remained running with four tasks. Old corpora and checkpoints were retained.

Run one uploader per run. This patch does not introduce a distributed writer
lease, so simultaneous uploaders can still race after reconciliation. Listing
and hashing the prior uploads is proportional to corpus size; this repair is
not a constant-cost streaming redesign. The patched uploader was used for the
bounded live recovery; no new recurring service or GPU fleet was launched.
