# Batch compact-orbit producer failure-status control

Status: frozen before the source edit or regression runs. The current
`examples/koblitz_orbit_dlp_s3_batch.rs` on main `5874ce8b3` reports every
target record with `exit_code:0` and returns process success even if
`group_verified` is false or null. This contract defect can hide a failed
one-target result from a shell runner. The hypothesis is that a minimal
status-only correction will preserve all rank/target decisions and raw
records while reporting any failed target through both the per-target field
and process exit status. This is a correctness fix, not a timing comparison,
new factor-base result, or speedup claim.

The frozen inputs are already public n41 points from the merged
`koblitz-base-size-cold-panel-20261004` study. `failed_pilot_q.jsonl`
is `[1119268120096,860021860491]`, SHA-256
`4db7c82489a620f5aaeee8835cd0cdfc076ee3295adbedd0feec10dc315e0c97`.
The prior **different** compact-orbit producer missed it at K=20; whether
this S3-batch producer also misses is unknown at freeze. The successful
control `success_holdout_q.jsonl` is
`[1071506060992,898053054019]`, SHA-256
`6af48a03b161422fa3ed8f58df4eef78e49e2abf21f808775df67c0817895db4`.
No known scalar is passed to either producer. Both use the same source
Koblitz n41 curve, `construct:41:0:20` for the failure control and
`construct:41:0:85` for the success control, rank seed `410041`, one thread,
`KIC_S3_BATCH_WINDOW=64`, `KIC_S3_PREFILTER=blocked`, and
`KIC_QUERY_BACKEND=point_sum`. Do not change inputs or `K` after outcomes.

Build the exact old main source and corrected source with the same offline
Cargo lock, save both binary/source hashes, and run one fresh process per
arm under a 45-second cap. Archive stdout summary, target JSONL, stderr,
exit status, source/lock hashes, and exact command. A completed failed-Q
control should have the same public point and failed target decision in
both arms; the old process may return 0, while the corrected process must
return 1 with `target.exit_code=1` and `targets_failed=1`. The success-Q
corrected process must return 0 with `target.exit_code=0`, one verified
scalar, and `targets_failed=0`. Preserve a timeout, unexpected success, or
other mismatch as its own outcome; do not replace the Q or claim that the
process-level failure path was exercised. An independent code-level test
must cover true, false and unknown verification states.

No rho run or CPU ratio belongs to this control. Later fresh-Q point-sum
evaluation must independently replay every relation/scalar and use a
separate preregistered full-cost same-Q rho comparison. A process exit alone
is necessary failure accounting, not cryptographic verification.
