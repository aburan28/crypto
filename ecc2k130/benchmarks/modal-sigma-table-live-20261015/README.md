# Modal σ-table campaign rollout, 2026-10-15

This is an operational receipt for selecting the previously verified σ-table
kernel in the live ECC2K-130 DP32 Modal campaign. The source change is commit
`0e62889d7`; the GPU validation, deployed worker sample, uploader excerpt,
call IDs, launcher scripts and public status snapshot are retained here.
`receipt.json` states their observation boundary and file hashes are in
`SHA256SUMS`.

The selected kernel preserves the campaign's σ iteration function, 32-byte
DP format and packed checkpoint state. Its prior matched RTX PRO 6000 panel
measured 15.884369 B/s against 15.513653 B/s for the fused control, a
1.023930 paired median ratio, with an A/A maximum drift of 0.099%, 300/300
reference report replay per arm and identical sorted corpora. Those are
benchmark observations under the frozen protocol in
[the headline result](../sigma-table/headline/result.json) and
[the interpretation](../../SIGMA-TABLE.md), not a measured campaign solve-cost
change. The faster 20 B/s table walk has a different iteration function and
cannot be substituted into this corpus.

`CURVE=131 ./run.sh validate` passed on a Modal RTX PRO 6000 with sm_120,
CUDA compiler 13.3.73, the selected B16/T512/minBlocks1 binary and packed
GF(2^131) report replay. The first successful validation is Modal app
[`ap-Wy8WgzfwdoGA7egD5l87KX`](https://modal.com/apps/a-buran28/main/ap-Wy8WgzfwdoGA7egD5l87KX).
The first attempt, app `ap-nXd0dTEaYoCjvzoaWYccQe`, was interrupted with
exit 130 before validation when the launcher still passed minBlocks2; that
launcher argument was corrected. The subsequent validation passed, and the
complete repeat output is retained as `modal-validate.log`. That repeat ran
under app [`ap-H7n0QUIUhHERcq2Iu1Ov4J`](https://modal.com/apps/a-buran28/main/ap-H7n0QUIUhHERcq2Iu1Ov4J).
Its source hash matches the first validation, while its executable SHA-256
differs; the cause has not been established. The live call printed the
selected kernel flags but did not emit its executable hash, so the live rate
is kept descriptive.

The deployed app is [`ecc2k130`](https://modal.com/apps/a-buran28/main/deployed/ecc2k130).
The supervisor reattached to original calls 8009 and 8010 without starting
duplicate owners and spawned fresh run 8011 as
`fc-01M511XE04ND5GK8CBY6HGYGHR`. Its log names campaign slot 98011,
DP weight 32, WITNESS=0, `sigma fused: 1`, `sigma table: 1` and
`inline polynomial: 3`. At 5m15s it reported 15,893.1 M iterations/s and
13.58k DPs. This is a live client sample on one card, not a paired comparison.
The uploader sent 3,343 records and a 40-byte checkpoint for slot 98011 at
00:33:00–00:33:01 UTC, then more batches. The public feed at 00:36:01 UTC
reported four walking slots, 532,964,467 total DPs and zero outstanding
objects.

Calls 8009 and 8010 were still running their previous image when the
supervisor was changed. They retain their checkpoints and are configured to
use the deployed σ-table image at their next successful four-hour pass.
That future resume has not been observed in this record. The two local tmux
sessions remain `ecc2k130-modal-fleet-20261015` and `ecc2k130-modal-sync`;
the launchers saved here use the existing call state and uploader offset state
under `/Volumes/SSD990-2/cryptanalysis/ecc2k130/runner/research/production/2026-10-15-modal-restart-runtime/`.

No attack-cost or end-to-end performance claim is made from this rollout.

## Local release test gate

The first `cargo test --release --lib` run finished with 4,116 passed, 13
failed and 117 ignored. Nine Pollard collaboration tests failed because this
sandbox refused local socket operations (`Operation not permitted`). Four
isogeny-cover tests timed out while the suite ran concurrently. The original
output is retained in `release-suite.log`; the full-suite outcome is failed.

The Pollard collaboration group was rerun with loopback access and one test
thread: 37 passed. The isogeny-cover group was rerun with one test thread:
14 passed and one intentionally ignored, including all four original timeout
failures. The exact rerun outputs are `release-pollard-collab-rerun.log` and
`release-jv-cover-rerun.log`. `release-test-receipt.json` records the commands,
outcomes and scope. These reruns explain the 13 failures but do not turn the
original parallel, sandboxed full-suite run into a green result.
