# N83 K3072 target-worker screen: capped local serial arm

The frozen K3072 serial arm reached its 1,200-second process cap during
target-independent construction on the shared macOS host. Its producer
status is `TIMEOUT`, with zero target-result bytes, no base dump, no online
interval, and no completed cold process. The cancellation arm was left
unrun because a paired policy comparison could not be formed from this
execution. This cap is a censored local diagnostic, not a factor-base or
target-policy selection.

| Field | Recorded value |
| --- | --- |
| Candidate | `IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0h49e5d0fe6271` |
| Source SHA-256 | `8109fb33374ce2aa66c6a5a0b6e144d3b0143168ac756f9e0bd70af60a911829` |
| Binary SHA-256 | `df8340eb76ac90d2520e0202159505426b9d7879978f10c597731c4bfe73045a` |
| Frozen plan SHA-256 | `d1a8489cdb646b9c2c9bd06d44847600ea7ef14a2d78ab14c62cf418f65d73b8` |
| Outer wall | 1,200.396926 seconds |
| Child max RSS | 3,290,841,088 bytes |
| Target rows | 0 |
| stdout/stderr | Both empty |
| Local host tier | Ordinary shared-host diagnostic |

The versioned [run receipt](local_serial_timeout_receipt.json) has SHA-256
`d2b236e57416b8306badf61d782e65867b798311fe06d6603c98007f947700f0`.
During this run, the shared host's load average was around 90–130 and its
system data volume repeatedly approached full capacity; these observations
limit interpretation of elapsed time. The raw process and empty target file
remain under `/Volumes/SSD990/crypto/artifacts/n83_k3072_target_cancel_20261017/serial`.
The sealed fresh-point holdout remains unused.
