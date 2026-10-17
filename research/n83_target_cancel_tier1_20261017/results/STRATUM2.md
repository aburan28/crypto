# N83 cancellation Tier 1 stratum 2

Cooperative two-worker target cancellation was faster than two-worker join
on all six sealed public points assigned to stratum 2. The geometric
target-online times were 39.286846 ms for cancellation and 125.404958 ms
for join, a paired join/cancellation ratio of 3.192034. The geometric
same-point rho/cancellation ratio was 64.131893 across the six separately
solved rho controls. This is an admitted six-point stratum; the predeclared
24-point exact selector waits for the other three strata.

| Point index | Join online | Cancellation online | Join/cancellation | Rho paired with cancellation |
| ---: | ---: | ---: | ---: | ---: |
| 12 | 112.056458 ms | 73.368028 ms | 1.5273 | 1,608.147632 ms |
| 13 | 41.785380 ms | 29.418040 ms | 1.4204 | 2,900.855846 ms |
| 14 | 554.866819 ms | 129.964157 ms | 4.2694 | 2,223.936485 ms |
| 15 | 228.160473 ms | 81.023166 ms | 2.8159 | 2,564.564154 ms |
| 16 | 150.335101 ms | 4.715396 ms | 31.882 | 2,803.713929 ms |
| 17 | 43.645295 ms | 34.309269 ms | 1.2721 | 3,429.301374 ms |

The cloud verifier admitted the host with `PASS_TIER1`,
`PASS_INDEPENDENT_REPLAY`, passing resource and event-order checks, zero
observed steal or OOM gates, and before/after rho calibration within 10%
for every pair. An independent local replay of the provider ZIP exited 0
with empty stderr and reproduced the cloud `analysis.json` byte for byte
(SHA-256
`6a1dbbf20bd5fe4c2c1ebba534da00bafd93db10c947fddfe5bedfc18cce79a5`).
The exact [analysis](stratum2_analysis.json),
[provider receipt](stratum2_provider_receipt.json), and
[local replay receipt](stratum2_local_replay_receipt.json) are versioned here.
The provider artifact ID is `12108835537`; its 49,901,327-byte ZIP has
SHA-256 `b1d0248f993b8fd6deefaec5cfa024d65db9ba1d46c3d1135860402f87115181`.
Full raw records remain under
`/Volumes/SSD990/crypto/artifacts/n83_target_cancel_tier1_20261017_raw`.

Both methods use the same K992/root11 factor base with 164,672 usable points
and 992 folded columns. Candidate IDs are
`IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h33c9c3b2f933`
for join and
`IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h58ec337f71f7`
for cancellation. Every target-online interval includes query,
decomposition, relation check, descent, and producer recovery check after
target-independent setup. The full host and per-point cold costs remain in
the analysis JSON; the 24-point aggregate will report the final exact
paired selector and same-point rho ratios.
