# Graded reuse attempt ledger

## Failed discovery 01: post-seal float readback

The first Linux ARM64 workflow used source commit
`1c9b658ff9b183176d7f0cae6e8d2fb05b423813` and frozen protocol SHA-256
`7c25cfa898c690c755c8d6a6377dada44314404774f0c30e037de457b8077862`.
Its reserved worker completed 200 cells, 13,600 cold samples and 291,040
oracle-checked outputs. The first native verifier wrote `results.json`; the
resource record was uncontended, with a 243.204428617-second worker wall,
2.17 other-process CPU seconds and preflight CPU PSI some avg10 3.85.

The workflow **failed in the second, post-seal verifier pass**. Re-reading
`results.json` with the default `serde_json` parser shifted two reported
paired-median ratios by one floating-point unit in the last digit:
`1.0247742717975405` became `1.0247742717975403`, and
`1.4431176749304875` became `1.4431176749304877`. The raw observations and
the serialized result were unchanged: an untimed local recomputation from the
downloaded raw file produced a byte-for-byte identical `results.json`. Building
the unchanged worker and verifier source with `serde_json/float_roundtrip`
then replayed all 26 sealed members and the mathematical result successfully.
This identifies a representation/readback defect, not evidence that the
algorithm or raw sample changed.

The original 27-file downloaded artifact is retained unmodified in
`failed_discovery_01`. Its manifest SHA-256 is
`fe542fed0ea4bee987f4828debb7283f86235d4029d2225e202bdae272b2daf2`,
raw SHA-256 is `3408bd4c545b54a7b1631949df9b47ed30435c5cd15e414275d436722692e86a`,
and result SHA-256 is `988aba9c172c11717d29c10554195e71eb888de5141c64ba30b16bcfa0d0e9ca`.
GitHub run [36960153463](https://github.com/aburan28/crypto/actions/runs/36960153463)
uploaded artifact ID `11206899826`; the upload log reports archive SHA-256
`6d54ac7bb33c835147e695b4c3e1ce769ed94aede5b91d6af5c750c5a41ee0f5`.
Every downloaded member hash matches the original manifest. The failure
happened after that manifest was written, so this historical bundle has no
`failure.json`; this ledger and the workflow status supply the missing refusal.

**Admission:** no timing from this failed attempt is promoted. The new parser
regression test checks exact bitwise JSON round trips for both affected
medians. The launcher now retains a second sealed failure record if any future
post-seal replay fails. A fresh discovery using the same frozen workload and
updated source is required. Neither holdout seed has been run.
