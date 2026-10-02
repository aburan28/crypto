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

## Discovery 02: complete replay, incomplete strongest-control roster

The corrected Linux ARM64 workflow [36982536433](https://github.com/aburan28/crypto/actions/runs/36982536433)
passed correctness, the reserved campaign, first verification and post-seal
replay at source `bee44bf317d6294a2b99db240f04ed644c551539`. The
downloaded 27-file artifact is preserved byte-for-byte in
`comparator_incomplete_discovery_01`. All 26 member hashes match manifest
SHA-256 `5a4594b1022270d49af8e6303b78a59b8f9f2cd4822261e0d04639eeda603afa`;
raw SHA-256 is `b400af9978144c8ab9ad1ca07b13b300d6429e4b82adcd5fcf8a5b9b2ee923d3`;
result SHA-256 is `bb8bd89bc582daee8fbd95e4cee19810cfcb4458a5fcd10800952bb0dbd18916`.
Artifact ID `11216762190` has archive digest
`sha256:b75fb2f0627247baeaf7e2efba8173aeb90740541ecb974e1e3c4312b95add8a`.
An untimed replay with the same source on another host also passed.

It verified 200 fixtures, 13,600 measured arm batches and 291,040 outputs.
The one-core receipt was uncontended: worker wall 242.548990231 seconds,
other-process CPU 2.36 seconds, and preflight PSI some avg10 3.64. Its stored
eight-group gate was false, with only the two n=24 groups above the 1.5x
discovery lower bound. No holdout was launched.

**Comparator correction:** the dense packed reference in that source stored
one 32-bit entry for each n-bit monomial mask. At n=24 its table alone used
`2^24 * 4 = 67,108,864` bytes, exactly the frozen 64 MiB context cap;
the basis pushed the retained context above the cap, so the worker marked it
inapplicable. Every ambient column index fits in 15
bits; a 16-bit entry with one high/low flag uses only `2^24 * 2 = 33,554,432`
bytes plus the basis and fits the same cap. Thus the original n=24 ratios
near 2.9x were against a weaker reference than the protocol's *fastest
applicable correct control*. They cannot be promoted as strongest-control
gains, despite successful resource and mathematical replay. No frozen file,
raw sample, limit or gate is rewritten. The source now uses exact 16-bit slots
and tests every ambient coordinate, including n=24; a third fresh discovery
must measure it in the same binary. The holdouts are still unused.
