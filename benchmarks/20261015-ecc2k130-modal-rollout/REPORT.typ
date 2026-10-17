#set page(paper: "us-letter", margin: 0.65in)
#set text(font: "Arial", size: 10pt)
#set heading(numbering: none)

= ECC2K-130 Modal σ-table rollout

_Descriptive visual record, 2026-10-15 UTC. This report records a matched kernel benchmark and a later operational rollout. It makes no campaign attack-cost claim._

== Matched kernel measurement

The frozen RTX PRO 6000 headline panel measured the fused-control median at *15.513653 billion scalar updates/s* and the σ-table median at *15.884369 billion scalar updates/s*. The median of five paired A/B ratios was *1.023930* (range 1.023803–1.024086); maximum A/A drift was 0.099%. Both arms passed 300/300 reference replay, with identical sorted distinguished-point corpora. The chart zooms the rate axis to display the measured difference; it does not depict a whole-campaign rate.

#image("gpu-rate.svg", width: 100%)

== Live deployment observation

Two Modal RTX PRO 6000 validations passed with the selected B16/T512/minBlocks1 build and packed GF(2^131) replay. The deployed image enabled inline polynomial mode 3, fused σ scheduling, and the shared σ table. Fresh call 8011 reported *15.8931 billion iterations/s* at 5m15s and uploaded distinguished points and a checkpoint. Calls 8009 and 8010 continued their prior image; the supervisor is configured to use the new image after each existing four-hour pass completes. That future switch was not observed in this record.

The public status snapshot at 00:36:01 UTC reported four walking slots, 532,964,467 total distinguished points, and zero outstanding ingest objects. This live client rate is descriptive: the executable hash was not emitted, and the two validation builds had different executable hashes despite equal source hashes.

== Evidence and limits

- Matched data: `ecc2k130/benchmarks/sigma-table/headline/result.json` and `samples.tsv`.
- Rollout and raw logs: `ecc2k130/benchmarks/modal-sigma-table-live-20261015/receipt.json`, `modal-validate.log`, worker and uploader excerpts, and `SHA256SUMS`.
- Test receipt: `ecc2k130/benchmarks/modal-sigma-table-live-20261015/release-test-receipt.json`. The parallel sandboxed full suite failed 13 tests; all 13 passed in targeted reruns with socket access or one test thread. The original full-suite result remains failed.

The matched panel establishes a local kernel throughput comparison on the specified card and workload. The rollout confirms that one new call ran the selected flags and produced uploads. It does not establish an end-to-end solve-cost change or the later image switch for the two older calls.
