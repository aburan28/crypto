# Width 6 paired complete-call result

The preregistered native one-thread Apple ARM64 follow-up completed all four
frozen seeds. Each seed had one warmup per arm, five width-8/width-8 A/A pairs,
and five alternating width-8/width-6 pairs. All 88 calls completed. Both
arms used the exact cut-19 support certificate and returned identical
non-timing output, including rank, canonical and raw fingerprints, and
original rows. The six smaller cases also matched. Width 6 used 86.15–86.32%
of width 8's counted reduction word XORs, passing the 90% work gate.

The table gives width-8 time divided by width-6 time for the n24 degree-4
**complete call**. Intervals are exact five-pair bootstrap percentile bounds.
A/A is the observed width-8/width-8 ratio range on the same seed.

| Seed | Width-6 work / width-8 | A/A range | Paired median | Paired 95% bootstrap interval |
| --- | ---: | ---: | ---: | ---: |
| frozen | 86.15% | 0.531–1.212 | 1.155 | 0.730–1.718 |
| holdout_a | 86.26% | 1.006–1.289 | 1.027 | 0.892–1.131 |
| holdout_b | 86.20% | 0.900–1.131 | 1.141 | 0.375–3.567 |
| holdout_c | 86.32% | 0.808–1.016 | **0.903** | 0.863–0.928 |

No primary paired median fell below its seed's A/A minimum, so width 6
passes the preregistered gate to a physical Linux test. No smaller-case
median fell below its A/A minimum. The holdout-c median nevertheless shows
a local slowdown, and these Apple timings do not establish a speedup.
Width 6's reliable finding here is lower counted work with exact output;
whether it improves complete-call wall time remains unverified.

The four full per-call receipts, including process status, source and binary
hashes, hardware, complete and phase times, and output, are
`PAIRED_<seed>.json.gz`. `PAIRED_RAW_SHA256.txt` hashes the uncompressed
JSON; `PAIRED_ARCHIVE_SHA256.txt` hashes the deterministic gzip archives.
The native harness source is `examples/f5_support_rank_bits_pair.rs`.

The unchanged requested acceptance gate is a physical isolated Linux x86-64
one-thread comparison against selective echelon: each of four seeds needs
a complete-call median and exact bootstrap lower bound above 2.00×, with
six smaller-case and two-thread nonregression controls. No such comparison
has run for this width-6 variant. These are F5 solver-stage diagnostics,
not an IC online or rho result.
