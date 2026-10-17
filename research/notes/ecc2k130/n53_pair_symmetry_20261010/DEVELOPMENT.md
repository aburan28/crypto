# Pair-swap quotient development controls

The corrected n53 index makes 1,282,710 S3 construction calls for the
23,320-point, 220-column base, compared with 2,565,200 ordered calls. It
restores every swapped state at its original scan position, preserves all
2,564,528 root-table entries, and skips the later member of a pair during
each circular query scan. The complete n53 root-table key/value test passes.

The 39 raw development files are preserved losslessly in
[`controls.tar.gz`](controls.tar.gz). `controls_manifest.json` gives each
file's SHA-256 and size; `python3 archive_controls.py` replays the archive
check. Archive SHA-256:
`92fe4d43aa552e91240c9aa5d8992c43dbf7712dd2498739e4134569291c9539`.
Each archived rank trace and target output has an independent replay receipt.
These controls preceded the frozen held-out point and are not part of the
preregistered six-seed decision.

| Development arm | Index S3 calls | Ordered states | Rank support probes | Rank S3 calls | Complete cold ms | Index ms | Result |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| n53 initial ordered | 2,565,200 | 2,565,200 | 18,413,613 | 9,206,863 | 6,399.634 | not separately recorded | rank 220; scalar replay PASS |
| n53 initial compacted pairs | 1,282,710 | 1,282,710 | 20,510,895 | 10,255,503 | 4,892.440 | not separately recorded | rank 220; scalar replay PASS |
| n53 corrected ordered | 2,565,200 | 2,565,200 | 18,413,613 | 9,206,863 | 4,809.736 | 925.318 | rank 220; scalar replay PASS |
| n53 corrected paired | 1,282,710 | 2,565,200 | 18,399,197 | 9,199,655 | 4,388.170 | 458.956 | rank 220; scalar replay PASS |

All n53 rows use the prior development point
`[7960849849661793,7443722527872608]` and rank seed 530053. The initial
compacted prototype changed the scan origin and increased rank work, so its
cold timing cannot isolate the index improvement. The corrected pair agrees
on all 220 rank relation rows, base digest, and independently replayed target
scalar `17385600002971`; its 14,416 fewer support probes and 7,208 fewer
rank S3 calls arise from skipping symmetric repeats after the same circular
origin. The corrected control and paired target online times were 4.960 and
5.029 ms, respectively, on a shared host. Host wall times are exploratory
under the repository CPU-isolation gate.

An n13 pair with public point `[384,1476]` passed independent rank and
target replay in both arms. Those archived n13 runs exercise the initial
compacted version; the current n13 and full n53 root-table tests cover the
corrected version. During manual setup, two attempted n13 commands supplied
an output JSONL path in place of the public-point input and exited 101 with
`ParseIntError`; those command failures were not captured as durable raw
files, and no performance or correctness conclusion uses them.

The next decision uses a fresh public point and six rank seeds under
[`PROTOCOL.md`](PROTOCOL.md). The candidate must preserve the full rank-row
sequence and target scalar and must not increase rank root calls in any pair.
