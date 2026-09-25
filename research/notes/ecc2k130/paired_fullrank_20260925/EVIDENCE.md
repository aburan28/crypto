# Provisional paired panel, 2026-09-25

This evidence is a correctness and cost control on a **source overlay**, not a clean-main performance verdict. It reconstructs exact pinned-main producer and arithmetic blobs on checkout `74cf8424`; `source_provenance.json` lists all 21 overlaid Git blob IDs and SHA-256 values. The reproducible clean-main workflow in this PR is the promotion gate.

Each arm received the same public `hash:SEED` target. Every pair emitted identical Q coordinates and the same recovered scalar; the independent Python verifier replayed every relation, row, rank transition and final solve, all base-point orbit labels, and `[d]G=Q`. All six paired runs passed, with 820 relations replayed. Altering one target key or rho scalar caused replay to fail.

| n / seed index | IC wall ms | rho wall ms | rho/IC wall | IC CPU s | rho CPU s | rho/IC CPU | IC peak RSS MiB | relations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 37 / 0 | 424.3 | 254.9 | 0.601 | 0.039 | 0.011 | 0.292 | 10.5 | 15 |
| 37 / 1 | 264.3 | 129.4 | 0.490 | 0.038 | 0.024 | 0.634 | 10.6 | 13 |
| 37 / 2 | 104.5 | 46.7 | 0.447 | 0.038 | 0.015 | 0.388 | 10.4 | 13 |
| 41 / 0 | 19112.8 | 811.4 | 0.042 | 9.914 | 0.450 | 0.045 | 1700.1 | 271 |
| 41 / 1 | 17768.8 | 600.6 | 0.034 | 9.625 | 0.331 | 0.034 | 1700.5 | 219 |
| 41 / 2 | 21848.8 | 649.6 | 0.030 | 10.075 | 0.297 | 0.029 | 1673.8 | 289 |

The ρ walk used the signed-Frobenius quotient (`A=74` at n37; `A=82` at n41) and 0.32–1.10 times the ideal collision-step estimate across these six runs. No rho failure or health outlier explains the direction. Wall at n37 varies substantially relative to CPU; it cannot justify a precision claim. All costs above are fresh-process elapsed/CPU usage and include JSON emission and verification, rather than a producer-only stage timer. The direct fixture recomputes full rank for every Q in its batch mode; this is not the compact producer's shared-log batch.

`parts_manifest.json` gives five Git-tracked binary parts, byte counts and SHA-256. Reassemble exactly as stated there, check the archive SHA-256, extract, then run `shasum -a 256 -c SHA256SUMS` inside it. The archive contains full raw JSONL, stderr, replay receipts, per-run manifests, input spec, generated Cargo.lock and build log. Neither a local path nor the expiring CI artifact is the sole evidence.

Classification: negative control for this point-defined baseline, provisional performance evidence. Do not update the canonical crossover scoreboard from this source overlay alone; require the clean-main panel and a separate fully charged shared-log comparison.