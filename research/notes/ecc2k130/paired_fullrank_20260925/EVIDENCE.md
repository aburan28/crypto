# Paired panel, 2026-09-25

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
## Clean main-derived panel

The clean-source GitHub Actions run [36109718979](https://github.com/aburan28/crypto/actions/runs/36109718979) built the PR merge checkout `a045a0a80d3a8af00fedc85c0342f4e7d026aacd`, proved its pinned-main ancestry and clean tracked source, then finished all six pairs with independent replay. Every IC/rho pair has the same Q and recovered d. The verifier replayed 820 relations and 38,334 point orbit labels in total; all six terminal ranks were full (13 at n37, 146 at n41). Source and executable SHA-256 values, full raw JSONL, stderr, receipts, input commands, process wall/CPU/RSS, and failure status are in the committed [clean evidence archive](../paired_fullrank_clean_evidence_20260925/clean_archive_manifest.json). The archive SHA-256 is `5e01c894c05919a573f2716bbdeace53952cdd009a85969af56aaad591c2c9ed` (1,168,806 bytes).

| n / seed index | IC wall ms | rho wall ms | rho/IC wall | IC CPU s | rho CPU s | IC peak RSS MiB | relations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 37 / 0 | 36.0 | 20.5 | 0.569 | 0.032 | 0.013 | 18.0 | 15 |
| 37 / 1 | 30.5 | 30.8 | 1.012 | 0.026 | 0.022 | 19.3 | 13 |
| 37 / 2 | 30.7 | 20.6 | 0.670 | 0.027 | 0.018 | 18.2 | 13 |
| 41 / 0 | 6736.7 | 575.0 | 0.085 | 6.729 | 0.568 | 1632.2 | 271 |
| 41 / 1 | 6851.3 | 434.1 | 0.063 | 6.847 | 0.431 | 1630.9 | 219 |
| 41 / 2 | 6857.9 | 383.8 | 0.056 | 6.851 | 0.383 | 1632.6 | 289 |

The median rho/IC wall ratio was **0.670 at n37** and **0.0634 at n41**. One n37 seed had near parity in the opposite direction (1.012), so three seeds do not justify a precise small-rung performance estimate. At n41 the direct arm used about 1.71 GB peak RSS and 6.7–6.9 s process wall; each matched rho run used 0.38–0.58 s. The measured point-defined full-rank method does not beat the frozen rho arm on this panel. This says nothing about an amortized shared-log method, which this fixture does not implement.

The earlier clean workflow [36109167304](https://github.com/aburan28/crypto/actions/runs/36109167304) completed all six producer/replay pairs but failed evidence sealing because the default shallow checkout could not prove the pinned ref was an ancestor. Its raw [failure artifact](https://github.com/aburan28/crypto/actions/runs/36109167304/artifacts/10852239614) has digest `sha256:6a8c808d1dcc8e9416c48fb599d9566d4bea16cd2d234691ebb87d3c8faade25`. The workflow now fetches complete ancestry and checks clean status before running pairs; the strict archive gate remained intact. The local source-overlay table above remains provisional and is not substituted for this clean panel.

The committed archive was produced before main advanced through PRs #739 and #738. The final PR-head CI must independently replay it on the current main-derived merge ref and confirm all source hashes before this panel is merged. If any measured source hash changed, rerun the six-pair measurement on that source instead of carrying forward the old timing.
