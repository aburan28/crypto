# N83 load75 root index: 464,800-point base within 14 GiB

The five-byte exact root index with a `ceil(4S/3)` table request completed a
**464,800-point N83 factor base at K2800** with 13.836 GiB observed peak RSS.
All five frozen arms recovered and verified the same public target, and the
independent Python checker replayed their complete bases, subgroup points,
four-point relations, scalars, and exclusive online phases on the cloud
runner.
The strongest same-point interval in this sweep was K2543 at 8.710240 ms;
its shorter search on this point motivates a secondary target panel before
the next factor base is selected.

| Arm | Usable points | Root keys | Target probes | Online ms | Setup ms | Observed peak GiB | 14-GiB envelope |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| K2830 dense | 469,780 | 664,725,722 | 177,836 | 107.203147 | 476,243.9 | 14.021 | Exceeded observed peak |
| K2800 dense | 464,800 | 650,707,253 | 83,076 | 45.773536 | 443,638.9 | 13.836 | Within |
| K2750 dense | 456,500 | 627,675,194 | 116,974 | 63.026135 | 424,738.4 | 13.510 | Within |
| K2543 dense | 422,138 | 536,737,856 | 16,174 | 8.710240 | 358,143.6 | 12.294 | Within |
| K2543 regular | 422,138 | 536,737,856 | 16,174 | 8.737751 | 357,822.2 | 12.294 | Within |

The same-point strong-rho reference recovered and replayed the scalar in
1,874.295374 ms online. These are individual cloud-instance observations
on a previously used public point. Their online intervals and exclusive
phase fields are retained in `measurements.json`; controlled timing will
use a fresh point and alternating Tier 1 pairs after base selection.

## Resource and verification receipts

The K2830 producer completed and its point, base, and scalar passed the
independent checker. The periodic RSS sampler saw 14,814,588,928 bytes;
the child's `wait4` peak recorded **15,055,405,056 bytes**, exceeding the
14-GiB guard of 15,032,385,536 bytes by 23,019,520 bytes. Its row therefore
retains mathematical verification and an observed resource-envelope failure
separately. K2800 is the largest completed arm inside that envelope in this
sweep. A future capacity runner must check the final `wait4` peak even when
the sampler did not stop the child.

[Run 38512134514](https://github.com/aburan28/crypto/actions/runs/38512134514)
executed source commit `e31ede13a5cdf131cd52c9b05e08b4c0f18658d7`,
source SHA-256 `ccaa44824ce3fa6587ea886ed0803cd180c641e05ebc3c6efe10ad665452e672`,
and one frozen public point. The source example suite passed 16/16 tests;
the separate release-library job passed. GitHub artifact ID `11872002940`
published SHA-256
`ada8236d645c48a9ad7f5bf31d0253769322800dd0f5de0527a759476d0a3725`.
Download that artifact from the linked run, extract it, then run
`python3 materialize.py --artifact-zip PATH_TO_DOWNLOADED_ZIP PATH_TO_EXTRACTED_RESULTS`
to check the raw archive before regenerating the candidate and workload
records. The runner's
`results/analysis.json` SHA-256 is
`c9dc35239d12c00058f3dcfb1380959c02e8df9702877b6144c627950afe9dfb`.
`materialize.py` binds completed rows to distinct canonical `IC1` manifests
and the one-point workload record using the sealed artifact and analysis digest.

The next selection experiment measures the archived secondary target panel
for K2543 and the largest replayed in-envelope bases from the 26-byte and
root-only sweeps. It will choose the factor-base policy before generating a
new public point for the five-pair Tier 1 online comparison.
