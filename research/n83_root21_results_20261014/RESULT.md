# N83 fixed-slot 21-byte roots: 537,840 usable points within 14 GiB

The 21-byte root-pair layout and fixed `2^30` exact five-byte root slots
completed **K3240 with 537,840 subgroup-usable factor-base points**,
3,240 folded columns and 871,283,595 exact root keys. The final `wait4`
peak was 14,954,881,024 bytes (13.928 GiB), **77,504,512 bytes below**
the 14-GiB observed-peak guard. The independent cloud checker replayed
every complete base, its subgroup points, four-point target relation,
scalar and exclusive online phases. All five frozen arms passed.

| Arm | Usable points | Folded columns | Root keys | Target probes | Online ms | Setup ms | `wait4` peak GiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K3240 root21 | 537,840 | 3,240 | 871,283,595 | 35,691 | 25.020425 | 693,871.641 | 13.928 |
| K3220 root21 | 534,520 | 3,220 | 860,560,226 | 93,347 | 57.618809 | 621,142.873 | 13.822 |
| K3190 root21 | 529,540 | 3,190 | 844,599,620 | 74,825 | 44.613202 | 598,232.887 | 13.660 |
| K2543 root21 | 422,138 | 2,543 | 536,737,856 | 16,174 | 8.661459 | 358,540.077 | 10.545 |
| K2543 root22 | 422,138 | 2,543 | 536,737,856 | 16,174 | 8.806532 | 357,788.318 | 10.795 |

The K2543 controls have identical base hash, exact root-key count,
target probe count and scalar. The same-point strong-rho arm recovered and
replayed the scalar in 1,691.535783 ms online. The listed online intervals
are individual cloud-instance observations on the previously used public
point. `measurements.json` retains the canonical candidate and workload
identities, exclusive phases, full process costs, resource counters and
correctness for each arm. K3240 adds 20,916 usable points to the earlier
in-envelope K3114 result. The archived five-point secondary panel will
select a base by target-search cost before the fresh-point Tier 1 paired run.

## Evidence and replay

[Cloud run 38520546294](https://github.com/aburan28/crypto/actions/runs/38520546294)
executed commit `d8e679f137dd13037cfe2b1471be3cfa76d52eb3`, source SHA-256
`accd99313204d0e3327968230ed3c64df49d661ebddb466147a9a8fc844bd1c9`,
and public point SHA-256
`21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
GitHub Actions artifact ID `11875926923` has SHA-256
`7bca69a792f0f686aab8f3a5b2c3fec0546bb98643807f36056d23ce5fb8da2b`.
Its `results/analysis.json` SHA-256 is
`bd6c1d65343903f584e2336372424aad4e65cb0a0d2b5177572396db9802b4e2`.
Download that artifact, extract it, then run
`python3 materialize.py --artifact-zip PATH_TO_DOWNLOADED_ZIP PATH_TO_EXTRACTED_RESULTS`
to check the raw archive and regenerate the canonical `IC1` candidate and
one-point workload records. The complete cloud verifier passed before the
artifact was uploaded.
