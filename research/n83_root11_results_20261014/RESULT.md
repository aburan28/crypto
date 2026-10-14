# N83 reconstructed second root: 564,400 usable factor-base points

The 11-byte retained state stores one canonical root and reconstructs its
partner from the indexed pair and root using the characteristic-two sum and
product identities. With fixed `2^30` exact five-byte root slots, the cloud
run completed **K3400 with 564,400 subgroup-usable factor-base points** and
3,400 folded columns. Its 959,460,734 exact root keys occupied 89.36% of
the allocated table. The final `wait4` process peak was 11,112,914,944
bytes (10.35 GiB), within the 14-GiB observed-peak guard. The independent
cloud checker replayed each complete base, subgroup point, four-point target
relation, scalar and exclusive online phases. All five frozen arms passed.

| Arm | Usable points | Folded columns | Root keys | Target probes | Online ms | Setup ms | `wait4` peak GiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K3400 root11 | 564,400 | 3,400 | 959,460,734 | 49,703 | 28.364702 | 484,870.239 | 10.350 |
| K3300 root11 | 547,800 | 3,300 | 903,852,075 | 41,137 | 17.561242 | 424,843.768 | 10.050 |
| K3200 root11 | 531,200 | 3,200 | 849,903,212 | 23,693 | 9.562454 | 400,991.080 | 9.762 |
| K2543 root11 | 422,138 | 2,543 | 536,737,856 | 16,174 | 6.281393 | 244,013.602 | 8.049 |
| K2543 root21 | 422,138 | 2,543 | 536,737,856 | 16,174 | 6.639159 | 256,560.031 | 10.544 |

The K2543 controls have identical base hash, exact root-key count, target
probe count and scalar. The root11 state allocation was 2,952,130,555
bytes versus 5,635,885,605 bytes for root21; the complete process peak
fell by 2,679,210,984 bytes. The same-point strong-rho arm recovered and
replayed the scalar in 1,225.661793 ms online. These online intervals are
individual cloud-instance observations on the previously used public point.
K3400 adds 26,560 usable points over K3240 in the preceding 21-byte run.
The five-point secondary panel will select the base by target-search cost
before the fresh-point Tier 1 paired run.

`materialize.py` checks the raw ZIP, independently replayed analysis digest,
executed source, workload and complete factor-base dumps. It emits immutable
`IC1` candidate manifests, the one-point workload and measurement rows with
exclusive phases, process costs, resource counters and correctness status.
The paired K2543 manifest records the 21-byte control as a separate candidate.

## Proof and evidence

`../n83_root11_20261014/RootSum.lean` checks the root reconstruction identity
and the cancellation step needed to recover the second root from the first
when the first root is nonzero. `PROOF_MAP.md` records exact hypotheses,
remaining links to the implementation, and the Lean command and axiom audit.

[Cloud run 38522628671](https://github.com/aburan28/crypto/actions/runs/38522628671)
executed commit `f7aff3a7d83b602df2024c28d596be37fb8c698f`, source SHA-256
`07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae`,
and public point SHA-256
`21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
GitHub Actions artifact ID `11876331881` has SHA-256
`ab8bf4014f6d4d4be185e8e0e6bb9bbeb2882896aacf5bd540a1afa8468633ad`.
Its `results/analysis.json` SHA-256 is
`e1fb189923f11d009978017fe3c8e8806a90ba6ffcd6de065aea0c6cfbabc39b`.
Download the artifact, extract it, then run
`python3 materialize.py --artifact-zip PATH_TO_ZIP PATH_TO_EXTRACTED_RESULTS`
to check the raw archive and regenerate the candidate and workload records.
