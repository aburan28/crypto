# N83 K3300: five fresh Tier 1 same-point pairs

The selected K3300/root11 factor base solved five frozen, previously unseen
subgroup points on the measured cloud instance. All five alternating IC/rho
pairs passed independent curve, relation, scalar, phase and Tier 1 noise
checks. The median paired online ratio was **113.011×** rho/IC; the
target-resampled 95% interval was **[40.222×, 148.824×]**. This is a
five-target uncertainty interval, not a search-wide guarantee.

## Exact instance and comparison

- Curve: `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)` with
  defining polynomial bitmask `9671406556917033397649557`; model
  `y²+xy=x³+x²+1`; subgroup order `8569786107849059`, cofactor
  `1128547018`.
- Candidate: `IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`.
  Its enumerated base has 547,800 usable points, 3,300 sign/Frobenius-folded
  columns, and source base hash
  `a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
- Five public points SHA-256:
  `2ccaac2baf367598e0d4e724b500ba66c5b4a711f0be8207565ce3acd6b1d0b4`.
  The answer sidecar SHA-256 is
  `5b170d9bc4d34fdd21c032990ceee3a62af335c9af7fcaae2596335b8699dd8c`.
  Each point has its own canonical one-target workload and run ID in
  [measurements.json](measurements.json).
- Allocation: GitHub Actions run
  [38660150265](https://github.com/aburan28/crypto/actions/runs/38660150265),
  Linux x86-64 on AMD EPYC 7763, exposed CPUs `0,1,2,3`, two visible
  sibling pairs, 16,373,444 KiB reported host memory, and a 14 GiB
  process RSS guard. This is **Tier 1 cloud instance**, with exposed
  topology, affinity, cgroup and calibration receipts in the artifact.

| Pair | Order | IC online ms | Same-point rho online ms | Rho/IC |
| ---: | :--- | ---: | ---: | ---: |
| 1 | IC/rho | 21.126 | 2,390.967 | 113.178× |
| 2 | rho/IC | 38.568 | 2,791.445 | 72.378× |
| 3 | IC/rho | 22.150 | 890.915 | 40.222× |
| 4 | rho/IC | 8.417 | 1,252.612 | 148.824× |
| 5 | IC/rho | 14.334 | 1,619.910 | 113.011× |

The online interval starts with target-dependent work after the reusable
base, index and factor-base logs are ready. It includes target query,
four-summand decomposition, relation check, descent and recovery check.
Rho starts with its first target-dependent walk and stops after independently
verified recovery. Fixture construction and process launch are outside both
intervals. The IC base/index/log setup took 599–628 seconds per fresh
process, and the recorded IC peak RSS was about 10.05 GiB. These costs and
all exclusive phase readings remain in the run rows and raw artifact.

## Receipt and replay

The frozen producer commit is `d0e4aa13074c9fea910592c6e26b60b5e97fc439`;
the IC source SHA-256 is
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
GitHub artifact `11942380164` has SHA-256
`536ed537b302bc16efca9f8b86ae26195d955bd22346da67f6ae316a150ab5c5`.
Its independently generated `analysis.json` has SHA-256
`a44849bef5f772cfb96b0a31f5227ba872a649b33f9351a33276d203ecd590a6`.
The cloud verifier replayed each full base, four-point relation, recovered
scalar, same-point rho answer, phase sum and Tier 1 gate; all five rows are
`PASS_TIER1_PAIR`. The sealer checked the raw ZIP, source, selected candidate,
fixture and analysis digests before writing
[workloads.json](workloads.json) and [measurements.json](measurements.json).
