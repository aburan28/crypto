# N83 K3300 root-table layout confirmation

The selected five-byte root11 table is **4.987 ms faster per target** than
the four-byte fixed31 table on 4,000 frozen public points with the **same
K3300 factor base and 903,852,075 exact root keys**. The predeclared paired
root11-minus-root4 95% contiguous-block bootstrap interval is
**[−5.322, −4.662] ms**. Both complete arms passed independent base,
four-point relation and scalar replay. Root11 also used less memory in this
allocation.

| Exact K3300 arm | Mean target online ms | Base/index/log setup s | Root table GiB | Root states GiB | Peak RSS GiB | Status |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| Four-byte fixed31, 11-byte states | 29.768 | 500.733 | 8.00 | 4.63 | 13.05 | `PASS_INDEPENDENT_REPLAY` |
| Five-byte fixed30, 11-byte states | 24.782 | 424.916 | 5.00 | 4.63 | 10.05 | `PASS_INDEPENDENT_REPLAY` |

Both arms used the identical enumerated base hash
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`,
547,800 usable points, 3,300 sign/Frobenius-folded columns, and exact root-key
count. Each ran in a fresh process with an empty base/index/log cache and a
14 GiB process RSS guard. Target-independent preparation is recorded
separately from the online interval. The 4,000-target comparison is a
secondary `cloud_instance_exploratory` layout diagnostic; the selected
five-byte candidate's fresh same-point Tier 1 rho result is recorded in its
own one-target run rows.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. The canonical candidates differ in table construction:

- Four-byte fixed31: `IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h158d17ece5f1`.
- Five-byte fixed30: `IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`.

The public point file SHA-256 is
`451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c`
and the answer sidecar SHA-256 is
`8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101`.
The workload ID is `b3757b0c7e01`. The frozen contrast was
`mean(root11_ms − root4_ms)`, with 80 contiguous 50-target blocks and 10,000
bootstrap resamples using seed `2026101599`. The observed mean was
`−4.986641173` ms and the upper interval endpoint is below zero.

## Receipt and replay

The executed source commit is `009a3025fe362c85e46170e30648edaaec11c6b3`;
the IC source SHA-256 is
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
[GitHub Actions run 38666101824](https://github.com/aburan28/crypto/actions/runs/38666101824)
completed on Linux x86-64 with an exposed AMD EPYC 9V45. Raw artifact
`11942743006` has SHA-256
`9abb1d09c3e52440bf0b39aeec04bf2653e6cbf998c1b2830f2d6868832db2cb`;
its independent `analysis.json` has SHA-256
`e3906ebcd2bc55ffcf613f114299433ee808570f65d4e8a4c2b826b442a432cf`.
The independent verifier checked both full bases, exact root-key parity,
all 8,000 recovered target scalars and four-point relations. The sealer
checked the raw ZIP, source, workload and analysis digests and wrote
[confirmation.json](confirmation.json), the canonical candidate manifests,
and lossless compressed [workload](workload.json.gz) and
[measurement rows](measurements.json.gz). [records_index.json](records_index.json)
gives the compressed and decoded SHA-256 values.
