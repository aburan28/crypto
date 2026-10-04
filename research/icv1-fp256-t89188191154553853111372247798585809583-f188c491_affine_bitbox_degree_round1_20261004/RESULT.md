# P-256 affine-bitbox degree round 1: result

Date run: 2026-10-04

The affine-bitbox factor base is **not** a solving-degree improvement over the
accepted Dickson factor base.  On every completed paired `S3` cell its solving
degree was higher: 5 instead of 4 at depths 3 and 4, and 6 instead of 4 at
depth 5.  The hypothesis frozen in `PROTOCOL.md` is falsified.

The useful result is a formulation result, not a new-factor-base result.  The
quadratic incidence lift lowers the Dickson solving degree from 4 to 3 at
depths 2 and 3, then reaches 4 at depths 4 and 5.  It is the next branch worth
pursuing with the Dickson trace structure.  This remains a small-prime solver
diagnostic; no P-256 Gröbner degree was measured.

## Degree table

All 72 systems completed within the 30-second per-system budget, all were
certified, and all certified verdicts agreed with exhaustive enumeration.
There were no timeouts, pairs above the degree bound, or correctness
mismatches.  `cols` and `ops` are medians to the last productive F4 step.

| depth | formulation | family | verdict | targets | correct | complete | solving degree min/median/max | cols | field ops | degree / matched Dickson | classification |
|---:|:--|:--|:--|---:|:--:|:--:|:--|---:|---:|---:|:--|
| 2 | S3 | Dickson | negative | 3 | yes | yes | 4/4/4 | 38 | 3,703 | 1.00 | reference |
| 2 | S3 | affine-bitbox | negative | 3 | yes | yes | 5/5/5 | 135 | 8,409 | 1.25 | regression |
| 2 | quadratic incidence | Dickson | negative | 3 | yes | yes | 3/3/3 | 112 | 23,736 | 1.00 | reference |
| 2 | quadratic incidence | affine-bitbox | negative | 3 | yes | yes | 3/3/3 | 218 | 108,678 | 1.00 | regression in secondary costs |
| 3 | S3 | Dickson | positive | 3 | yes | yes | 4/4/4 | 186 | 90,360 | 1.00 | reference |
| 3 | S3 | affine-bitbox | positive | 3 | yes | yes | 5/5/5 | 350 | 356,996 | 1.25 | regression |
| 3 | S3 | Dickson | negative | 3 | yes | yes | 4/4/4 | 169 | 67,442 | 1.00 | reference |
| 3 | S3 | affine-bitbox | negative | 3 | yes | yes | 5/5/5 | 350 | 179,554 | 1.25 | regression |
| 3 | quadratic incidence | Dickson | positive | 3 | yes | yes | 3/3/3 | 189 | 286,851 | 1.00 | reference |
| 3 | quadratic incidence | affine-bitbox | positive | 3 | yes | yes | 3/3/3 | 370 | 4,010,913 | 1.00 | regression in secondary costs |
| 3 | quadratic incidence | Dickson | negative | 3 | yes | yes | 3/3/3 | 196 | 203,702 | 1.00 | reference |
| 3 | quadratic incidence | affine-bitbox | negative | 3 | yes | yes | 3/3/3 | 370 | 734,136 | 1.00 | regression in secondary costs |
| 4 | S3 | Dickson | positive | 3 | yes | yes | 4/4/4 | 309 | 1,162,960 | 1.00 | reference |
| 4 | S3 | affine-bitbox | positive | 3 | yes | yes | 5/5/5 | 971 | 24,587,681 | 1.25 | regression |
| 4 | S3 | Dickson | negative | 3 | yes | yes | 4/4/4 | 309 | 767,028 | 1.00 | reference |
| 4 | S3 | affine-bitbox | negative | 3 | yes | yes | 5/5/5 | 999 | 5,259,134 | 1.25 | regression |
| 4 | quadratic incidence | Dickson | positive | 3 | yes | yes | 4/4/4 | 952 | 3,743,183 | 1.00 | reference |
| 4 | quadratic incidence | affine-bitbox | positive | 3 | yes | yes | 4/4/4 | 1,690 | 218,589,170 | 1.00 | regression in secondary costs |
| 4 | quadratic incidence | Dickson | negative | 3 | yes | yes | 4/4/4 | 964 | 3,969,470 | 1.00 | reference |
| 4 | quadratic incidence | affine-bitbox | negative | 3 | yes | yes | 4/4/4 | 1,690 | 8,689,822 | 1.00 | regression in secondary costs |
| 5 | S3 | Dickson | positive | 3 | yes | yes | 4/4/4 | 648 | 30,719,829 | 1.00 | reference |
| 5 | S3 | affine-bitbox | positive | 3 | yes | yes | 6/6/6 | 5,117 | 1,767,777,263 | 1.50 | regression |
| 5 | quadratic incidence | Dickson | positive | 3 | yes | yes | 4/4/4 | 1,305 | 55,302,796 | 1.00 | reference |
| 5 | quadratic incidence | affine-bitbox | positive | 3 | yes | yes | 4/4/4 | 2,989 | 9,690,106,874 | 1.00 | regression in secondary costs |

Depth 2 had no common positive target.  Depth 5 had no common negative target;
the frozen selection rule therefore emitted only the available verdict class
in those two depths.

Against the ordinary Dickson `S3` boundary specifically, Dickson quadratic
incidence has degree ratios `0.75, 0.75, 1.00, 1.00` at depths 2 through 5.
That is a degree advance at the first two depths, not a demonstrated
asymptotic reduction: the advantage disappears at depths 4 and 5.

## Direct P-256 factor base

The candidate does improve cardinality slightly.  The frozen eight-window
selection chose offset `0x180000`:

| FB1 | family | bit width | columns | signed points | m=17 Poisson success | versus Dickson columns |
|:--|:--|---:|---:|---:|---:|---:|
| `FB1ha8fa8dae90f9` | affine-bitbox | 18 | 131,440 | 262,880 | 96.372082631861% | +201 (+0.1532%) |
| `FB1h2ea06bef7f7a` | Dickson-torus | 18 | 131,239 | 262,478 | 96.049526975874% | reference |

The success-model increase is 0.32256 percentage points.  It does not offset
the observed algebraic regression.  The selected FB1 has full digest
`a8fa8dae90f9b2293103476492b39685842f0e049294b11a23cde29fa55a321d`
and sorted point-key digest
`a0c0a105ee4096bc0a1181839ff8e1023fff207bd290009598e9d99b66c84147`.
The native rebuild verified every emitted point row and all eight frozen
candidate counts.

## Reproduction and evidence

```bash
cargo run --release --bin p256_factor_base_degree -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_bitbox_degree_round1_20261004/degree-result.json

cargo run --release --bin p256_bitbox_factor_base -- \
  --curve icv1-fp256-t89188191154553853111372247798585809583-f188c491 \
  --verify \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_bitbox_degree_round1_20261004/factor-base-result.json
```

| artifact | bytes | SHA-256 |
|:--|---:|:--|
| `degree-result.json` | 53,054 | `efe17e77cd23affe03b7b2ee8146490409afbd779f694058606caf798909c3be` |
| `factor-base-result.json` | 1,911 | `bf4c0f95ca46bf92b236eea7b8f37a3902f744b07a898a41f6b120faab2e35a4` |

The next bounded iteration should retain the Dickson chain and sweep
nondegenerate terminal trace constants/cosets, first on a toy prime with a
larger torus quotient.  Its success condition should be a completed common
positive solving degree below 4 at depths 4 and 5; cardinality and matrix cost
must remain separate gates.
