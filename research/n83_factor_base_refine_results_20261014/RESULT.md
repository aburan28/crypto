# K3300 with 11-byte root states leads the frozen N83 factor-base refinement

The 1,000-point paired panel selected K3300 with 547,800 distinct usable
factor-base points and 3,300 folded columns. Its mean target-online time was
40.638 ms on this cloud allocation. K3100/root22 was the next fastest at
42.513 ms; the paired runner-up-minus-selected mean was 1.875 ms, with a
target-block bootstrap 95% interval of −1.467 to 5.232 ms. The selected
candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`.
Its exact base digest is
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
The result follows the frozen lowest-mean rule; the close K3100 comparison
motivates a fresh paired confirmation.

## Instance and measurement

The curve is `EC1N83Ckb1a1haacb69ee9f37` over the polynomial-basis
`GF(2^83)` with defining polynomial bitmask
`9671406556917033397649557`, model `y² + xy = x³ + x² + 1`, and prime
subgroup order `8569786107849059`. All arms used rank seed 29, symmetric
S3 traversal, four rank threads, one target thread, fresh factor logs and
index per arm, a 1,200-second child cap, and a 14-GiB observed RSS guard.
The 1,000 common public targets have SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`;
the independently checked scalar sidecar has SHA-256
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
The cloud host exposed four CPUs on an AMD EPYC 9V74 processor and 16 GB of
memory. These online times are a secondary, cloud-instance exploratory
factor-base selection measurement. Each online interval includes target
query, point decomposition, relation check, descent, and recovery check;
setup is recorded separately.

| K / layout | IC1 hash suffix | Usable points | Mean online ms | 95% target-block interval ms | Setup s | Peak GiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3300 / root11 | `0f8ab8bfcff7` | 547,800 | **40.638** | 38.328–43.010 | 724.4 | 10.05 |
| 3100 / root22 | `aa5b1411e6a0` | 514,600 | 42.513 | 40.032–45.069 | 632.1 | 13.56 |
| 3000 / root22 | `66fcddf8918c` | 498,000 | 42.743 | 40.066–45.517 | 591.2 | 13.03 |
| 3596 / root4 | `e2fd524905ce` | 596,936 | 43.418 | 40.662–46.264 | 948.1 | 13.96 |
| 3060 / root22 | `ea2149500b31` | 507,960 | 43.684 | 40.949–46.531 | 634.9 | 13.35 |
| 3020 / root22 | `8c4bc59d910f` | 501,320 | 44.032 | 41.385–46.861 | 601.9 | 13.14 |
| 3220 / root21 | `5e3935b37a0a` | 534,520 | 44.685 | 42.094–47.420 | 741.4 | 13.82 |
| 2980 / root22 | `f44e030a52e4` | 494,680 | 45.182 | 42.473–48.087 | 583.7 | 12.93 |
| 2850 / root22 | `55a7c98718cf` | 473,100 | 47.343 | 44.511–50.297 | 533.4 | 12.26 |
| 2900 / root22 | `e77e660a76b0` | 481,400 | 48.323 | 45.298–51.435 | 549.3 | 12.52 |
| 2800 / root22 | `631aee3b1da7` | 464,800 | 49.341 | 46.341–52.455 | 515.0 | 12.02 |
| 2700 / root22 | `1974695d312a` | 448,200 | 51.546 | 48.347–54.886 | 483.7 | 11.52 |
| 2940 / root22 | `751130834947` | 488,040 | 55.865 | 52.420–59.485 | 635.0 | 12.72 |
| 2600 / root22 | `7ca0dce174a8` | 431,600 | 57.564 | 54.033–61.195 | 439.1 | 11.05 |

The suffix and usable-point count identify each full IC1 candidate in the
adjacent `candidate_K*.json` manifests; `measurements.json.gz` retains the
full candidate, workload and run IDs, every one of the 1,000 online
intervals and phase costs per arm, setup, final rank, memory, and process
status. `records_index.json` binds the compressed and decoded JSON digests.
All 14
arms passed full-base and 1,000-target independent replay. The incumbent
K2940 base reproduced its previously sealed base digest. The matched
K2940/K3300 online-mean ratio on this allocation is 1.375; the K3100/K3300
ratio is 1.046. These ratios compare factor bases within this frozen panel.

## Evidence and next validation

[Cloud run 38591919093](https://github.com/aburan28/crypto/actions/runs/38591919093)
executed source commit `b1a4131558b3121ed562ca830abce972aa1dc580` with
IC source SHA-256
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
Its raw artifact ID is `11921089967`, SHA-256
`bd8fc10e07e640f3c4e3afbfaabb21128e84364d7a12a2f44d330d1dd71720ce`.
The cloud analysis and independent local full-base replay both produced
`analysis.json` SHA-256
`5c9d85e3029f85dabe5d91fd92881edea8bc1e8c9fd91d4eedfd0dcf10e046c2`.
The canonical selected-candidate record and all 14 measurement rows are
sealed beside this report. Five separate points frozen before this panel
feed the alternating same-point Tier 1 IC/rho validation in
[PR #1855](https://github.com/aburan28/crypto/pull/1855).
