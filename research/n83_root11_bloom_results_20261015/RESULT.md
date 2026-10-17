# N83 root11 Bloom prefilter speeds four-summand lookup

An optional three-hash, 512 MiB Bloom prefilter reduced K3300/root11 mean
target online time from **34.766 ms** to **31.324 ms** on 1,000 paired N83
points. The paired Bloom-on minus Bloom-off mean was **−3.441 ms** with a
predeclared 95% contiguous-block bootstrap interval of **[−3.617, −3.266]
ms**. Both arms used the same 547,800-point factor base and 903,852,075-key
root index, completed every target, and passed independent base, scalar, and
four-point relation replay under the 14 GiB guard. The prefilter skips an
exact index lookup after a negative test; a positive test still uses the
exact lookup, preserving the lookup result.

| Arm | Mean target online ms | Setup s | Peak RSS GiB | Bloom bytes | Status |
| :--- | ---: | ---: | ---: | ---: | :--- |
| Off | 34.766 | 595.341 | 10.052 | 0 | `PASS_INDEPENDENT_REPLAY` |
| On | 31.324 | 629.328 | 10.552 | 536,870,912 | `PASS_INDEPENDENT_REPLAY` |

The online saving is **9.90%** relative to the off arm. The target probe
count was identical in both arms (55,787.783 mean); the change acts on the
cost of each indexed lookup. The target-independent setup increase was
33.988 s and the observed peak RSS increase was 0.500 GiB. The target
interval includes query, four-summand lookup, relation check, descent and
recovery check. Each arm prepared its base, index and logs in a fresh process.
The 20 contiguous blocks of 50 targets were bootstrapped 10,000 times with
seed `2026101507`; the frozen promotion gate required the upper interval
endpoint to be below zero.

This was a cloud-instance stage comparison on the previously used 1,000-point
refinement panel, recorded as `cloud_instance_exploratory`. The next promotion
step is a fresh same-point Tier 1 IC/rho run with the Bloom-on candidate. The
curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)` with
modulus bitmask `9671406556917033397649557`, model `y²+xy=x³+x²+1`,
subgroup order `8569786107849059`, and cofactor `1128547018`. The base hash
is `a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
The exact Bloom-on candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0hc3b6db4ff8c0`.

## Receipts

The executed source commit is `a9b44ee1454beb2bbe8614b681714df07fd11f3c`,
with IC source SHA-256
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.
[GitHub Actions run 38674728464](https://github.com/aburan28/crypto/actions/runs/38674728464)
completed the producer, exact-source tests, and independent verifier. Raw
artifact `11945754917` has SHA-256
`098c1ea1383c3ddaf78aaf1ac28bc32c80cb2caa5c69ccfb155bb060f5d2bda2`;
its analysis has SHA-256
`2976f932e8eb2fd2d285b1f40f065c43dd6573fb6879845bd915f9b7ee675a6a`.
The [selection](selection.json), two canonical candidate manifests,
[measurement rows](measurements.json.gz), and [digest index](records_index.json)
seal the exact result. The shared canonical 1,000-point
[workload](../n83_factor_base_refine_results_20261014/workload.json.gz) has ID
`0d7b05d9fca6` and public-point SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`.
