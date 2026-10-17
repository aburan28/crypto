# Exact support ceilings for the primary cross-stream pilot pairs

The [saved screen](screen.json) covers all **648** ordered distinct-stream
pairs among the 27 primary-arm pilot archives at 64, 256 and 600 signed
Frobenius orbits. It uses the exact pool-minus-base residual orbit counts from
the [independently audited overlap census](../cross-stream-pilot-overlap-20261015/RESULT.md),
which binds every archived object to the existing S3 factor-base manifest.
There are 33 empty and 57 singleton residual pairs; the current worker
requires at least two residual orbits, leaving **558 source-admissible pairs**.
The saved result has SHA-256
`82ff0dcbd4dbb9fd7ee4269cade509abd1f2d1eefc973dcc64900bea3af36e9e`.

For a base with `B` points, an exact disjoint residual with `R` points, `m`
summands and `t` residual summands, the number of unordered point multisets
with repetition is exactly
`C(B+m-t-1,m-t) × C(R+t-1,t)`. The count is an **upper bound** on distinct
subgroup sums: collisions can only reduce that number. A target uniform over
the whole primary subgroup of order `2417851639230796216685689` therefore
hits one of these sums with probability at most `min(1,count/order)`. The
worker's fixed public query coefficients are deterministic, so this sampling
premise is a planning screen, not a measured relation probability.

The following bins count only the 558 pairs admitted by the current worker.
“Vacuous” means the multiset count reaches or exceeds the subgroup order;
it does not establish a likely successful solver call.

| Total / residual summands | Under 1% | 1–10% | 10–50% | 50–100% | Vacuous |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5 / 1 | 414 | 102 | 42 | 0 | 0 |
| 5 / 2 | 432 | 84 | 42 | 0 | 0 |
| 6 / 1 | 62 | 144 | 3 | 0 | 349 |
| 6 / 2 | 23 | 46 | 99 | 30 | 360 |

The retained hosted K600 Gray pair has 97 residual orbits. Its five-summand,
two-residual count is exactly `21349930494938690639600`, giving a **0.883%**
uniform-target hit ceiling. The matched disjoint hash pair has 600 residual
orbits and count `816832271335643493360000`, a **33.783%** ceiling. Their
5,000-conflict searches were both censored with zero models; those solver
observations neither attain nor estimate these ceilings. The six-summand
two-residual count is large enough to make this ceiling vacuous for 360
admissible pairs, which supports testing solver access at that shape after
the fixed-witness gate. It does not show that six-summand solving is faster.

The [generator](screen.py) pins the overlap, independent-audit and
prior-comparison receipt SHA-256 values; reconstructs the full 648-pair
universe; checks every residual partition; and recomputes all four worker
shapes. The [focused tests](test_screen.py) independently enumerate small
multisets, verify the disjoint-pool partition identity, check the hosted pair
against its frozen plan, and reject a changed residual count. Reproduce with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-support-20261017/screen.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-support-20261017 -p test_screen.py -v
```

This stage read committed files only, launched no solver, and made no S3
write. It does not measure natural partial yield, large-prime graph closure,
independent rank, linear algebra, individual-log extraction or complete IC
runtime. The wider v4 grid still contains 26,939,278,080 declared addresses;
this screen covers the 648 exact archived primary cross-stream pairs and four
currently shaped mixed-orbit worker settings, with other addresses retaining
their source dispositions. A lowest-total-runtime base remains unselected.
