# P-256 Dickson residual-depth scaling round 6: result

Date run: 2026-10-04

**PASS on all three frozen hypotheses.**  The residual-depth law survived four
original depths, both deterministic terminal fibres, and matched positive and
negative targets.  Every one of 114,240 component systems completed and agreed
with exhaustive signed-point addition.  Residual depths 1 and 2 stayed at
maximum solving degree 3 in all 32 cells; residual depth 3 reached degree 4 in
all 16 controls.

The result strengthens round 5's degree reduction from a single depth into a
measured scaling law through depth 8.  It also makes the remaining obstruction
explicit: total work grows by almost exactly four whenever the original depth
grows by one, because the cost per component is essentially constant while the
ordered branch-pair count quadruples.

## Boundary summary

The algebraic boundary is degree 4.  The enumeration boundary is exactly
`4^(D-r)` component pairs per target.  Rows below aggregate both terminals and
both target classes, so each contains four frozen cells.

| D | residual r | charged components | positive components | max degree | degree / 4 | max columns | total field ops | ops / component | classification |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 5 | 1 | 1,024 | 4 | 3 | 0.75 | 24 | 807,800 | 788.867 | degree advance |
| 5 | 2 | 256 | 4 | 3 | 0.75 | 44 | 1,133,326 | 4,427.055 | degree advance |
| 5 | 3 | 64 | 2 | 4 | 1.00 | 196 | 6,602,032 | 103,156.750 | scaling control |
| 6 | 1 | 4,096 | 4 | 3 | 0.75 | 24 | 3,230,804 | 788.771 | degree advance |
| 6 | 2 | 1,024 | 4 | 3 | 0.75 | 44 | 4,499,895 | 4,394.429 | degree advance |
| 6 | 3 | 256 | 4 | 4 | 1.00 | 196 | 26,342,038 | 102,898.586 | scaling control |
| 7 | 1 | 16,384 | 6 | 3 | 0.75 | 24 | 12,925,397 | 788.904 | degree advance |
| 7 | 2 | 4,096 | 6 | 3 | 0.75 | 44 | 17,972,479 | 4,387.812 | degree advance |
| 7 | 3 | 1,024 | 6 | 4 | 1.00 | 196 | 105,319,482 | 102,851.057 | scaling control |
| 8 | 1 | 65,536 | 19 | 3 | 0.75 | 24 | 51,702,629 | 788.920 | degree advance |
| 8 | 2 | 16,384 | 19 | 3 | 0.75 | 44 | 71,883,372 | 4,387.413 | degree advance |
| 8 | 3 | 4,096 | 19 | 4 | 1.00 | 196 | 421,257,286 | 102,846.017 | scaling control |

At depth 8, residual depth 1 costs 0.7193 times residual depth 2 and
0.1227 times residual depth 3 despite enumerating respectively four and sixteen
times as many components.  The additional branching wins on these toy systems
because a one-level residual component is about 5.56 times cheaper than a
two-level component and about 130 times cheaper than a three-level component.
This is a solver-stage engineering result; the branch exponent remains four per
added chain level.

## Frozen inputs and correctness

The deterministic scan selected terminals 1 and 272 over `F_7681`.

| terminal | depth-8 roots | liftable x | signed points |
|---:|---:|---:|---:|
| 1 | 256 | 114 | 228 |
| 272 | 256 | 128 | 256 |

The toy curve has 7,520 affine points and uses
`b_P256 mod 7681 = 3506`.  All 48 cells completed with zero timeouts, zero
degree-cap stops, zero classification mismatches, and zero exceptional records.
The 114,240 component records charge 723,676,540 finite-field operations in
total.

## Mechanical P-256-depth projection

Holding the depth-8 average component cost fixed and substituting original
depth 18 gives the following arithmetic projection for one two-summand target:

| residual r | ordered branch pairs | depth-8 ops / pair | projected field ops |
|---:|---:|---:|---:|
| 1 | `4^17 = 17,179,869,184` | 788.920 | `1.355e13` |
| 2 | `4^16 = 4,294,967,296` | 4,387.413 | `1.884e13` |
| 3 | `4^15 = 1,073,741,824` | 102,846.017 | `1.104e14` |

This is an extrapolation of a stage counter, not a P-256 measurement.  It omits
P-256 coefficient effects, longer-relation construction, relation collection,
matrix work, and logarithm recovery.  It establishes neither a practical
P-256 decomposition oracle nor an ECDLP speedup.  It does show that residual
depth 1, not 2, is the lowest measured fully charged branch strategy, while the
unpruned branch frontier remains the dominant scaling problem.

## Complete cell table

`correct / complete` counts every ordered component pair in the cell.  Each
digest covers the canonical ordered component record stream, including the
expected sign, solver verdict, completion flags, degree, columns, and field
operations.

| terminal | D | r | target | pairs | correct / complete | max degree | degree / 4 | max columns | field ops | ops / pair | stream SHA-256 | class |
|---:|---:|---:|:--|---:|:--:|---:|---:|---:|---:|---:|:--|:--|
| 1 | 5 | 1 | positive | 256 | 256/256 | 3 | 0.75 | 24 | 202028 | 789.172 | `6760a0b5f6670734977a0e725d299ddb8e40b8b9edbbc9937698d044c4ec5ac5` | degree advance |
| 1 | 5 | 1 | negative | 256 | 256/256 | 3 | 0.75 | 24 | 201872 | 788.563 | `98d0ea7706eefe018e622545c88cda8e44a30d34f442c793643f3a5254ba2e1b` | degree advance |
| 1 | 5 | 2 | positive | 64 | 64/64 | 3 | 0.75 | 44 | 292605 | 4571.953 | `11fe3ee3a9e9cf0ba43c0a5571a02b0e6474affed01b422aab3a29ab36c8897a` | degree advance |
| 1 | 5 | 2 | negative | 64 | 64/64 | 3 | 0.75 | 40 | 274072 | 4282.375 | `9b48b740b40791af13020d4036dc8a69abd0bd7c9ee201686b3b024c748798b0` | degree advance |
| 1 | 5 | 3 | positive | 16 | 16/16 | 4 | 1.00 | 196 | 1906927 | 119182.938 | `cc8f8200a597fe843d91896950d88f8d2f9bbb94a130a4fe64f72e68867f8d79` | scaling control |
| 1 | 5 | 3 | negative | 16 | 16/16 | 4 | 1.00 | 186 | 1395205 | 87200.313 | `714d4eef2c2310828614bb45f7762e85b9f550a463358df43f4b658f9b36142e` | scaling control |
| 1 | 6 | 1 | positive | 1024 | 1024/1024 | 3 | 0.75 | 24 | 807862 | 788.928 | `f3dec2282f0acf6db22ec1c4945cafb046975283f3a8bcb52a71328e139f81fd` | degree advance |
| 1 | 6 | 1 | negative | 1024 | 1024/1024 | 3 | 0.75 | 24 | 807698 | 788.768 | `05229a87776a70b4b1965a4c02ac627803104c05d6e5c27c6462196dcfbf0bca` | degree advance |
| 1 | 6 | 2 | positive | 256 | 256/256 | 3 | 0.75 | 44 | 1153347 | 4505.262 | `3fa2ddd7b2be757931981cbbeceecee9cd5691b8add794e16c334161a705b055` | degree advance |
| 1 | 6 | 2 | negative | 256 | 256/256 | 3 | 0.75 | 40 | 1096607 | 4283.621 | `febb77302f87ed9b82f379aa185feb5457389771aae21d93c5dbe5aa1557d60b` | degree advance |
| 1 | 6 | 3 | positive | 64 | 64/64 | 4 | 1.00 | 196 | 7590411 | 118600.172 | `a9a4e970955e4012d166638e9faaa04992b70c7b7948c37fe71ea8f9b990a8cc` | scaling control |
| 1 | 6 | 3 | negative | 64 | 64/64 | 4 | 1.00 | 186 | 5580924 | 87201.938 | `8ce2cac6821674d1b548cf01a21a853eae0f84302138b89b1ea7ff2a3f2bde98` | scaling control |
| 1 | 7 | 1 | positive | 4096 | 4096/4096 | 3 | 0.75 | 24 | 3231527 | 788.947 | `ef9d68822de7c24407ed973da42e0d88af13b5563d80a781ee8f22f128e2bd8e` | degree advance |
| 1 | 7 | 1 | negative | 4096 | 4096/4096 | 3 | 0.75 | 24 | 3231252 | 788.880 | `1e6774c1aae7255d1f3787a96f07eac16cbe314a3121d6850cad724b842f9555` | degree advance |
| 1 | 7 | 2 | positive | 1024 | 1024/1024 | 3 | 0.75 | 44 | 4602517 | 4494.646 | `2f66efe91e6106045cfd89daacaaa985d66fd940e3e9a21a96f4c8d0265f1718` | degree advance |
| 1 | 7 | 2 | negative | 1024 | 1024/1024 | 3 | 0.75 | 40 | 4387070 | 4284.248 | `45c2bd408a74794211970953833a1309a36d3b91e2a7a9f91f4f778c89cd0297` | degree advance |
| 1 | 7 | 3 | positive | 256 | 256/256 | 4 | 1.00 | 196 | 30347431 | 118544.652 | `1df4d4c9caf2faef186dcfb06d55e89945df183688dd7134dc80df3fca81628b` | scaling control |
| 1 | 7 | 3 | negative | 256 | 256/256 | 4 | 1.00 | 186 | 22322885 | 87198.770 | `fbf6c51708e2d273b2ec702e60ef0d585abca274734a4241f8be412d307eaa5e` | scaling control |
| 1 | 8 | 1 | positive | 16384 | 16384/16384 | 3 | 0.75 | 24 | 12925155 | 788.889 | `b5859aafeada25598a817d27879a334aa3d3790d72383a7b8bb4e0bb586b1471` | degree advance |
| 1 | 8 | 1 | negative | 16384 | 16384/16384 | 3 | 0.75 | 24 | 12925653 | 788.919 | `c45b2c2f6ad494d43dc2cf9f15394d8015b16972035afd0bab6ab06de0e22f49` | degree advance |
| 1 | 8 | 2 | positive | 4096 | 4096/4096 | 3 | 0.75 | 40 | 17572701 | 4290.210 | `3e796cbbd55dc03900a2975af8b025d30a399b15c0a13b18395d556b88d29d59` | degree advance |
| 1 | 8 | 2 | negative | 4096 | 4096/4096 | 3 | 0.75 | 44 | 18364805 | 4483.595 | `0a514062cea20a2d7b8188f7d2ea7d9a7505a238570a591537ab1b3c1490357e` | degree advance |
| 1 | 8 | 3 | positive | 1024 | 1024/1024 | 4 | 1.00 | 186 | 89353745 | 87259.517 | `91abe1db13c85e99f5b9a10d96787bcf6ab2aa25b9b1069e53f113c548c8df44` | scaling control |
| 1 | 8 | 3 | negative | 1024 | 1024/1024 | 4 | 1.00 | 196 | 121247315 | 118405.581 | `01cdf192b7f2e9979eac2b25000e9745e01f2b54f559634856a36ab04b47b9cc` | scaling control |
| 272 | 5 | 1 | positive | 256 | 256/256 | 3 | 0.75 | 24 | 202028 | 789.172 | `84f0b5df8dfb3b55ae95deb97e943cf7d348e9e15c37f38816c7b1661e5fe00b` | degree advance |
| 272 | 5 | 1 | negative | 256 | 256/256 | 3 | 0.75 | 24 | 201872 | 788.563 | `dee9ce37a2ecaf71bfc23b184683553eebbc90b70833dbf822951f51dde0a86a` | degree advance |
| 272 | 5 | 2 | positive | 64 | 64/64 | 3 | 0.75 | 40 | 279846 | 4372.594 | `d0fff0a47020d9370dc79d5d44b319820d121009d6f0db00af80dbaadc1112d5` | degree advance |
| 272 | 5 | 2 | negative | 64 | 64/64 | 3 | 0.75 | 44 | 286803 | 4481.297 | `3bf455ac4569d160e7c48590f082b9e810563814584f6a5a7ce272e7dd5f3313` | degree advance |
| 272 | 5 | 3 | positive | 16 | 16/16 | 4 | 1.00 | 186 | 1410324 | 88145.250 | `638667844eb582219e7d2f277ee9287ee0931e6902733238b9683528e2476cfc` | scaling control |
| 272 | 5 | 3 | negative | 16 | 16/16 | 4 | 1.00 | 196 | 1889576 | 118098.500 | `3d1a07c83c18af423741993c92d4ffdde3e24de2de74784bc9709142371edad0` | scaling control |
| 272 | 6 | 1 | positive | 1024 | 1024/1024 | 3 | 0.75 | 24 | 807540 | 788.613 | `f84f6cb68bb9c0bafe45584781d8fa4ff6f537690c409735e12bae48cebfadf5` | degree advance |
| 272 | 6 | 1 | negative | 1024 | 1024/1024 | 3 | 0.75 | 24 | 807704 | 788.773 | `05b4dd9e5339e2c35e67a667fe4a3381f23f14cb08e64d05f8d44d1342d7634c` | degree advance |
| 272 | 6 | 2 | positive | 256 | 256/256 | 3 | 0.75 | 44 | 1153337 | 4505.223 | `c8bf859a451615355fc8c6b19555040b9582ea493b2043b7f19de252ca48c03c` | degree advance |
| 272 | 6 | 2 | negative | 256 | 256/256 | 3 | 0.75 | 40 | 1096604 | 4283.609 | `b4aeee076e105ed638f45a314ebf02895add12dcffc89c68f19603744e93444f` | degree advance |
| 272 | 6 | 3 | positive | 64 | 64/64 | 4 | 1.00 | 196 | 7589952 | 118593.000 | `b178a3a0b11ad70ed08a6ac8bf21bf33be008e1db056b14f382ed27fcd4667c1` | scaling control |
| 272 | 6 | 3 | negative | 64 | 64/64 | 4 | 1.00 | 186 | 5580751 | 87199.234 | `5feac366144e954fcbf70f92a16166dbb8971795d7d9a7693e52f23352ecb294` | scaling control |
| 272 | 7 | 1 | positive | 4096 | 4096/4096 | 3 | 0.75 | 24 | 3231403 | 788.917 | `fadf1e6f8703f0290644f2d83337392c3f3cea0b7f7051995142460eac5d8ce0` | degree advance |
| 272 | 7 | 1 | negative | 4096 | 4096/4096 | 3 | 0.75 | 24 | 3231215 | 788.871 | `90add01b84ed8022de6377a9acd144c8748b5e33a6d2db15bd41d9f55318a94b` | degree advance |
| 272 | 7 | 2 | positive | 1024 | 1024/1024 | 3 | 0.75 | 44 | 4596223 | 4488.499 | `11f42314365bd6f578926d4a409097f4fb4b6893797843ae1a297191df40f2db` | degree advance |
| 272 | 7 | 2 | negative | 1024 | 1024/1024 | 3 | 0.75 | 40 | 4386669 | 4283.856 | `10af6da6e882a99c5aea1ab36b684b0f693e929c6f4e7ed95e9e11d772c536f6` | degree advance |
| 272 | 7 | 3 | positive | 256 | 256/256 | 4 | 1.00 | 196 | 30327241 | 118465.785 | `18b2d57b230733d932b20980ab09dc999f81ce4cd36e13d2a2df3d455d7a5e54` | scaling control |
| 272 | 7 | 3 | negative | 256 | 256/256 | 4 | 1.00 | 186 | 22321925 | 87195.020 | `fa9875a8745df34fb7d014a7c552849c141f09351fa9b14c1212072f3d070d30` | scaling control |
| 272 | 8 | 1 | positive | 16384 | 16384/16384 | 3 | 0.75 | 24 | 12925855 | 788.932 | `7c963e4f1c8c1d944cb8f601d9ff90713fba88b12eec898d1156a4c9853433d9` | degree advance |
| 272 | 8 | 1 | negative | 16384 | 16384/16384 | 3 | 0.75 | 24 | 12925966 | 788.938 | `3443f81f703d573ffe42e29590daf4f1a0f79de9427188f7406016ab50d612f5` | degree advance |
| 272 | 8 | 2 | positive | 4096 | 4096/4096 | 3 | 0.75 | 40 | 17581170 | 4292.278 | `b766bb7f102fbe528432e1f7c025e3643bb4e71795ce35798ecb86e3a316140b` | degree advance |
| 272 | 8 | 2 | negative | 4096 | 4096/4096 | 3 | 0.75 | 44 | 18364696 | 4483.568 | `e93e9775e9005fee32ae773c7e3870e5b34553b93c55dc119265f67be3ef2490` | degree advance |
| 272 | 8 | 3 | positive | 1024 | 1024/1024 | 4 | 1.00 | 186 | 89379987 | 87285.144 | `70db528fd319cf75fa0cf9c816f3e4a7ca220696fe426950d68e12dff4e3f4f0` | scaling control |
| 272 | 8 | 3 | negative | 1024 | 1024/1024 | 4 | 1.00 | 196 | 121276239 | 118433.827 | `e4d5b7d6168034ee6922ef5d92633920a8776e8a196745d68675e7016fe48020` | scaling control |

## Reproduction and artifact custody

```bash
cargo run --release --bin p256_dickson_residual_scaling -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json
```

The canonical `degree-result.json` is 129,547 bytes with SHA-256
`71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4`.
An immediate independent replay produced the identical byte hash.

The first format included accumulated solver milliseconds.  Its replay kept
every algebraic field identical but naturally changed those timings, so the
whole-file hashes differed.  Both artifacts are retained rather than hidden:

- `degree-result-first-with-wall.json`: 132,133 bytes,
  `13b69b0fe567ab727d387429da467a32663c0d0572870d4db20e51777685f878`;
- `degree-result-replay-with-wall.json`: 132,125 bytes,
  `09fcb589a8e44b2baee2be357a8cfbe3264a58693f8601f271a6abf08efca800`.

Deleting `total_solver_milliseconds` from those two JSON documents makes them
identical to each other and to the canonical result.  The canonical schema
therefore excludes nondeterministic wall time and retains the counted field
operations.

## Decision

Accept the residual-depth law through depth 8.  Retain residual depth 1 as the
best fully charged degree-3 toy strategy, but classify it as solver-stage
engineering because its branch count still scales as `4^(D-1)`.  The next
credible round is not another terminal sweep: it is a frozen compatibility
filter that must reject branch pairs before F4 while proving zero false
negatives against this complete 114,240-component corpus.
