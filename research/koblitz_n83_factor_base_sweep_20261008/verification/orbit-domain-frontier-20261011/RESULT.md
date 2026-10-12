# N83 stored-base size/arity screen (static)

The [reproducible receipt](static-screen.json) covers all 30 stored `a=0` primary factor-base objects in the four committed construction manifests: 27 pilot objects at K=64/256/600 and one each at K=1,182, K=2,001 structured, and K=2,048. The last three designed sizes, K=4,096/8,192/16,627, have no stored object. Each listed object has 166K point records. This screen performed **zero new SAT construction, solver search, factor-base construction or S3 transfer**; the 9,000-second charged local experiment allowance is unchanged.

The table applies the exact `chain_s3_max_variables` and `frobenius_orbit_domain_size` source formulas to the all-affine chained-S3 model for five or six summands. `V` is the projected total variable count, `C` is the exact orbit-domain clause count, and `E` is the number of unordered point multisets divided by the primary subgroup order. `E` is a uniform-target counting mean, **not** a measured probability, fixed-target relation yield, search rate or runtime. The source and manifest SHA-256 identities are in the receipt. The [generator](generate_static_screen.py) asserts the executed K=2,001 five-summand and previously admitted K=2,048 five-summand counts.

| K | Stored primary objects | m=5 V | m=5 C | m=5 E | m=6 V | m=6 C | m=6 E |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 9 | 91,638 | 39,350 | 4.67e-7 | 114,182 | 47,220 | 8.27e-4 |
| 256 | 9 | 93,558 | 121,910 | 4.78e-4 | 116,486 | 146,292 | 3.38 |
| 600 | 9 | 96,998 | 269,830 | 0.0338 | 120,614 | 323,796 | 561 |
| 1,182 | 1 | 102,818 | 520,090 | 1.00 | 127,598 | 624,108 | 3.28e4 |
| 2,001 | 1 | 111,008 | 872,260 | 13.9 | 137,426 | 1,046,712 | 7.72e5 |
| 2,048 | 1 | 111,478 | 892,470 | 15.7 | 137,990 | 1,070,964 | 8.87e5 |
| 4,096 | 0 | 131,958 | 1,773,110 | 501 | 162,566 | 2,127,732 | 5.68e7 |
| 8,192 | 0 | 172,918 | 3,534,390 | 1.60e4 | 211,718 | 4,241,268 | 3.63e9 |
| 16,627 | 0 | 257,268 | 7,161,440 | 5.52e5 | 312,938 | 8,593,728 | 2.54e11 |

Under the *proposed*, still unauthorised K=2,048 gate's 150,000-variable and 1,000,000-domain-clause caps, every stored primary size is statically admissible at m=5. At m=6, K=1,182 and smaller are statically admissible; K=2,001 and K=2,048 exceed the domain-clause cap. A pass here says only that two count thresholds permit attempting construction. It gives no memory, wall time, solver, relation, rank, linear algebra or target-log result. The six-summand K=1,182 counting mean is much larger than the five-summand K=2,048 mean, but its larger SAT core and search difficulty could outweigh that difference.

The outstanding empirical step remains the separately proposed guarded K=2,048 five-summand capacity gate in [the admission receipt](../k2048-orbit-domain-admission-20261010/RESULT.md), conditional on an additional 900-second budget. Once construction admits a candidate, compare same-target cold search and ultimately complete precompute, relation collection, independent-rank growth, sparse linear algebra and individual logarithm against matched rho. **No factor base can yet be selected by minimum complete IC runtime.**

Validation: the generator reproduced the receipt from the committed manifests. The study Python suite passed 29/29, the boundary suite 39/39, and the exporter example 19/19. The required unfiltered release library command was run with two threads and reported 1,261 passing, 17 ignored and no failures before interruption at two unrelated long JV cases. Excluding exactly those two cases, the broad library run completed with 3,829 passing, 110 ignored and nine failures; all nine were Pollard collaboration tests blocked from opening local sockets by the sandbox (`Operation not permitted`). The exact Pollard collaboration group passed 37/37 with loopback access. The [validation receipt](validation.json) hashes the three compressed raw Rust logs. Thus these checks provide broad relevant coverage, but do not constitute an unqualified full-suite pass.
