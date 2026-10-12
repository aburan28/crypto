# Complete standalone conditional-linear-fiber comparison

168 fixed cells, 37632 observations; all completed and verified: True.

The current reference includes all 28 predecessor methods, including initial restriction and all SIMD treatments. Every cold solve and validation cost is charged. Ratios below use matched repetitions; censored costs remain null.

| Candidate | Split | Variables | Family | Reference | Median ratio | 95% interval | Pass |
|---|---|---:|---|---|---:|---|---|
| fiber_simd | regression | 16 | planted | retained_frontier | 0.433256041455865 | [0.40280877193487347, 0.4470072652278834] | False |
| fiber_simd | regression | 16 | planted | same_policy | 1.7805891902985844 | [1.743705413767634, 1.819334975369458] | True |
| fiber_simd | regression | 16 | cross_planted | retained_frontier | 0.4920072147881469 | [0.45350210569840654, 0.5127152661662131] | False |
| fiber_simd | regression | 16 | cross_planted | same_policy | 1.8132669934218997 | [1.7103664353262686, 1.8898763814777664] | True |
| fiber_simd | regression | 16 | unplanted | retained_frontier | 0.4090250877514354 | [0.3965602009733128, 0.42615384030288883] | False |
| fiber_simd | regression | 16 | unplanted | same_policy | 1.954916609166824 | [1.9245741736979278, 2.008654607462306] | True |
| fiber_simd | regression | 20 | planted | retained_frontier | 0.580806359251663 | [0.5663648038303277, 0.5930952401495589] | False |
| fiber_simd | regression | 20 | planted | same_policy | 2.4692296347299623 | [2.4154574069997796, 2.528997027751] | True |
| fiber_simd | regression | 20 | cross_planted | retained_frontier | 0.4423277755108741 | [0.3744910731153991, 0.4637629200113425] | False |
| fiber_simd | regression | 20 | cross_planted | same_policy | 2.578752341395748 | [2.5218866513234532, 2.6771376950087777] | True |
| fiber_simd | regression | 20 | unplanted | retained_frontier | 0.6467188090760396 | [0.6377469346495148, 0.6634600665806789] | False |
| fiber_simd | regression | 20 | unplanted | same_policy | 2.640408054036181 | [2.557695007035985, 2.704952461413872] | True |
| fiber_simd | regression | 24 | planted | retained_frontier | 0.854444274787231 | [0.7042272880999911, 1.0980976927493111] | False |
| fiber_simd | regression | 24 | planted | same_policy | 2.4310865030000928 | [2.3960083181484793, 2.527191255336584] | True |
| fiber_simd | regression | 24 | cross_planted | retained_frontier | 0.5954280879774547 | [0.5128278482511632, 0.7852516245110972] | False |
| fiber_simd | regression | 24 | cross_planted | same_policy | 2.6684583499102796 | [2.5979718363603137, 2.715509900776838] | True |
| fiber_simd | regression | 24 | unplanted | retained_frontier | 1.1018129936702312 | [1.0001564939561502, 1.2361652443608389] | False |
| fiber_simd | regression | 24 | unplanted | same_policy | 2.7013441571033336 | [2.624289790675006, 2.7656368734098224] | True |
| fiber_simd | holdout | 16 | planted | retained_frontier | 0.4419650349650349 | [0.4036697247706422, 0.4883720930232558] | False |
| fiber_simd | holdout | 16 | planted | same_policy | 1.7919832703242236 | [1.7106415094339622, 1.8741224489795918] | True |
| fiber_simd | holdout | 16 | cross_planted | retained_frontier | 0.41524329706667806 | [0.3047358170236462, 0.544981904109773] | False |
| fiber_simd | holdout | 16 | cross_planted | same_policy | 1.8405678625386648 | [1.643835616438356, 2.019370273190307] | True |
| fiber_simd | holdout | 16 | unplanted | retained_frontier | 0.4409631370477951 | [0.37592906998679493, 0.5491277226455586] | False |
| fiber_simd | holdout | 16 | unplanted | same_policy | 2.0361212275586853 | [1.9568621192566693, 2.172693296779913] | True |
| fiber_simd | holdout | 20 | planted | retained_frontier | 0.4877585392000466 | [0.22433386679806244, 0.7839146030511258] | False |
| fiber_simd | holdout | 20 | planted | same_policy | 2.5482309086874846 | [2.1798969957081544, 2.8119857651245552] | True |
| fiber_simd | holdout | 20 | cross_planted | retained_frontier | 0.5080353229219867 | [0.3795966583664069, 0.6957023843069959] | False |
| fiber_simd | holdout | 20 | cross_planted | same_policy | 2.7843586777183895 | [2.239718867651933, 3.487571984707718] | True |
| fiber_simd | holdout | 20 | unplanted | retained_frontier | 0.40396194637328575 | [0.2397391637897967, 0.5683070649178883] | False |
| fiber_simd | holdout | 20 | unplanted | same_policy | 2.750972998473917 | [2.4473607876924053, 2.944261205287515] | True |
| fiber_simd | holdout | 24 | planted | retained_frontier | 0.5411125737193314 | [0.5249885222140228, 0.5890848126232742] | False |
| fiber_simd | holdout | 24 | planted | same_policy | 2.5507872863604577 | [2.325826913456728, 2.7667302246117953] | True |
| fiber_simd | holdout | 24 | cross_planted | retained_frontier | 0.9114295514722195 | [0.804432506941, 1.0246855812290547] | False |
| fiber_simd | holdout | 24 | cross_planted | same_policy | 2.2046523666438023 | [1.0246965490688769, 3.404951259291757] | False |
| fiber_simd | holdout | 24 | unplanted | retained_frontier | 1.1274090487210275 | [0.98331150738931, 1.287099768733177] | False |
| fiber_simd | holdout | 24 | unplanted | same_policy | 2.942116068341957 | [2.7316299869163654, 3.212671150879358] | True |

All gates: `{"fiber_simd": {"incremental": "REJECTED", "retained_frontier": "REJECTED", "same_policy": "REJECTED"}}`

These are bounded generated Boolean solves. No calibrated-operation, production index-calculus, asymptotic or rho claim is established.
