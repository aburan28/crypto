# Stage 174: size-gated contiguous rows in native Boolean F4

This stage keeps the single-target work on the repository-native Boolean F4 backend. It adds an opt-in contiguous row arena for large M4RI matrices while retaining the segmented representation below an explicit packed-word threshold. The exact public target remains the already-opened `n=59, ell=9, m=3` true-negative, so this is same-target engineering rather than a fresh holdout.

The selected setting is `PQ_F4_FLAT_M4RI=1` with `PQ_F4_FLAT_M4RI_MIN_WORDS=1048576`, twelve Rayon workers, and one complete 512-value X1 batch. It routes 241 of 723 M4RI matrices through the contiguous arena. Across three interleaved exact-commit pairs, the median candidate/control ratios are `0.901566` wall, `0.868772` CPU, and `1.009054` peak RSS. The representative selected run takes `26.358218` wall seconds, `143.569806` core-seconds, and `2676736000` bytes peak RSS.

Every candidate and control visits all 512 masks, completes all 242 rational fixed-X1 systems, performs exactly 99,199,976,264 charged F4 word XORs, preserves equation fingerprint `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`, finds no algebraic roots, and returns exhaustive UNSAT. The same flat output path also returns an exact polynomial-valid and curve-group-valid witness on the earlier SAT fixture. Thirteen Boolean-F4 tests and three backend tests pass.

The all-flat policy is rejected for single-core UNSAT. The 8 MiB hybrid is neutral/slightly slower there, so single-core retains segmented rows. Variable-order search, native/LTO code generation, and PGO are also rejected and charged in the development ledger.

The exact clean build costs `204.886368` wall seconds and `190.145197` core-seconds. Direct MITM remains decisively faster: the selected F4 run is `8.07x` slower by wall, `44.38x` more expensive by CPU, and `59.60x` larger by RSS than the same-binary direct median.

The seven-gate status is unchanged in substance: licensed Magma is missing, full cost does not beat automorphism-optimized rho, and unaffiliated reproduction and novelty review are missing. This is not Koblitz index-calculus SOTA.
