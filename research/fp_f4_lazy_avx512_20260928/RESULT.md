# Prime F4 lazy AVX-512 experiment

The frozen paired Linux result is pending. This branch adds an explicit `F4_FP_LAZY_AVX512=1` path and a direct wide/current/candidate complete-call comparison on an AVX-512F host. The existing scalar default and opt-in AVX2 path remain available. No AVX-512 gain or 2× complete-call result is claimed before the exactness tests and paired receipt are reviewed.
