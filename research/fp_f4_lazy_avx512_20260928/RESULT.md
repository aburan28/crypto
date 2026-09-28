# Prime F4 lazy AVX-512 experiment

The frozen paired Linux result is pending. This branch adds an explicit `F4_FP_LAZY_AVX512=1` path and a direct wide/current/candidate complete-call comparison on an AVX-512F host. The existing scalar default and opt-in AVX2 path remain available. No AVX-512 gain or 2× complete-call result is claimed before the exactness tests and paired receipt are reviewed.

[CI run 36469467083, attempt 1](https://github.com/aburan28/crypto/actions/runs/36469467083/attempts/1) passed the AVX-512 boundary and prime F4 test step and built both binaries, but its AMD EPYC 7763 runner lacked AVX-512F. The benchmark recorded `unsupported_host` and zero process calls. The complete [zero-call receipt](runs/36469467083-attempt1-unsupported.json.gz) has compressed SHA-256 `85fd394c082154fff951848dd23483754bb0f7764e05c330850bfc5c70bfec28` and uncompressed SHA-256 `ca205921d1fac42ccbc3001b5709dcc943de6fd2a84797e1582295307de86089`. No timing sample was discarded.

Attempt 2 of the same run was canceled during the exactness-test step, before either binary build or any benchmark process, for the pre-measurement AVX2 dispatch-guard amendment recorded in the protocol. It produced no timing sample.
