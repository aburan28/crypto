# Stage 185: selected dense-pair default replay

An exact detached checkout of selection commit `8014149a2` passed ten F4 tests and three backend tests in both default-dense and explicit quadratic-control modes. With no selector environment variable, the target replay reported 1,011,275 dense selections, zero quadratic selections, zero full-M4RI matrices, and exact exhaustive UNSAT.

The load-affected default replay took 179.680324 wall seconds, 348.510060 core-seconds, and 3335372800 bytes RSS. Same-binary direct MITM took 2.343044 wall seconds, 0.689407 core-seconds, and 45367296 bytes RSS. The resulting stage ratios are 76.69x wall, 505.52x CPU, and 73.52x RSS; they are not a full-method comparison.

The first exact `--locked` build failed because the selection commit did not track `Cargo.lock`. The successful retry supplied lock SHA-256 `4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5`; restoring that file to the PR is required.

The failed build, successful build, four validation commands, and two replay processes charge 786.160829 wall seconds, 1452.597108 core-seconds, and 3966713856 bytes maximum RSS across 8 components.

The dense pair default is validated on this one target. The result remains implementation engineering, not a full attack or SOTA evidence.
